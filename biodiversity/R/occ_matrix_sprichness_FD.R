################################################################################
##
## Script to determine the list of species in each surveys of RLS data  
## (occurrence matrix), taxonomic entropy, functional diversity, functional 
## entropy, and trophic and size structures
##
## occ_matrix_sprichness_FD.R
##
## 19/10/2022
##
## Ulysse Flandrin, adapted from Sebastien Villéger
##
################################################################################

## cleaning memory
rm(list=ls())

##------------------- loading datasets-------------------
load(here::here("data","data_species.Rdata"))
load(here::here("data","data_surveys.Rdata"))

#phylogeny of ray-finned fishes
trees <- ape::read.tree(file = here::here("biodiversity", "data","Code&Data_Siquiera2020",
                                          "TACT", "Reef_fish_all_combined.trees")) 
tree <- trees[[1]]; rm(trees)

# -------------------biomass of species in each survey as a classical matrix -------------------
surveys_sp_biom <- data_surveys %>% 
  select(SurveyID, species, biomass) %>%
  group_by(SurveyID, species) %>%
  summarize( sp_biom=sum(biomass) ) %>%
  pivot_wider(names_from = species, values_from = sp_biom, values_fill = 0) %>%
  column_to_rownames(var="SurveyID") %>%
  as.matrix()

dim(surveys_sp_biom) # OK

#Change species names:  Abudefduf_luridus should be Similiparma_lurida and
#                       Rhinesomus_triqueter should be Lactophrys_triqueter
colnames(surveys_sp_biom)[which(colnames(surveys_sp_biom) == "Abudefduf_luridus")] <- "Similiparma_lurida"
colnames(surveys_sp_biom)[which(colnames(surveys_sp_biom) == "Rhinesomus_triqueter")] <- "Lactophrys_triqueter"

# occurrence of species in surveys ----
surveys_sp_occ <- surveys_sp_biom
surveys_sp_occ[surveys_sp_occ!=0] <- 1

#relative biomass of species in surveys
surveys_sp_pbiom<-surveys_sp_biom/apply(surveys_sp_biom,1,sum)

save(surveys_sp_occ, file=here::here("biodiversity", "outputs", "occurrence_matrix_sp_survey_01.Rdata"))
save(surveys_sp_pbiom, file=here::here("biodiversity", "outputs", "occurrence_matrix_sp_survey_relative_biomass.Rdata"))


# -------------------relative biomass of families in each survey as a classical matrix -------------------
data_surveys_fam <- data_species %>% 
  dplyr::select(species, family) %>%
  dplyr::right_join(data_surveys)

surveys_family_biom <- data_surveys_fam %>% 
  select(SurveyID, family, biomass) %>%
  group_by(SurveyID, family) %>%
  summarize( family_biom = sum(biomass) ) %>%
  pivot_wider(names_from = family, values_from = family_biom, values_fill = 0) %>%
  column_to_rownames(var="SurveyID") %>%
  as.matrix()

dim(surveys_family_biom) # OK

#relative biomass of families in surveys
surveys_family_pbiom<-surveys_family_biom/apply(surveys_family_biom,1,sum)

save(surveys_family_pbiom, file=here::here("biodiversity", "outputs", "occurrence_matrix_family_survey_relative_biomass.Rdata"))


#-------------------same process but only with species known in our phylogeny-------------------
surveys_sp_biom <- surveys_sp_biom[, colnames(surveys_sp_biom) %in% tree[["tip.label"]] ]
dim(surveys_sp_biom) # 10 species missing in the phylogeny:
#[1] "Pareques_acuminatus"   "Odontoscion_dentex"    "Equetus_punctatus"     "Mugil_cephalus"        "Crenimugil_crenilabis"
#[6] "Bothus_mancus"         "Odontoscion_xanthops"  "Pareques_viola"        "Mugil_galapagensis"    "Neomyxus_leuciscus" 

surveys_sp_occ <- surveys_sp_biom
surveys_sp_occ[surveys_sp_occ!=0] <- 1

surveys_sp_pbiom<-surveys_sp_biom/apply(surveys_sp_biom,1,sum)

save(surveys_sp_occ, file=here::here("biodiversity", "outputs", "occurrence_matrix_spinphylogeny_survey_01.Rdata"))
save(surveys_sp_pbiom, file=here::here("biodiversity", "outputs", "occurrence_matrix_spinphylogeny_survey_relative_biomass.Rdata"))


### -------------------taxonomic diversity in surveys -------------------####
#Select species in the phylogeny
data_species <- data_species %>% filter(species_corrected %in% tree[["tip.label"]] )
data_surveys <- data_surveys %>% filter(species %in% data_species$species ) #loss of 115 rls observation


# taxo richness and Shannon-like entropy (eq number of species, exp(H))
surveys_biodiversity <- data.frame(
  taxo_richness = apply(surveys_sp_occ, 1,  sum),
  taxo_entropy = apply(surveys_sp_pbiom, 1, 
                       function(x) { exp( (-1) * sum( x[x!=0]*log(x[x!=0]) ) ) }
  )
)

##------------------- functional diversity in surveys -------------------####

# dataframe with species as row names and traits as variables
sp_traits<-data_species
row.names(sp_traits)<-sp_traits$species_corrected
sp_traits<-sp_traits[,c("Size", "Diet", "Position", "Activity")]

# type of traits
traits_cat<-data.frame(trait_name=c("Size", "Diet", "Position", "Activity"),
                       trait_type=c("Q", "N", "O", "N")
)

# Gower distance between species ----
sp_gower<-funct.dist(sp_tr = sp_traits, tr_cat = traits_cat, metric = "gower")
summary(as.numeric(sp_gower)) # most distances are small (Q3=0.32) 
# warning message means some species have same trait values
# => not an issue for computation of Chao et al 2019 indices


# computing functional distinctiveness
funct_distinctiveness_sp <- apply( as.matrix(sp_gower), 1, sum) / (ncol(as.matrix(sp_gower))-1) #Grenié et al. 2018

funct_distinct_surveys_raw <- lapply(1:nrow(surveys_sp_occ),function(i){
    Sp <- colnames(surveys_sp_occ) [which(surveys_sp_occ[i,]==1)]
    mean(funct_distinctiveness_sp[Sp])
})

surveys_biodiversity$funct_distinctiveness <-  do.call(rbind,funct_distinct_surveys_raw)[,1]

# computing functional richness and functional entropy ----
# applying Chao 2019 framework with package mFD

# richness on species occurrences with q=0
surveys_biodiversity$funct_richness<-mFD::alpha.fd.hill(asb_sp_w = surveys_sp_occ, 
                                                       sp_dist = sp_gower,
                                                       q=0, 
                                                       tau="mean",
                                                       details_returned =FALSE)[,1]

# richness on species relative biomass with q=1
surveys_biodiversity$funct_entropy<-mFD::alpha.fd.hill(asb_sp_w = surveys_sp_pbiom, 
                                                      sp_dist = sp_gower,
                                                      q=1, 
                                                      tau="mean",
                                                      details_returned =FALSE)[,1]

summary(surveys_biodiversity)




##-------------------biomass belonging to low and high trophic levels -------------------####
summary(data_species$Diet)

# species with low TL = herbivores_microvores_detritivores diet
species_lowTL<- data_species %>%
  filter(Diet=="herbivores_microvores_detritivores") %>%
  pull(species_corrected) 
length(species_lowTL) # 196 species

biom_lowTL <- rowSums(surveys_sp_biom[,species_lowTL]) 
summary(biom_lowTL) # from 0 to 1123316 , Q1=3643, median=21833, Q3=24280

# species with medium TL = all type of invertivorous diet (including planktivores and corrallivores)
species_mediumTL<- data_species %>%
  filter(Diet %in% c("corallivores", "sessile_invertivores", 
                     "microinvertivores", "macroinvertivores", "crustacivores", 
                     "planktivores")  ) %>%
  pull(species_corrected) 
length(species_mediumTL) # 743 species

biom_mediumTL <- rowSums(surveys_sp_biom[,species_mediumTL]) 
summary(biom_mediumTL) # from 0 to 5629533, Q1=5653, median=13216, Q3=28245


# species with high TL = piscivores diet
species_highTL<- data_species %>%
  filter(Diet=="piscivores") %>%
  pull(species_corrected) 
length(species_highTL) # 75 species

biom_highTL <- rowSums(surveys_sp_biom[,species_highTL])
summary(biom_highTL) # from 0 to 474937.1, Q1=7.7, median=685.6, Q3=2511.0

# merging
surveys_biomTL<-data.frame(biom_lowTL, biom_mediumTL, biom_highTL) %>%
  rownames_to_column("SurveyID")
summary(surveys_biomTL)


# size_structure of fishes (weighted by number of individuals)
surveys_size<- data_surveys %>%
  select(SurveyID, species, size_class, number) %>%
  uncount( number ) %>%
  select( - species) %>%
  group_by(SurveyID) %>%
  summarize( size_Q1=quantile(size_class,0.25) ,
             size_median=quantile(size_class,0.5) ,
             size_Q3=quantile(size_class, 0.75),
             size_p90=quantile(size_class, 0.90)
  )
summary(surveys_size)  

# merging 
surveys_biodiversity <- surveys_biodiversity %>%
  rownames_to_column("SurveyID")  %>%
  left_join(surveys_size) %>%
  left_join(surveys_biomTL)
summary(surveys_biodiversity)

# filtering 0.1% outliers
surveys_biodiversity_without_outliers <- surveys_biodiversity %>%
  filter(taxo_richness < quantile(surveys_biodiversity$taxo_richness,0.999)) %>% 
  filter(funct_distinctiveness < quantile(surveys_biodiversity$funct_distinctiveness,0.999)) %>% 
  #filter(taxo_entropy < quantile(surveys_biodiversity$taxo_entropy,0.999)) %>% 
  #filter(funct_richness < quantile(surveys_biodiversity$funct_richness,0.999)) %>% 
  filter(biom_lowTL < quantile(surveys_biodiversity$biom_lowTL ,0.999)) %>%
  filter(biom_mediumTL < quantile(surveys_biodiversity$biom_mediumTL ,0.999)) %>%
  filter(biom_highTL < quantile(surveys_biodiversity$biom_highTL ,0.999))

## saving as Rdata ##
save(surveys_biodiversity_without_outliers, file=here::here("biodiversity","outputs", "surveys_biodiversity.Rdata") )
