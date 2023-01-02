################################################################################
##
## Assess the total biomass of fishery species in each transect
##
## fishery_species_biomass.R
##
## 24/11/2022
##
## Ulysse Flandrin
##
################################################################################
#-----------------Loading data---------------------
load(here::here("data", "data_surveys.Rdata"))
load(here::here("data", "data_species.Rdata"))
load(here::here("recycling", "outputs","flux_final_data_surveys.Rdata"))

#----------------- fishery species ---------------------
data_surveys_fishery <- data_surveys %>%
  left_join(data_species[,c("species", "family", "Size")]) %>%
  rename(max_size = Size)
table(data_surveys_fishery$family)

#all species targeted in the family:
all_sp <- c("Acanthuridae", "Caesionidae", "Carangidae", "Ephippidae", "Haemulidae", "Kyphosidae",
            "Labridae", "Lethrinidae", "Lutjanidae", "Mullidae", "Nemipteridae", "Scaridae",
            "Scombridae", "Serranidae", "Siganidae", "Sparidae", "Sphyraenidae")
#families with species targeted if larger than 20cm:
target_larger_20cm <- c("Balistidae", "Holocentridae", "Pomacanthidae", "Priacanthidae")

#non-targeted families:
setdiff(unique(data_surveys_fishery$family), c(all_sp, target_larger_20cm))
# [1] "Pomacentridae"  "Tetraodontidae" "Chaetodontidae" "Scorpaenidae"   "Cirrhitidae"    "Ostraciidae"    "Zanclidae"     
# [8] "Monacanthidae"  "Pempheridae"    "Fistulariidae"  "Sciaenidae"     "Mugilidae"      "Bothidae"  

#keep only targeted fishes
data_surveys_fishery <- data_surveys_fishery %>% 
  filter(family %in% all_sp | (family %in% target_larger_20cm & max_size > 20))

# -------------------biomass of fishery species in each survey -------------------
surveys_fishery_sp_biom <- data_surveys_fishery %>% 
  select(SurveyID, species, biomass) %>%
  group_by(SurveyID, species) %>%
  summarize( fishery_biomass=sum(biomass) ) %>%
  pivot_wider(names_from = species, values_from = fishery_biomass, values_fill = 0) %>%
  column_to_rownames(var="SurveyID") 


surveys_fishery_biom <- as.data.frame(rowSums(surveys_fishery_sp_biom))
surveys_fishery_biom <- tibble::rownames_to_column(surveys_fishery_biom)
colnames(surveys_fishery_biom) <- c("SurveyID", "fishery_biomass")

#remove 0.1% outliers
surveys_fishery_biom <- surveys_fishery_biom %>%
  filter(fishery_biomass < quantile(surveys_fishery_biom$fishery_biomass,0.999))
  
# -------------------save data -------------------
save(surveys_fishery_biom, file = here::here("nutrients", "outputs", "fishery_tot_biomass.Rdata"))

###look co-variations with total biomass
biomass <- left_join(task3_data_surveys,surveys_fishery_biom)
plot(biomass$fishery_biomass ~ biomass$Btot)
