################################################################################
##
## Assess evolutionnary distinctivness of reef species from trees of Siqueira 
## et al. 2020, and mean these ED at the surveys scale
## Assess phylogenetic diversity and phylogenetic endemism of surveys
##
## phylogenetic_diversity.R
##
## 19/10/2022
##
## Ulysse Flandrin
##
################################################################################
## cleaning memory
rm(list=ls())

##-------------loading data-------------
#datasets
load(here::here("data","data_species.Rdata"))
load(here::here("data","data_surveys.Rdata"))
#trees
trees <- ape::read.tree(file = here::here("biodiversity", "data","Code&Data_Siquiera2020",
                                          "TACT", "Reef_fish_all_combined.trees")) #chronogram of ray-finned fishes, Siquiera 2020 from Rabosky 2018
#occurrence matrix
load(here::here("biodiversity", "outputs","occurrence_matrix_sp_survey_01.Rdata"))

#relative biomass matrix
load(here::here("biodiversity", "outputs","occurrence_matrix_sp_survey_relative_biomass.Rdata"))

## -------------extract trees with RLS reef species-------------
## Change names in matrix to fit with names in phylogeny
treesp <- trees[[1]][["tip.label"]]
names <- rep(NA, ncol(surveys_sp_occ))
sp_wrong_name <- c()
for( i in 1:ncol(surveys_sp_occ)){
  if(colnames(surveys_sp_occ)[i] %in% treesp){
    names[i] <- colnames(surveys_sp_occ)[i]
    
  }else if(data_species$species_corrected[
    which(data_species$species == colnames(surveys_sp_occ)[i])] %in% treesp){
    names[i] <- data_species$species_corrected[which(data_species$species == colnames(surveys_sp_occ)[i])]
    
  }else{
    sp_wrong_name <- c(sp_wrong_name, colnames(surveys_sp_occ)[i])
    cat(colnames(surveys_sp_occ)[i], " not in phylogeny \n")
  }
}

#Is there a synonym names in the tree?
for( sp in sp_wrong_name){
  sp <- stringr::str_replace(sp, "_", " ")
  tab_syn <- taxize::synonyms(sp , db="worms")
  if( nrow(tab_syn[[sp]]) > 0){
    sp_syn <- stringr::str_replace(tab_syn[[sp]][["scientificname"]], " ", "_")
    cat(length(sp_syn), "synonyms", "\n" )
    cat( "new name : ", sp_syn[which(is.element(sp_syn, treesp)==T)], "\n" )
  }
} #No synonyms found in the tree

colnames(surveys_sp_occ) <- colnames(surveys_sp_pbiom) <- names

surveys_sp_occ <- as.data.frame(surveys_sp_occ) |>
  dplyr::select( dplyr::contains( "_")) #Remove species absent form phylogeny

surveys_sp_pbiom <- as.data.frame(surveys_sp_pbiom) |>
  dplyr::select( dplyr::contains( "_")) |> #Remove species absent form phylogeny
  as.matrix()


#Extract trees
phylo_100<-list()
for (i in 1:100) {
  phylo_100[[i]]<-picante::match.phylo.comm(trees[[i]], surveys_sp_occ)
}

#Occurrence matrix in RLS surveys
occ_matrix <- as.matrix(phylo_100[[1]][["comm"]])

##-------------Compute phylogenetic indices-------------

## Evolutionary distinctivness: ED  Isaac et al. method
#by species
ED_species_raw<-parallel::mclapply(phylo_100, mc.cores=20, function(x) {
  picante::evol.distinct(x$phy, type = c("fair.proportion"), scale = FALSE, use.branch.lengths = TRUE)})

ED_species <- as.data.frame(do.call(cbind, lapply(ED_species_raw,function(y){y[,2]})))
rownames(ED_species) <- phylo_100[[1]]$phy$tip.label

# save(ED_species, file = here::here("biodiversity", "outputs", "evolutionary_distinctivness_species_100trees.Rdata"))

ED_species_summary <- t(apply(ED_species, 1, summary))
ED_species_summary <- cbind( ED_species_summary, sd = apply(ED_species, 1, sd))
colnames(ED_species_summary) <- paste0("ED_", colnames(ED_species_summary))

save(ED_species_summary, file = here::here("biodiversity", "outputs", "evolutionary_distinctivness_species.Rdata"))

#by surveys
mean_ED_sp <- apply(ED_species, 1, mean)
ED_surveys_raw <- parallel::mclapply(1:nrow(occ_matrix), mc.cores=15 ,function(i){
  if(sum(occ_matrix[i,])==0){
    rep(0,6)
  }else{
    Sp <- colnames(occ_matrix) [which(occ_matrix[i,]==1)]
    summary(mean_ED_sp[Sp])
  }
})

ED_surveys <- do.call(rbind,ED_surveys_raw)
row.names(ED_surveys) <- row.names(occ_matrix)
colnames(ED_surveys) <- paste0("ED_", colnames(ED_surveys))

save(ED_surveys, file = here::here("biodiversity", "outputs", "evolutionary_distinctivness_surveys.Rdata"))



## PD and SES.PD quantification with Phyloregion
#PD
sparse_occ_matrix <- methods::as( as.matrix(occ_matrix), "sparseMatrix")
PD_surveys_raw <- lapply(phylo_100, function (x){ phyloregion::PD(sparse_occ_matrix, x$phy) })
PD_surveys_100 <- t( do.call(rbind, PD_surveys_raw ) )

PD_surveys_summary <- t(apply(PD_surveys_100, 1, summary))
colnames(PD_surveys_summary) <- paste0("PD_", colnames(PD_surveys_summary))


#Residuals of PD ~ taxonomic richness
Mean_PD <- apply(PD_surveys_100,1,mean)
taxo_richness <- apply(occ_matrix, 1, sum)
residuals_PD_richness <- residuals(lm(Mean_PD ~ taxo_richness)) # Approach using residuals

#SES.PD
SES_PD_raw <- parallel::mclapply(phylo_100, mc.cores= 20, function (x){
  phyloregion::PD_ses(sparse_occ_matrix, x$phy, model = c("tipshuffle"), reps= 1000)$zscore #Change reps for a quickier analyse
}) # Approach using SES.PD based on a null model shuffling tip labels                        

SES_PD_surveys_100 <- do.call(cbind, SES_PD_raw)
SES_PD_surveys_summary <- t(apply(SES_PD_surveys_100, 1, summary))
colnames(SES_PD_surveys_summary) <- paste0("SES_PD_", colnames(SES_PD_surveys_summary))


PD_surveys <- cbind( PD_surveys_summary, residuals_PD_richness = residuals_PD_richness, SES_PD_surveys_summary)
save(PD_surveys, file = here::here("biodiversity", "outputs", "phylogenetic_diversity_surveys.Rdata"))

## PE: Phylogenetic endemism (pkg Pyloregion)
PE_surveys_raw <- parallel::mclapply(phylo_100, mc.cores=15, function(x) {
  phyloregion::phylo_endemism(sparse_occ_matrix, x$phy, weighted = TRUE)
  })

PE_surveys_100 <- do.call(cbind, PE_surveys_raw)
PE_surveys_summary <- t(apply(PE_surveys_100, 1, summary))
colnames(PE_surveys_summary) <- paste0("PE_", colnames(PE_surveys_summary))

save(PE_surveys_summary, file = here::here("biodiversity", "outputs", "phylogenetic_endemism_surveys.Rdata"))


## phylogenetic entropy: Marcon 2015, from Allen 2009  -> very long time to run: run only on 10 trees due to negligible variations
library('entropart')
phylo_entropy_raw <- lapply( phylo_100[1:10], function(x) {
  cat("Start computation for one tree... \n")
  list <- parallel::mclapply(rownames(surveys_sp_pbiom), mc.cores = 21, function(survey){
    community_biom <- surveys_sp_pbiom[survey,]
    phylo_entropy <- entropart::PhyloEntropy(Ps=community_biom,
                                             q = 1,
                                             Tree = x[["phy"]], #must be written x[[1]][["phy"]] without lapply
                                             Normalize = T) #phylogenetic entropy of order 1 of relative biomass vector.
    phylo_entropy[["Total"]]
  })
  cat("one more tree done.  \n")
  do.call(rbind, list)
})


phylo_entropy_surveys_10 <- do.call(cbind, phylo_entropy_raw)
  # test variability
  sd <- apply(phylo_entropy_surveys_10, 1, sd)
  summary(sd) # Median 0.0000918      Mean 0.0011747  3rd Qu. 0.0010619      Max. 0.0405119
  #
  
phylo_entropy_summary <- t(apply(phylo_entropy_surveys_10, 1, summary))
colnames(phylo_entropy_summary) <- paste0("phylo_entropy_", colnames(phylo_entropy_summary))
rownames(phylo_entropy_summary) <- rownames(surveys_sp_pbiom)

save(phylo_entropy_surveys_10, file = here::here("biodiversity", "outputs", "phylogenetic_entropy_surveys_10trees.Rdata"))
save(phylo_entropy_summary, file = here::here("biodiversity", "outputs", "phylogenetic_entropy_surveys.Rdata"))

##-------------Save phylogenetic indices-------------
phylo_indices_surveys_all <- cbind(ED_surveys, PD_surveys, PE_surveys_summary, phylo_entropy_summary )
phylo_indices_surveys_all <- as.data.frame(phylo_indices_surveys_all)
phylo_indices_surveys_all <- tibble::rownames_to_column(phylo_indices_surveys_all,var="SurveyID")

# # filtering 0.1% outliers
# phylo_indices_surveys <- phylo_indices_surveys_all %>%
#   filter(ED_Mean < quantile(phylo_indices_surveys_all$ED_Mean,0.999)) %>% 
#   #filter(residuals_PD_richness < quantile(phylo_indices_surveys_all$residuals_PD_richness,0.999)) %>% 
#   #filter(SES_PD_Mean < quantile(phylo_indices_surveys_all$SES_PD_Mean,0.999)) %>% 
#   #filter(PE_Mean < quantile(phylo_indices_surveys_all$PE_Mean,0.999)) %>%
#   filter(phylo_entropy < quantile(phylo_indices_surveys_all$phylo_entropy, 0.999)) 
# 

save(phylo_indices_surveys_all, file = here::here("biodiversity", "outputs", "phylogenetic_indices_surveys.Rdata"))
# save(phylo_indices_surveys, file = here::here("biodiversity", "outputs", "phylogenetic_indices_surveys_without_outliers.Rdata"))

