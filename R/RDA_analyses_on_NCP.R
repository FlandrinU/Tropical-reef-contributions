################################################################################
##
## Makes RDA analyses on all NCPs. /!\ the rda analysis with rdacca.hp function
## requires a very large virtual memory and a long time to run.
##
## RDA_analyses_on_NCP.R
##
## 03/01/2023
##
## Ulysse Flandrin
##
################################################################################

#----------------- Loading packages -------------------
pkgs <- c("here", "tidyverse","FactoMineR", "tibble", "questionr", "corrplot",
          "factoextra", "ggpubr", "scico", "RColorBrewer", "plotly", "fishualize", 
          "ggplot2", "patchwork", "colormap", "grDevices", "ggnewscale", "sf",
          "rdacca.hp", "dplyr")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##------------- loading data -------------
load(here::here("outputs","all_NCP_site.Rdata"))
load(here::here("biodiversity", "outputs", "occurrence_matrix_family_survey_relative_biomass.Rdata"))
load(here::here("biodiversity", "outputs", "occurrence_matrix_family_survey_01.Rdata"))
load(here::here("biodiversity", "outputs", "occurrence_matrix_nb_sp_per_family_survey.Rdata"))
load(here::here("data","metadata_surveys.Rdata"))

##------------- preping data -------------
rownames(NCP_site) <- NCP_site$SiteCode
NCP_site_clean <- subset(NCP_site, select = -c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                                 SiteMeanSST, SiteLatitude, SiteLongitude, Btot,
                                                 HDI, gravtot2, MarineEcosystemDependency,
                                                 coral_imputation))


##------------- studying family importance -------------
family_tot_presence <- colSums(surveys_family_occ)
family_tot_presence[order(family_tot_presence, decreasing = T)] 
# 7 rarest families: Ostraciidae > Pempheridae > Fistulariidae > Scorpaenidae > Sciaenidae > Mugilidae > Bothidae 

colSums(surveys_nb_sp_per_family)[order(colSums(surveys_nb_sp_per_family), decreasing = T)]
colSums(surveys_family_pbiom)[order(colSums(surveys_family_pbiom), decreasing = T)]


##------------- computing rda with family occurrence-------------
#Sum up relative biomass of each family at the site scale
surveys_family_occ <- as.data.frame(surveys_family_occ) %>%
  tibble::rownames_to_column("SurveyID") %>%
  dplyr::left_join(metadata_surveys[,c("SurveyID", "SiteCode")]) %>%
  dplyr::select(- SurveyID) 

site_family_occ <- surveys_family_occ %>%
  group_by( SiteCode) %>%
  summarise(across(.cols = everything(),
                   .fns = max, 
                   .names = "{.col}")) %>%
  dplyr::filter(SiteCode %in% NCP_site$SiteCode) %>%
  dplyr::select( - SiteCode)

#summary(site_family_occ)


#reduce data
NCP_site_for_rda <- scale(NCP_site_clean[,] )
site_family_occ <- subset(site_family_occ, 
                          select = -c(Ostraciidae, Pempheridae, Fistulariidae,
                                      Scorpaenidae, Sciaenidae, Mugilidae, 
                                      Bothidae)) #without rare families


#compute RDA
cat("Computing RDA with family occurences...", "\n")
rda_family_occ <- rdacca.hp::rdacca.hp(dv = NCP_site_for_rda, 
                     iv = site_family_occ,
                     method = "RDA", type = "adjR2",
                     scale = F)

save(rda_family_occ, file = here::here("outputs", "RDA_NCP_explained_by_19_occ_family.Rdata") )
cat("End of computation for RDA with family occurences", "\n")


##------------- computing rda with family richness-------------
#Sum up relative biomass of each family at the site scale
surveys_family_richness <- as.data.frame(surveys_nb_sp_per_family) %>%
  tibble::rownames_to_column("SurveyID") %>%
  dplyr::left_join(metadata_surveys[,c("SurveyID", "SiteCode")]) %>%
  dplyr::select(- SurveyID) 

site_family_richness <- surveys_family_richness %>%
  group_by( SiteCode) %>%
  summarise(across(.cols = everything(),
                   .fns = max, 
                   .names = "{.col}")) %>%
  dplyr::filter(SiteCode %in% NCP_site$SiteCode) %>%
  dplyr::select( - SiteCode)

#summary(site_family_richness)


# #reduce data
# NCP_site_for_rda <- scale(NCP_site_clean[,] )
# site_family_richness <- site_family_richness[,]

#compute RDA
cat("Computing RDA with family richness...", "\n")
rda_family_richness <- rdacca.hp::rdacca.hp(dv = NCP_site_for_rda, 
                                       iv = site_family_richness,
                                       method = "RDA", type = "adjR2",
                                       scale = F)

save(rda_family_richness, file = here::here("outputs", "RDA_NCP_explained_by_19_richness_family.Rdata") )
cat("End of computation for RDA with family richness", "\n")


##------------- computing rda with relative biomass-------------

#Sum up relative biomass of each family at the site scale
surveys_family_pbiom <- as.data.frame(surveys_family_pbiom) %>%
  tibble::rownames_to_column("SurveyID") %>%
  dplyr::left_join(metadata_surveys[,c("SurveyID", "SiteCode")]) %>%
  dplyr::select(- SurveyID) 

site_family_pbiom <- surveys_family_pbiom %>%
  group_by( SiteCode) %>%
  summarise(across(.cols = everything(),
                   .fns = mean, 
                   .names = "{.col}")) %>%
  dplyr::filter(SiteCode %in% NCP_site$SiteCode) %>%
  dplyr::select( - SiteCode)

#summary(site_family_pbiom)


# #reduce data
# NCP_site_for_rda <- scale(NCP_site_clean[,] )
# site_family_pbiom <- site_family_pbiom[,]

#compute RDA
cat("Computing RDA with family relative biomass...", "\n")
rda_family_pbiom <- rdacca.hp::rdacca.hp(dv = NCP_site_for_rda, 
                                         iv = site_family_pbiom,
                                         method = "RDA", type = "adjR2",
                                         scale = F)

save(rda_family_pbiom, file = here::here("outputs", "RDA_NCP_explained_by_19_pbiom_family.Rdata") )
cat("End of computation for RDA with family relative biomass", "\n")

##------------- computing rda with relative biomass only on NN or NS NCP-------------
#reduce data
NN <- c("recycling_N", "recycling_P", "taxo_richness", "funct_entropy", "funct_distinctiveness",
        "phylo_entropy", "ED_Mean", "iucn_species", "elasmobranch_diversity", "low_mg_calcite",
        "high_mg_calcite", "aragonite", "monohydrocalcite", "amorphous_carbonate", "biom_lowTL", 
        "biom_mediumTL", "biom_highTL") 
NS <- c("Productivity", "Selenium_C", "Zinc_C", "Omega_3_C", "Calcium_C", "Iron_C", "Vitamin_A_C",
        "aesthe_survey", "fishery_biomass")

NN_site_for_rda <- scale(NCP_site_clean[,NN] )
NS_site_for_rda <- scale(NCP_site_clean[,NS] )
#site_family_pbiom <- site_family_pbiom[,]

#compute RDAs
cat("Computing RDA with family relative biomass on NN NCPs...", "\n")
rda_family_pbiom_NN <- rdacca.hp::rdacca.hp(dv = NN_site_for_rda, 
                                         iv = site_family_pbiom,
                                         method = "RDA", type = "adjR2",
                                         scale = F)

save(rda_family_pbiom_NN, file = here::here("outputs", "RDA_NN_explained_by_19_pbiom_family.Rdata") )
cat("End of computation for RDA with family relative biomass on NN NCPs", "\n")

cat("Computing RDA with family relative biomass on NS NCPs...", "\n")
rda_family_pbiom_NS <- rdacca.hp::rdacca.hp(dv = NS_site_for_rda, 
                                            iv = site_family_pbiom,
                                            method = "RDA", type = "adjR2",
                                            scale = F)

save(rda_family_pbiom_NS, file = here::here("outputs", "RDA_NS_explained_by_19_pbiom_family.Rdata") )
cat("End of computation for RDA with family relative biomass on NS NCPs", "\n")