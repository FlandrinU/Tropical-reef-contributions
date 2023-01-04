################################################################################
##
## Makes RDA analyses on all NCPs
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
          "rdacca.hp")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##------------- loading data -------------
load(here::here("outputs","all_NCP_site.Rdata"))
load(here::here("biodiversity", "outputs", "occurrence_matrix_family_survey_relative_biomass.Rdata"))
load(here::here("data","metadata_surveys.Rdata"))

##------------- preping data -------------
rownames(NCP_site) <- NCP_site$SiteCode
NCP_site_for_rda <- subset(NCP_site, select = -c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                                 SiteMeanSST, SiteLatitude, SiteLongitude, Btot,
                                                 HDI, gravtot2, MarineEcosystemDependency,
                                                 coral_imputation))

##------------- computing rda -------------

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

summary(site_family_pbiom)


#reduce data
NCP_site_for_rda <- NCP_site_for_rda[c(1:1000), c(1:12)]
site_family_pbiom <- site_family_pbiom[c(1:1000),c(1:15)]

#compute RDA
rdacca.hp::rdacca.hp(dv = NCP_site_for_rda, 
                     iv = site_family_pbiom,
                     method = "RDA", type = "adjR2",
                     scale = TRUE)



library(vegan)
data(mite)
data(mite.env)
data(mite.xy)
data(mite.pcnm)
#Hellinger-transform the species dataset for RDA
mite.hel <- decostand(mite, "hellinger")
rdacca.hp(mite.hel,mite.env,method="RDA",type="adjR2")
rdacca.hp(vegdist(mite),mite.env,method="dbRDA",type="adjR2")
rdacca.hp(mite,mite.env,method="CCA",type="adjR2")
iv <- list(env=mite.env,xy=mite.xy,pcnm=mite.pcnm)
rdacca.hp(mite.hel,iv,method="RDA",var.part = TRUE)
rdacca.hp(vegdist(mite),iv,method="dbRDA",var.part = TRUE)
rdacca.hp(mite,iv,method="CCA",var.part = TRUE)
