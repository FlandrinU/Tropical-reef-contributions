################################################################################
##
## Makes t-SNE analyses on all sites, with all NCPs. 
##
## tSNE_analyses_on_NCP.R
##
## 09/01/2023
##
## Ulysse Flandrin
##
################################################################################

#----------------- Loading packages -------------------
pkgs <- c("here", "tidyverse", "tibble", "questionr", "dplyr", "Rtsne")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##------------- loading data -------------
load(here::here("outputs","all_NCP_site.Rdata"))
load(here::here("data","metadata_surveys.Rdata"))

##------------- preping data -------------
#NCP_site <- dplyr::filter(NCP_site, SiteCountry != "Australia")
NCP_site_clean <- subset(NCP_site, select = -c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                                 SiteMeanSST, SiteLatitude, SiteLongitude, Btot,
                                                 HDI, gravtot2, MarineEcosystemDependency,
                                                 coral_imputation))
NCP_unique <- unique (NCP_site_clean)
NCP_site_scaled <- Rtsne::normalize_input( as.matrix(NCP_unique) )

##------------- computing t-SNE -------------

NCP_tsne <- Rtsne::Rtsne(NCP_site_scaled, dims = 2)
plot(NCP_tsne$Y[,c(1,2)], col=NCP_site$SiteCountry, asp=1)

