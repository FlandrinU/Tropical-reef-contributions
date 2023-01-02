################################################################################
##
## Give aesthetic scores of each RLS surveys
##
## aesthetic_analyses.R
##
## 25/10/2022
##
## Ulysse Flandrin
##
################################################################################
#-----------------cleaning memory-------------------
rm(list=ls())

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse", "stats")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

#-----------------set directory---------------------
path <- here::here()
setwd(path)

#-----------------load data---------------------
survey_aesth_all <- read.csv(here::here("aesthetic", "outputs", "survey_aesth.csv"))

#-----------------run aesthetic analyses---------------------
survey_aesth_all$SurveyID <- as.character(survey_aesth_all$SurveyID)

#filter 0.01 outliers
survey_aesth <- survey_aesth_all %>%
  filter(aesthe_survey < quantile(survey_aesth_all$aesthe_survey,0.999))

#save
save(survey_aesth, file = here::here("aesthetic", "outputs", "survey_aesth_without_outliers.Rdata"))

