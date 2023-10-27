################################################################################
##
## Give aesthetic scores of each RLS surveys, and others cultural contributions
##  (public interest and academic knowledge)
##
## cultural_analyses.R
##
## 25/10/2022
##
## Ulysse Flandrin
##
################################################################################

rm(list=ls())

# #-----------------Loading packages-------------------
# pkgs <- c("here", "tidyverse", "stats")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))


#-----------------load data---------------------
survey_aesth_all <- read.csv(here::here("cultural_contributions", "data", "survey_aesth.csv"))

#----------------- aesthetic score---------------------
survey_aesth_all$SurveyID <- as.character(survey_aesth_all$SurveyID)

#save
save(survey_aesth_all, file = here::here("cultural_contributions", "outputs", "survey_aesth.Rdata"))

#-----------------run public interest and academic knowledge analyses---------------------
source(here::here("cultural_contributions", "R", "public_and_scientific_interest_surveys.R"))
