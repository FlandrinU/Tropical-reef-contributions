#-----------------cleaning memory-------------------
rm(list=ls())

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse", "rfishbase", "dplyr", "DHARMa", "gtools", 
          "parallel", "car", "PresenceAbsence",  "igraph", "NetIndices",
          "ggplot2", "ggsignif" )
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

#-----------------set directory---------------------
path <- here::here() 
setwd(path)

#-----------------run trophic web analyses---------------------
# Extract species data from fishbase to infer common size o species
source(here::here("trophic_web", "R", "common_length_and_TL_of_species.R"))

# Construct the metaweb from species size (and diet), with model calibrated with Barnes 2008 data
source(here::here("trophic_web", "R", "metaweb_construction.R"))

# Test the confidence of the size niche model with a Boyce index
source(here::here("trophic_web", "R", "cross_validation_niche_model.R"))

# Extract local trophic web and trophic indicators of each sites
source(here::here("trophic_web", "R", "extract_local_web_trophic_indicators.R"))

# Figures of trophic web analysis
source(here::here("trophic_web", "R", "metaweb_figures.R"))



