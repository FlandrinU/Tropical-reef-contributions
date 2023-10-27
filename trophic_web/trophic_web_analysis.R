#################################################################################
#'
#'This script infers the trophic interaction metaweb among all the RLS species
#' and extract different trophic indicators of the local webs, defined at the 
#' survey level.
#' 
#'
#'@author Ulysse Flandrin, \email{ulysse.flandrin@@gmail.com}
#'
#' @date 2023/02/08
#'
################################################################################
#-----------------cleaning memory-------------------
rm(list=ls())

# #-----------------Loading packages-------------------
# pkgs <- c("here", "tidyverse", "rfishbase", "dplyr", "DHARMa", "gtools", 
#           "parallel", "car", "PresenceAbsence",  "igraph", "NetIndices",
#           "ggplot2", "ggsignif" )
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))


#-----------------run trophic web analyses---------------------
# Extract species data from fishbase to infer common size of species
source(here::here("trophic_web", "R", "common_length_and_TL_of_species.R")) 

# Construct the metaweb from species size (and diet), with model calibrated with 
 # observed interactions from Barnes et al. 2008 data
source(here::here("trophic_web", "R", "metaweb_construction.R")) 

# Test the confidence of the size niche model with a Boyce index
source(here::here("trophic_web", "R", "cross_validation_niche_model.R")) 

# Extract local trophic web and trophic indicators of each sites
source(here::here("trophic_web", "R", "extract_local_web_trophic_indicators.R")) 

# Figures of trophic web analysis -> not necessary, here to evaluate the metaweb
source(here::here("trophic_web", "R", "metaweb_figures.R")) 



