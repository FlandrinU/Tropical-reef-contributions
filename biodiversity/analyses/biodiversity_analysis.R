#################################################################################
#'
#'This script runs all scripts of `biodiversity/R/` to obtain the biodiversity
#' metrics (taxonomy, fonctional, phylogenetic, endemism) at the survey scale 
#' 
#'
#'@author Ulysse Falndrin, \email{ulysse.flandrin@@gmail.com}
#'
#'
################################################################################

#-----------------cleaning memory-------------------
rm(list=ls())
cat("Run biodiversity analysis... \n")

# #-----------------Loading packages-------------------
# pkgs <- c("here", "tidyverse", "mFD", "picante", "ape", "worrms", "taxize",
#           "phyloregion", "parallel", "FactoMineR", "tibble", "questionr",
#           "factoextra", "rredlist", "forcats", "entropart")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

#-----------------run biodiversity analyses---------------------

##Extract occurrence matrix, Assess taxonomic richness and entropy, functional richness
## and entropy, and trophic and size structure
source(here::here("biodiversity", "R", "occ_matrix_sprichness_FD.R")) #OK


## Assess phylogenetic diversty
source(here::here("biodiversity", "R", "phylogenetic_diversity.R")) #long time to run
  # Analyse phylogenetic indices
  source(here::here("biodiversity", "R", "phylogenetic_indices_correlations.R")) 

## Diversity of elasmobranchii and iucn species
source(here::here("biodiversity", "R", "iucn_and_elasmobranch_indices.R"))

## Assess mean endemism
source(here::here("biodiversity", "R", "endemism_survey_index.R"))

