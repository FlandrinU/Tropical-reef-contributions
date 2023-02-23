#-----------------cleaning memory-------------------
rm(list=ls())

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse", "brms", "readr", "picante", "tidybayes", 
          "ggridges", "patchwork", "ggfortify", "corrmorant", "stats")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

#-----------------set directory---------------------
path <- here::here() #should arrive into the R_NCPs file
setwd(path)

#-----------------run recycling analyses---------------------
##fishflux data of every row of rls surveys
source(here::here("recycling", "R", "0g_run_fishflux.R"))   # /!\ long time to run (around 40')

##agregate flux estimations at the species and surveys scale
source(here::here("recycling", "R", "1b_surveys_fluxes.R")) 

##filter surveys with 80% of species computed by fishflux
source(here::here("recycling", "R", "1d_surveys_merging_filtering.R")) 

##summarize species fluxes
source(here::here("recycling", "R", "2_species_fluxes.R")) 

##agregate surveys data at the site scale
source(here::here("recycling", "R", "5_sites_fluxes.R")) 


