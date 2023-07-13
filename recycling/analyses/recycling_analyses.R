#################################################################################
#'
#'This script runs all scripts of `recycling/R/` to obtain N and P recycling 
#' flows per species and per communities from RLS data and the rfishflux package
#' 
#'
#'@author Ulysse Flandrin, \email{ulysse.flandrin@@gmail.com}
#'
#'
#' @date 2022/10/01 first created
################################################################################

#-----------------cleaning memory-------------------
rm(list=ls())
cat("Run recycling analysis... \n")

#-----------------run recycling analyses---------------------
## fishflux data of every row of rls surveys
source(here::here("recycling", "R", "0g_run_fishflux.R"))   # /!\ long time to run (around 40')

## agregate flux estimations at the species and surveys scale
source(here::here("recycling", "R", "1b_surveys_fluxes.R")) 

## filter surveys with 80% of species computed by fishflux
source(here::here("recycling", "R", "1d_surveys_merging_filtering.R")) 

## summarize species fluxes
source(here::here("recycling", "R", "2_species_fluxes.R")) 

## agregate surveys data at the site scale
source(here::here("recycling", "R", "5_sites_fluxes.R")) 


