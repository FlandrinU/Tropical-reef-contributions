################################################################################################
#'
#'This script run all scripts of `R/` to obtain cleaned data
#'
#'@author Ulysse Flandrin, \email{ulysse.flandrin@@gmail.com}
#'
#'
#' @date 2022/10/01 first created
################################################################################################

rm(list=ls())

#-----------------run data preparation---------------------
##import raw data from google drive
#source(here::here("R", "0a_import_data.R"))    #-> need S. Villeger's permission, save data in "data_raw/source"

##merge rls data
cat("Merge and filter data... \n")
source(here::here("R", "0b_merge_datasets.R"))

##keep only tropical actinopterygians
source(here::here("R", "0c_filtering_tropical_fish.R"))

##remove small fishes and outlier sizes
source(here::here("R", "0d_filtering_size.R"))

##fill NAs in biomass estimation 
cat("Fill biomass estimation with fishflux package... \n")
source(here::here("R", "0e_filling_fish_biomass.R")) #takes few minutes due to the fishflux package (needs internet connection)

##rename and save final data
source(here::here("R", "0f_fish_datasets.R")) #keep only species existing in fishflux package

##extract data for elamsobranch species
source(here::here("R", "0g_elasmobranchii_datasets.R")) 

cat("Data are ready! \n")

