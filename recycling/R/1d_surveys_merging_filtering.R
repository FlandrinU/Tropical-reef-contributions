#################################################################################
#'
#' Cleans data at the survey scale. Surveys for which less than 80% of 
#'  individuals are considered in the nutrient flow assessement are deleted.
#' 
#'
#'@author Sebastien Villéger
#'
#'
################################################################################

## cleaning memory
rm(list=ls())
library(ggplot2)

# data
load(here::here("recycling", "outputs", "surveys_fluxes.Rdata"))
load(here::here("data","metadata_surveys.Rdata"))


# merging survey metadata with biodiversity and nutrient fluxes
surveys_merged <- metadata_surveys |> 
  dplyr::left_join( surveys_fluxes , by="SurveyID")
  
summary(surveys_merged)


# filtering to keep only surveys
# >80% of biomass belonging to species with flux estimates
# >80% of abundance belonging to species with flux estimates
# total biomass < 500 kg 
# total abundance < 10 000
# without outliers 0.1%

task3_data <- surveys_merged |>
  dplyr::filter(Btot < 500000) |>
  dplyr::filter(Ntot < 10000) |>
  dplyr::filter (p_biom_nutflux > 0.8 ) |>
  dplyr::filter (p_abund_nutflux > 0.8 ) 

nrow(task3_data) # 2 402 surveys
dplyr::n_distinct(task3_data$SiteCode ) # 1 505 sites

# #Removing outliers, N and P recycling values superior to 99.9% of values
# task3_data_without_outliers <-  task3_data |>
#   dplyr::filter( recycling_N < quantile(task3_data$recycling_N, 0.999) ) |>
#   dplyr::filter( recycling_P < quantile(task3_data$recycling_P, 0.999) )


# log-transforming abundance and biomass
task3_data_surveys <- task3_data |>
  dplyr::mutate(log10_biomass=log10(Btot_fishflux) ) |>
  dplyr::mutate(log10_density=log10(Ntot_fishflux) )

summary(task3_data_surveys)

dplyr::n_distinct(task3_data_surveys$SurveyID ) # 2399 unique surveys

save(task3_data_surveys, file=here::here("recycling", "outputs", "flux_final_data_surveys.Rdata")  )

list_surveys_recycling <- task3_data_surveys$SurveyID
#save(list_surveys_recycling, file=here::here("data", "list_surveys_recycling.Rdata")  )