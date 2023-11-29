#################################################################################
#'
#'This script uses the results of nutrient concentration of fishes from 
#' Maire et al. 2021 (file: `nutrients/outputs/nutrient_sp_data.RData`) to
#' assess the mean concentration of nutrients per community. In addition, it 
#' uses the list of fished families in the world (Cinner et al. 2016) to extract
#' the fishery biomass.
#' 
#'
#'@author Ulysse Flandrin, \email{ulysse.flandrin@@gmail.com}
#'
#' @date 2022/11/24
#'
################################################################################

#-----------------cleaning memory-------------------
rm(list=ls())

#-----------------Loading all data---------------------
load(here::here("data", "metadata_surveys.Rdata"))
load(here::here("data", "data_surveys.Rdata"))
load(here::here("nutrients", "outputs", "nutrient_sp_data.RData"))


#----------------- biomass of exploited species ---------------------
source(here::here("nutrients", "R", "fishery_species_biomass.R"))


#-----------------Aggregating nutrients concentration at the surveys level---------------------
#aggregate by species in each transect
RLS_data <- data_surveys |>
  dplyr::group_by(SurveyID, species) |>
  dplyr::summarise(biomass_sp = sum(biomass)) |> 
  dplyr::filter(species %in% unique(data_surveys_fishery$species)) # keep only species targeted by fisheries

#nutrients per species per transect
RLS_nut_sp_surv <- RLS_data |>
  dplyr::left_join(nutrientdata) |>
  dplyr::select(SurveyID, species, biomass_sp, Selenium_mu, Zinc_mu, Omega_3_mu, Calcium_mu, Iron_mu, Vitamin_A_mu) |>
  dplyr::mutate(Selenium_tot  = Selenium_mu  * biomass_sp,
                Zinc_tot      = Zinc_mu      * biomass_sp,
                Omega_3_tot   = Omega_3_mu   * biomass_sp,
                Calcium_tot   = Calcium_mu   * biomass_sp,
                Iron_tot      = Iron_mu      * biomass_sp,
                Vitamin_A_tot = Vitamin_A_mu * biomass_sp)

#sum of nutrients in each transect, and concentration of it  
RLS_nut_surv <- RLS_nut_sp_surv |>
  dplyr::group_by(SurveyID) |>
  dplyr::summarise(biom_tot    = sum(biomass_sp),
                   Selenium_C  = sum(Selenium_tot) / biom_tot,
                   Zinc_C      = sum(Zinc_tot)     / biom_tot,
                   Omega_3_C   = sum(Omega_3_tot)  / biom_tot,
                   Calcium_C   = sum(Calcium_tot)  / biom_tot,
                   Iron_C      = sum(Iron_tot)     / biom_tot,
                   Vitamin_A_C = sum(Vitamin_A_tot)/ biom_tot )


save(RLS_nut_surv, file = here::here("nutrients", "outputs", "nutrient_concentration_surveys.Rdata"))
