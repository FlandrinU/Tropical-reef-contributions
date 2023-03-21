#-----------------cleaning memory-------------------
rm(list=ls())

#-----------------Loading packages-------------------
pkgs <- c("here", "dplyr", "stats")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

#-----------------Loading all data---------------------
load(here::here("data", "metadata_surveys.Rdata"))
load(here::here("data", "data_surveys.Rdata"))
load(here::here("nutrients", "outputs", "nutrient_sp_data.RData"))



#----------------- biomass of exploited species ---------------------
source(here::here("nutrients", "R", "fishery_species_biomass.R"))



#-----------------Aggregating nutrients concentration at the surveys level---------------------
#aggregate by species in each transect
RLS_data <- data_surveys %>%
  group_by(SurveyID, species) %>%
  summarise(biomass_sp = sum(biomass))

#nutrients per species per transect
RLS_nut_sp_surv <- RLS_data %>%
  left_join(nutrientdata) %>%
  select(SurveyID, species, biomass_sp, Selenium_mu, Zinc_mu, Omega_3_mu, Calcium_mu, Iron_mu, Vitamin_A_mu) %>%
  mutate(Selenium_tot  = Selenium_mu  * biomass_sp,
         Zinc_tot      = Zinc_mu      * biomass_sp,
         Omega_3_tot   = Omega_3_mu   * biomass_sp,
         Calcium_tot   = Calcium_mu   * biomass_sp,
         Iron_tot      = Iron_mu      * biomass_sp,
         Vitamin_A_tot = Vitamin_A_mu * biomass_sp)

#sum of nutrients in each transect, and concentration of it  
RLS_nut_surv <- RLS_nut_sp_surv %>%
  group_by(SurveyID) %>%
  summarise(biom_tot    = sum(biomass_sp),
            Selenium_C  = sum(Selenium_tot) / biom_tot,
            Zinc_C      = sum(Zinc_tot)     / biom_tot,
            Omega_3_C   = sum(Omega_3_tot)  / biom_tot,
            Calcium_C   = sum(Calcium_tot)  / biom_tot,
            Iron_C      = sum(Iron_tot)     / biom_tot,
            Vitamin_A_C = sum(Vitamin_A_tot)/ biom_tot )

# save(RLS_nut_surv, file = here::here("nutrients", "outputs", "nutrient_concentration_all_surveys.Rdata"))
# 
# 
# #Remove outliers 0.1%
# out_S <- which((RLS_nut_surv$Selenium_C < quantile(RLS_nut_surv$Selenium_C, 0.999)) == F)
# out_Z <- which((RLS_nut_surv$Zinc_C < quantile(RLS_nut_surv$Zinc_C, 0.999)) == F)
# out_O <- which((RLS_nut_surv$Omega_3_C < quantile(RLS_nut_surv$Omega_3_C, 0.999)) == F)
# out_C <- which((RLS_nut_surv$Calcium_C < quantile(RLS_nut_surv$Calcium_C, 0.999)) == F)
# out_I <- which((RLS_nut_surv$Iron_C < quantile(RLS_nut_surv$Iron_C, 0.999)) == F)
# out_V <- which((RLS_nut_surv$Vitamin_A_C < quantile(RLS_nut_surv$Vitamin_A_C, 0.999)) == F)
# outliers <- unique(c(out_S, out_Z, out_O, out_C, out_I, out_V))
# 
# RLS_nut_without_outliers <- RLS_nut_surv[ -outliers, ]
save(RLS_nut_surv, file = here::here("nutrients", "outputs", "nutrient_concentration_surveys.Rdata"))

#-----------------Aggregating at the site level---------------------
RLS_nut_site <- RLS_nut_surv %>%
  left_join(metadata_surveys)%>%
  group_by(SiteCode) %>%
  summarise(biom_tot    = mean(biom_tot),
            Selenium_C  = mean(Selenium_C),
            Zinc_C      = mean(Zinc_C),
            Omega_3_C   = mean(Omega_3_C),
            Calcium_C   = mean(Calcium_C),
            Iron_C      = mean(Iron_C),
            Vitamin_A_C = mean(Vitamin_A_C))

save(RLS_nut_site, file = here::here("nutrients", "outputs", "nutrient_concentration_sites.Rdata"))

