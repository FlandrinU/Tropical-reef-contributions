## cleaning memory
rm(list=ls())

# libraries
library(here)
library(tidyverse)
library(ggridges)
library(patchwork)

library(ggfortify)
library(corrmorant)

# data
load(here::here("results","surveys_biodiversity.Rdata"))
load(here::here("results","surveys_fluxes.Rdata"))
load(here::here("data","metadata_surveys.Rdata"))


# merging survey metadata with biodiversity and nutrient fluxes
surveys_merged <- metadata_surveys %>% 
  left_join( surveys_fluxes , by="SurveyID") %>%
  left_join( surveys_biodiversity, by="SurveyID" )
  
summary(surveys_merged)


# filtering to keep only surveys with
# >80% of biomass belonging to species with flux estimates
# and >80% of abundance belonging to species with flux estimates
# to not account for surveys with likely underestimates fluxes

surveys_merged <- surveys_merged %>%
  filter (p_biom_nutflux > 0.8 & p_abund_nutflux > 0.8)
  
# keeping only key variables about 
# total abundance and total biomass (of secies included in fluxes and biodiversiy estimates)
# biodiversity (taxonomic and functional richness and entropy)
# dominance of extreme trophic groups to total biomass
# size structure
# nutrient recycling, contribution of excretion to recycling, ratio between recycling and storing

  task3_dataforDAG <- surveys_merged %>%
  select(SurveyID, SiteCode,
         Ntot_fishflux, Btot_fishflux, 
         taxo_richness, taxo_entropy, funct_richness, funct_entropy,    
         size_Q1, size_median, size_Q3,    
         pbiom_lowTL, pbiom_highTL, 
         recycling_N, recycling_P,
         pexcr_recycling_N, pexcr_recycling_P, 
         recyc_stor_N, recyc_stor_P
         )

nrow(task3_dataforDAG) # 2 432 surveys
n_distinct(task3_dataforDAG$SiteCode ) # 1 522 sites


save(task3_dataforDAG, file=here::here("results", "task3_dataforDAG.Rdata")  )


  

