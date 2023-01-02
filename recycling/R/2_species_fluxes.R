## cleaning memory
rm(list=ls())

# datasets from surveys
load(here::here("recycling", "outputs", "surveys_species_fluxes.Rdata"))
load(here::here("recycling", "outputs", "flux_final_data_surveys.Rdata"))

#list of filtered surveys
surveys <- task3_data_surveys$SurveyID
length(surveys) #2042

flux_sp <- data.frame(species = colnames(surveys_species_fluxes[["recycling_P"]]))
variables <- c("recycling_P", "recycling_N")

#mass of nutrients recycled per gramme of each species
for ( k in variables){
  mat <- surveys_species_fluxes[[k]]
  mat <- mat[surveys,] #only selected surveys
  biom <- surveys_species_fluxes[["biomass_estimated"]][surveys,]
  
  recyc_per_g <- mat / biom
  
  m <- data.frame(colMeans(recyc_per_g, na.rm=T)) ; colnames(m) <- k 
  flux_sp <- cbind(flux_sp, m)
}

summary(flux_sp) #24 NAs (species absent in the selected surveys)

#Intrinsic ability of species in recycling N and P per unit of biomass, overlooking environnemental context
# recycling of P and N in g of nutrients per g of fish
save(flux_sp, file = here::here("recycling", "outputs", "flux_final_data_species.Rdata") )
