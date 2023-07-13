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

# data
load(here::here("recycling", "outputs", "surveys_fluxes.Rdata"))
load(here::here("data","metadata_surveys.Rdata"))


# merging survey metadata with biodiversity and nutrient fluxes
surveys_merged <- metadata_surveys |> 
  left_join( surveys_fluxes , by="SurveyID")
  
summary(surveys_merged)

# keeping only key variables about nutrient recycling and storing
task3_allsurveys_keyvariables <- surveys_merged |>
  select(SurveyID, SurveyDate, SurveyDepth, 
         SiteCode,SiteLatitude, SiteLongitude, SiteCountry, SiteEcoregion, SiteMeanSST,
         Ntot, Btot, Ntot_fishflux,Btot_fishflux, p_abund_nutflux, p_biom_nutflux,
         recycling_N, recycling_P,
         pexcr_recycling_N, pexcr_recycling_P, 
         recyc_stor_N, recyc_stor_P
         )

nrow(task3_allsurveys_keyvariables) # 3 627 surveys
head(task3_allsurveys_keyvariables) # 3 627 surveys

save(task3_allsurveys_keyvariables, file=here::here("recycling", "outputs", "flux_allsurveys_keyvariables.Rdata")  )


# plot of total biomass and abundance
var_plot<-c("Btot_fishflux",  "Ntot_fishflux") 

for (k in 1:length(var_plot) ) {
  
  var_k<-var_plot[k]
  
  plot_k<- surveys_merged |>
    mutate(type="survey") |>
    ggplot( aes(x = .data[[var_k]], y = type, fill = factor(stat(quantile))) ) +
    ggridges::stat_density_ridges( geom = "density_ridges_gradient", calc_ecdf = TRUE,
                         quantiles = c(0.01, 0.25, 0.5, 0.75, 0.99), quantile_lines = TRUE) +
    scale_fill_viridis_d(name = "Deciles") +
    scale_x_log10() +
    scale_y_discrete(expand = c(0,0)) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none") 
  # plot_k
  if(k==1) {to_plot <- plot_k} else { to_plot <- to_plot + plot_k }
  
}# end of k

to_plot <- to_plot + 
  patchwork::plot_layout(guides="collect", ncol=1)
ggsave(here::here("recycling", "figures", "histogram_surveys_abundbiom_raw.png"))



# filtering to keep only surveys
# >80% of biomass belonging to species with flux estimates
# >80% of abundance belonging to species with flux estimates
# total biomass < 500 kg 
# total abundance < 10 000
# without outliers 0.1%

task3_data <- surveys_merged |>
  filter(Btot < 500000) |>
  filter(Ntot < 10000) |>
  filter (p_biom_nutflux > 0.8 ) |>
  filter (p_abund_nutflux > 0.8 ) 

nrow(task3_data) # 2 402 surveys
n_distinct(task3_data$SiteCode ) # 1 505 sites

# #Removing outliers, N and P recycling values superior to 99.9% of values
# task3_data_without_outliers <-  task3_data |>
#   filter( recycling_N < quantile(task3_data$recycling_N, 0.999) ) |>
#   filter( recycling_P < quantile(task3_data$recycling_P, 0.999) )


# log-transforming abundance and biomass
task3_data_surveys <- task3_data |>
  mutate(log10_biomass=log10(Btot_fishflux) ) |>
  mutate(log10_density=log10(Ntot_fishflux) )

summary(task3_data_surveys)

n_distinct(task3_data_surveys$SurveyID ) # 2399 unique surveys

save(task3_data_surveys, file=here::here("recycling", "outputs", "flux_final_data_surveys.Rdata")  )

list_surveys_recycling <- task3_data_surveys$SurveyID
#save(list_surveys_recycling, file=here::here("data", "list_surveys_recycling.Rdata")  )