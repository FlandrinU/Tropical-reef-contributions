#################################################################################
#'
#' This script means the surveys per site to obtain nutrient flows at the 
#'  community scale.
#' 
#'
#'@author Sebastien Villéger
#'        Ulysse Flandrin, \email{ulysse.flandrin@@gmail.com}
#'
#'
################################################################################

## cleaning memory
rm(list=ls())

# datasets from surveys
load(here::here("data","metadata_surveys.Rdata"))
load(here::here("data","data_species.Rdata"))
load(here::here("data","data_surveys.Rdata"))
load(here::here("recycling", "outputs", "surveys_fluxes.Rdata"))


# merging survey metadata with nutrient fluxes
surveys_merged <- metadata_surveys |> 
  left_join( surveys_fluxes , by="SurveyID")

summary(surveys_merged)


# sites metadata
sites_metadata <- surveys_merged |>
  select(SiteCode, SiteLatitude, SiteLongitude, SiteCountry, SiteEcoregion, SiteMeanSST) |>
  distinct(SiteCode, .keep_all = TRUE )

             
## number of transects per Site
sites_nbsurveys <- surveys_merged |>
  distinct(SiteCode, SurveyID) |>
  group_by(SiteCode) |>
  summarize(n_surveys=n())
nrow(sites_nbsurveys) # 2 011
summary(sites_nbsurveys$n_surveys) # Q1=1, median= 2, Q3=2, max=7
length(which(sites_nbsurveys$n_surveys==1)) # 627 sites with only one survey
length(which(sites_nbsurveys$n_surveys>=5)) # only 10 sites with at least 5 surveys

sites_metadata <- sites_metadata |>
  left_join(sites_nbsurveys)

summary(sites_metadata)



# averaging fluxes among surveys from each site
sites_fluxes <- surveys_merged |>
  group_by(SiteCode) |>
  select(SiteCode,
         Ntot, Btot, Ntot_fishflux, Btot_fishflux,
         storage_N, storage_P, 
         excretion_N, excretion_P, 
         egestion_N, egestion_P,
         recycling_N, recycling_P) |>
  summarise(across(.cols = everything(), .fns = mean, .names = "{.col}"))
  
  summary(sites_fluxes)
  

# computing quality metrics about fluxes estimates and ratio between fluxes
sites_fluxes <- sites_fluxes |>
  mutate(p_abund_nutflux = Ntot_fishflux / Ntot ) |>
  mutate(p_biom_nutflux = Btot_fishflux / Btot ) |>
  mutate(pexcr_recycling_N = excretion_N / recycling_N ) |>
  mutate(pexcr_recycling_P = excretion_P / recycling_P ) |> 
  mutate(recyc_stor_N = recycling_N / storage_N) |>
  mutate(recyc_stor_P = recycling_P / storage_P) |>
  mutate(log10_biomass=log10(Btot_fishflux) ) |>
  mutate(log10_density=log10(Ntot_fishflux) )



# size_structure of fishes (weighted by number of individuals)
sites_sizeQ<- data_surveys |>
  left_join( select(metadata_surveys, SurveyID, SiteCode) ) |>
  group_by(SiteCode) |>
  uncount( number ) |>
  summarize( size_Q1=quantile(size_class,0.25) ,
             size_median=quantile(size_class,0.5) ,
             size_Q3=quantile(size_class, 0.75),
             size_p90=quantile(size_class, 0.90)
  )
summary(sites_sizeQ)  


# trophic structure  (weighted by biomass of species)
sites_trophic<- data_surveys |>
  left_join( select(metadata_surveys, SurveyID, SiteCode) ) |>
  left_join( select(data_species, species, Diet) ) |>
  group_by(SiteCode) |>
  mutate(biom_dether = case_when(Diet=="herbivores_microvores_detritivores" ~ biomass,
                                 TRUE ~ 0) ) |>
  mutate(biom_plankt = case_when(Diet=="planktivores" ~ biomass,
                                 TRUE ~ 0) ) |>
  mutate(biom_pisci = case_when(Diet=="piscivores" ~ biomass,
                                TRUE ~ 0) ) |>
  summarize( across( starts_with("biom"), sum, .names = "tot_{.col}" ) ) |>
  mutate(pbiom_dether = tot_biom_dether / tot_biomass) |>
  mutate(pbiom_plankt = tot_biom_plankt / tot_biomass) |>
  mutate(pbiom_pisci = tot_biom_pisci / tot_biomass)

summary(sites_trophic)  




# merging metadata, fluxes and indices
sites_merged <- left_join(sites_metadata, sites_fluxes) |>
  left_join(sites_sizeQ) |>
  left_join(sites_trophic)



# exploratory plot about biomass vs number of transects ->uncomment to see the plot
nrow(sites_merged)
toplot <- sites_merged |>
  mutate(nb_surveys = case_when ( n_surveys == 1 ~ "1",
                                  n_surveys == 2 ~ "2",
                                  n_surveys >= 3 ~ "3-6" )  )
summary(as.factor(toplot$nb_surveys) )

# ggplot(toplot, aes(x = nb_surveys, y = log10_biomass, fill = nb_surveys ) ) +
#   geom_violin( draw_quantiles = c(0.5, 0.75, 0.95),  trim=TRUE, scale="count")  



# filtering as for surveys
task3_data_sites <- sites_merged |>
  filter(Btot < 500000) |>
  filter(Ntot < 10000) |>
  filter (p_biom_nutflux > 0.8 ) |>
  filter (p_abund_nutflux > 0.8 )  

n_distinct(task3_data_sites$SiteCode ) # 1 343 sites
list_sites_recycling <- task3_data_sites$SiteCode

# saving
save(list_sites_recycling, file = here::here("data", "list_sites_recycling.Rdata"))
save(task3_data_sites, file=here::here("recycling", "outputs", "flux_final_data_sites.Rdata")  )


