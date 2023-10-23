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
  dplyr::left_join( surveys_fluxes , by="SurveyID")

summary(surveys_merged)


# sites metadata
sites_metadata <- surveys_merged |>
  dplyr::select(SiteCode, SiteLatitude, SiteLongitude, SiteCountry, SiteEcoregion, SiteMeanSST) |>
  dplyr::distinct(SiteCode, .keep_all = TRUE )


## number of transects per Site
sites_nbsurveys <- surveys_merged |>
  dplyr::distinct(SiteCode, SurveyID) |>
  dplyr::group_by(SiteCode) |>
  dplyr::summarize(n_surveys=dplyr::n())
nrow(sites_nbsurveys) # 2 011
summary(sites_nbsurveys$n_surveys) # Q1=1, median= 2, Q3=2, max=7
length(which(sites_nbsurveys$n_surveys==1)) # 627 sites with only one survey
length(which(sites_nbsurveys$n_surveys>=5)) # only 10 sites with at least 5 surveys

sites_metadata <- sites_metadata |>
  dplyr::left_join(sites_nbsurveys)

summary(sites_metadata)



# averaging fluxes among surveys from each site
sites_fluxes <- surveys_merged |>
  dplyr::group_by(SiteCode) |>
  dplyr::select(SiteCode,
                Ntot, Btot, Ntot_fishflux, Btot_fishflux,
                storage_N, storage_P, 
                excretion_N, excretion_P, 
                egestion_N, egestion_P,
                recycling_N, recycling_P) |>
  dplyr::summarise(across(.cols = everything(), .fns = mean, .names = "{.col}"))

summary(sites_fluxes)


# computing quality metrics about fluxes estimates and ratio between fluxes
sites_fluxes <- sites_fluxes |>
  dplyr::mutate(p_abund_nutflux = Ntot_fishflux / Ntot ) |>
  dplyr::mutate(p_biom_nutflux = Btot_fishflux / Btot ) |>
  dplyr::mutate(pexcr_recycling_N = excretion_N / recycling_N ) |>
  dplyr::mutate(pexcr_recycling_P = excretion_P / recycling_P ) |> 
  dplyr::mutate(recyc_stor_N = recycling_N / storage_N) |>
  dplyr::mutate(recyc_stor_P = recycling_P / storage_P) |>
  dplyr::mutate(log10_biomass=log10(Btot_fishflux) ) |>
  dplyr::mutate(log10_density=log10(Ntot_fishflux) )



# size_structure of fishes (weighted by number of individuals)
sites_sizeQ<- data_surveys |>
  dplyr::left_join( dplyr::select(metadata_surveys, SurveyID, SiteCode) ) |>
  dplyr::group_by(SiteCode) |>
  tidyr::uncount( number ) |>
  dplyr::summarize( size_Q1=quantile(size_class,0.25) ,
                    size_median=quantile(size_class,0.5) ,
                    size_Q3=quantile(size_class, 0.75),
                    size_p90=quantile(size_class, 0.90)
  )
summary(sites_sizeQ)  


# trophic structure  (weighted by biomass of species)
sites_trophic<- data_surveys |>
  dplyr::left_join( dplyr::select(metadata_surveys, SurveyID, SiteCode) ) |>
  dplyr::left_join( dplyr::select(data_species, species, Diet) ) |>
  dplyr::group_by(SiteCode) |>
  dplyr::mutate(biom_dether = dplyr::case_when(Diet=="herbivores_microvores_detritivores" ~ biomass,
                                        TRUE ~ 0) ) |>
  dplyr::mutate(biom_plankt = dplyr::case_when(Diet=="planktivores" ~ biomass,
                                        TRUE ~ 0) ) |>
  dplyr::mutate(biom_pisci = dplyr::case_when(Diet=="piscivores" ~ biomass,
                                       TRUE ~ 0) ) |>
  dplyr::summarize( dplyr::across( starts_with("biom"), sum, .names = "tot_{.col}" ) ) |>
  dplyr::mutate(pbiom_dether = tot_biom_dether / tot_biomass) |>
  dplyr::mutate(pbiom_plankt = tot_biom_plankt / tot_biomass) |>
  dplyr::mutate(pbiom_pisci = tot_biom_pisci / tot_biomass)

summary(sites_trophic)  




# merging metadata, fluxes and indices
sites_merged <- dplyr::left_join(sites_metadata, sites_fluxes) |>
  dplyr::left_join(sites_sizeQ) |>
  dplyr::left_join(sites_trophic)



# exploratory plot about biomass vs number of transects ->uncomment to see the plot
nrow(sites_merged)
toplot <- sites_merged |>
  dplyr::mutate(nb_surveys = dplyr::case_when ( n_surveys == 1 ~ "1",
                                         n_surveys == 2 ~ "2",
                                         n_surveys >= 3 ~ "3-6" )  )
summary(as.factor(toplot$nb_surveys) )

# ggplot(toplot, aes(x = nb_surveys, y = log10_biomass, fill = nb_surveys ) ) +
#   geom_violin( draw_quantiles = c(0.5, 0.75, 0.95),  trim=TRUE, scale="count")  



# filtering as for surveys
task3_data_sites <- sites_merged |>
  dplyr::filter(Btot < 500000) |>
  dplyr::filter(Ntot < 10000) |>
  dplyr::filter (p_biom_nutflux > 0.8 ) |>
  dplyr::filter (p_abund_nutflux > 0.8 )  

dplyr::n_distinct(task3_data_sites$SiteCode ) # 1 343 sites
list_sites_recycling <- task3_data_sites$SiteCode

# saving
save(list_sites_recycling, file = here::here("data", "list_sites_recycling.Rdata"))
save(task3_data_sites, file=here::here("recycling", "outputs", "flux_final_data_sites.Rdata")  )


