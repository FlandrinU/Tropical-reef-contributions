#################################################################################
#'
#' This script summarize nutrients flows of each individuals at the survey 
#'  scale (nutrient flows in each surveys)
#' 
#'
#'@author Sebastien Villéger
#'
#'
################################################################################

## cleaning memory
rm(list=ls())

## loading  datasets
load(here::here("data", "metadata_surveys.Rdata"))
load(here::here("data", "data_species.Rdata"))
load(here::here("data", "data_surveys.Rdata"))

# details about species: length-weight, age, nutrient contents, diet cat
species_par <- readr::read_csv(here::here("data", "parameters_sp_sst.csv")) |>
  dplyr::select(family, species, lwa_m, lwb_m, k_m, linf_m, sst, Qc_m, Qn_m, Qp_m, diet_cat) 



# fluxes from fishes
data_fishflux <- readr::read_csv(here::here("recycling", "outputs", "cnpflux_sp_size_sst.csv")) |> unique()
head(data_fishflux)
names(data_fishflux)

# merging datasets and computing fluxes for each species*size_class given abundance and sst
data_surveys_fluxes<- dplyr::left_join(data_surveys, 
                                dplyr::select(metadata_surveys, SurveyID, sst = SiteMeanSST ),
                                by="SurveyID") |>
  dplyr::mutate(sst = round(sst)) |>
  dplyr::left_join(data_fishflux) |>
  dplyr::left_join(unique(species_par))  |>
  dplyr::mutate(fishflux = as.factor(dplyr::case_when(is.na(Gc_median) ~ FALSE, 
                                        Gc_median >= 0 ~ TRUE) )
  ) |>
  dplyr::mutate(storage_C = Gc_median * number,
                storage_N = Gn_median * number,
                storage_P = Gp_median * number,
                excretion_N = Fn_median * number,
                excretion_P = Fp_median * number,
                egestion_N = Wn_median * number,
                egestion_P = Wp_median * number,
                egestion_C = Wc_median * number,
                biomass_estimated = lwa_m * (size_class^lwb_m) * number)

head(data_surveys_fluxes)
# summary(data_surveys_fluxes)
nrow(data_surveys_fluxes)

# computing for each flux, matrix survey*species (as for biomass) ----

# variables of interest with new names
fluxes_var<-c("storage_C", "storage_N", "storage_P",
              "excretion_N", "excretion_P", 
              "egestion_C", "egestion_N", "egestion_P",
              "biomass_estimated"
              )

# list to store matrices
surveys_species_fluxes<-list()

# loop on variables
for (k in fluxes_var ) {
  
  # merging size_classes per species in each survey
  mat_k<-data_surveys_fluxes |> 
    dplyr::select(SurveyID, species, !!k ) |>
    dplyr::rename( flux_k = !!k ) |>
    dplyr::group_by(SurveyID, species) |>
    dplyr::summarize( sum_flux_k=sum(flux_k) ) |>
    tidyr::pivot_wider(names_from = species, values_from = sum_flux_k, values_fill = 0) |>
    tibble::column_to_rownames(var="SurveyID") |>
    as.matrix()  
  
  # storing
  surveys_species_fluxes[[k]]<-mat_k
  
}# end of k

# lapply(surveys_species_fluxes, dim)

# total of fluxes per survey ----
surveys_fluxes<-unlist( sapply(surveys_species_fluxes, rowSums) ) |> 
  as.data.frame() |>
  tibble::rownames_to_column("SurveyID")
head(surveys_fluxes)
dim(surveys_fluxes)

## computing stoichiometric ratios of excretion and egestion
# for each species*size then median at survey level
mC<-12
mN<-14
mP<-31
median_ratio <- data_surveys_fluxes |>
  tidyr::uncount(number) |>
  dplyr::group_by(SurveyID) |>
  dplyr::summarise(egestion_CN = median( (egestion_C/mC) / (egestion_N/mN) ),
            egestion_CP = median( (egestion_C/mC)/ (egestion_P/mP) ),
            egestion_NP = median( (egestion_N/mN) / (egestion_P/mP) ),
            excretion_NP = median( (excretion_N/mN) / (excretion_P/mP) )
            )

# merging and ratio of quality of excretion and of egestion
surveys_fluxes <- surveys_fluxes |>
  dplyr::left_join(median_ratio) |>
  dplyr::mutate(excrNP_egesNP= excretion_NP / egestion_NP)


## computing recycling as sum of excretion plus egestion
surveys_fluxes <- surveys_fluxes |>
  dplyr::mutate(recycling_C = egestion_C) |>
  dplyr::mutate(recycling_N = excretion_N + egestion_N) |>
  dplyr::mutate(recycling_P = excretion_P + egestion_P)
  
surveys_species_fluxes[["recycling_N"]] <- surveys_species_fluxes[["excretion_N"]] + surveys_species_fluxes[["egestion_N"]]
surveys_species_fluxes[["recycling_P"]] <- surveys_species_fluxes[["excretion_P"]] + surveys_species_fluxes[["egestion_P"]]


## computing contribution of excretion to recycling of N and P
surveys_fluxes <- surveys_fluxes |>
  dplyr::mutate(pexcr_recycling_N = excretion_N / recycling_N) |>
  dplyr::mutate(pexcr_recycling_P = excretion_P / recycling_P)
  

## computing ratio between recycling and storing
surveys_fluxes <- surveys_fluxes |>
  dplyr::mutate( recyc_stor_C = recycling_C / storage_C ) |>
  dplyr::mutate( recyc_stor_N = recycling_N / storage_N ) |>
  dplyr::mutate( recyc_stor_P = recycling_P / storage_P )
  

# summary
summary(surveys_fluxes)
nrow(surveys_fluxes)


# saving
save(surveys_species_fluxes, file=here::here("recycling", "outputs", "surveys_species_fluxes.Rdata") )
save(surveys_fluxes, file=here::here("recycling", "outputs", "surveys_fluxes.Rdata") )
