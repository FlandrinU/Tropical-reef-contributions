################################################################################
#' 
#' Filling missing value for fish biomass in surveys
#' Functions used in the scripts of the script `analysis/preping_data_rls.R`
#'
#' @author Sebastien Villéger
#' 
#' 
################################################################################
## preparing environment ######

#Fixing error in curl package for internet conexion
np = !is.null(curl::nslookup("r-project.org", error = FALSE))
assign("has_internet_via_proxy", np, environment(curl::has_internet))

# cleaning memory
rm(list=ls())

## loading data from RLS and from fishflux ####
load(here::here("data_raw", "rls_trop_fish_sizeok.Rdata"))
head(rls_trop_fish_sizeok)

## looking at size and biomass data ####

# missing biomass but size known
biom_na<-rls_trop_fish_sizeok |>
  dplyr::filter( is.na(biomass) ) 
nrow(biom_na) # 1 320 cases
summary(biom_na$number) # from 1 to 400 individuals
taxa_biomna<-unique(biom_na$taxa) 
length(taxa_biomna) # 56 taxa
sp_biomna<-taxa_biomna[! stringr::str_detect(taxa_biomna,"spp")]
length(sp_biomna) # 50 species (= not Genus spp)

## => missing biomass not an issue for estimating nutrient fluxes because model based on length
# => minor issue just to compute proportion of biomass belonging to species with fishflux estimates

## estimating missing based on length-weight allometric relationships ####

# retrieving parameters of length_weight relationship from Fishbase  ----
# thanks to fishflux, takes few minutes
lw_sp_biomna<-list()
for (k in sp_biomna) {
  lw_sp_biomna[[k]]<-tryCatch( fishflux::find_lw(sub("_"," ",k), mirror = "se"), 
                               error = function(err) "NA" )
  print(k)
} # end of k


lw_sp_biomna_df<-do.call(rbind, lw_sp_biomna)

# checking
lw_sp_biomna_df |> dplyr::filter(lwa_m=="NA") 
# => NA for 1 species Alectis indicus => unaccepted name
lw_sp_biomna_df["Alectis_indicus",] <- tryCatch( fishflux::find_lw(sub("_"," ","Alectis_indica"), mirror = "se"), 
          error = function(err) "NA" )

# computing biomass for observations of each species (number*a*size^b) ----
for (k in sp_biomna) {
  pos_k <- which(rls_trop_fish_sizeok$taxa==k)
  if (lw_sp_biomna_df[k,"lwa_m"]!="NA" ) {
    rls_trop_fish_sizeok[pos_k, "biomass"] <- rls_trop_fish_sizeok[pos_k,"number"] *
      as.numeric(lw_sp_biomna_df[k,"lwa_m"])  * 
      (rls_trop_fish_sizeok[pos_k, "size_class"]^as.numeric(lw_sp_biomna_df[k,"lwb_m"]) ) 
  }# end of if computable
}# end of k

# remaining missing values ----
biom_na<-rls_trop_fish_sizeok |>
  dplyr::filter( is.na(biomass) ) 
nrow(biom_na) # 8 cases
dplyr::n_distinct(biom_na$SurveyID) # 9 surveys
summary(biom_na$number) # from 1 to 23 individuals

# removing transects with missing fish biomass  ----
rls_trop_fish_sizeok_biomok <- rls_trop_fish_sizeok |>
  dplyr::filter ( !(SurveyID %in% unique(biom_na$SurveyID) ) )


# diversity remaining
dplyr::n_distinct(rls_trop_fish_sizeok_biomok$SurveyID) # 3 632 surveys
dplyr::n_distinct(rls_trop_fish_sizeok_biomok$family) # 81 families
dplyr::n_distinct(rls_trop_fish_sizeok_biomok$taxa) # 1 679 taxa

## saving as Rdata ####
save( rls_trop_fish_sizeok_biomok, file=here::here("data_raw", "rls_trop_fish_sizeok_biomok.Rdata")  )
