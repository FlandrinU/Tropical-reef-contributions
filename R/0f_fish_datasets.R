################################################################################
#' 
#' Organize data, rename variables, ... and save all necessary information
#' into 3 files : - `data/metadata_surveys.Rdata`
#'                - `data/data_species.Rdata`
#'                - `data/data_surveys.Rdata`
#'                
#' Functions used in the scripts of the script `analysis/preping_data_rls.R`
#'
#' @author Sebastien Villéger
#' 
################################################################################
## cleaning memory
rm(list=ls())

## preparing environment ######

# loading data from RLS and from fishflux
load(here::here("data_raw","rls_trop_fish_sizeok_biomok.Rdata"))
head(rls_trop_fish_sizeok_biomok)
RLS_socio_data <- readRDS(here::here("data_raw", "source", "RLS_socio_withoutNA.rds"))

# checking taxonomic resolution
unique(rls_trop_fish_sizeok_biomok$taxo_level) # both species and genus
rls_trop_fish_sizeok_biomok |>
  dplyr::filter(taxo_level=="genus") |>
  dplyr::distinct(taxa) |>
  nrow() # 81 genera

## data about species for running fishflux ####

# loading
species_parameters<-read.csv(here::here("data_raw", "source", "species_parameters.csv") )
head(species_parameters)
nrow(species_parameters) # 1 110 species
sort(unique(species_parameters$family)) # 25 families (parrotfishes within Labridae)

# species names----
fishflux_species <- data.frame(species=species_parameters$species, fishflux="yes")

# merging fishflux status with data from RLS
data_rls_fishflux<-rls_trop_fish_sizeok_biomok |>
  dplyr::left_join( fishflux_species, by=c("taxa"="species")  ) |>
  dplyr::mutate(fishflux = tidyr::replace_na(fishflux, "no") )
summary(data_rls_fishflux)


## computing abundance, biomass metrics at survey level #####

# computing total number of individuals per survey ----
surveys_Ntot<-data_rls_fishflux |>
  dplyr::group_by(SurveyID) |>
  dplyr::select(SurveyID, number) |>
  dplyr::summarize(Ntot=sum(number) )
summary(surveys_Ntot$Ntot) # from 5 to 42 376 fish per survey


# computing total biomass ----
surveys_Btot<-data_rls_fishflux |>
  dplyr::group_by(SurveyID) |>
  dplyr::select(SurveyID, biomass) |>
  dplyr::summarize(Btot=sum(biomass) )
summary(surveys_Btot$Btot) # from 232 to 7 876 307g (7.8 tons) of fish


# computing total biomass of species with fishflux estimates ----
surveys_Btot_fishflux<-data_rls_fishflux |>
  dplyr::group_by(SurveyID) |>
  dplyr::filter(fishflux=="yes") |>
  dplyr::select(SurveyID, biomass) |>
  dplyr::summarize(Btot_fishflux=sum(biomass) )
summary(surveys_Btot_fishflux$Btot_fishflux) # from 9 to 6 062 096g of fish


# computing total abundance of species with fishflux estimates ----
surveys_Ntot_fishflux<-data_rls_fishflux |>
  dplyr::group_by(SurveyID) |>
  dplyr::filter(fishflux=="yes") |>
  dplyr::select(SurveyID, number) |>
  dplyr::summarize(Ntot_fishflux=sum(number) )
summary(surveys_Ntot_fishflux$Ntot_fishflux) # from  1 to 32 637 fish


# merging these 2 metrics in a single table and adding 2 quality metrics ----
# 'p_biom_nutflux' 'p_abund_nutflux' = ratio between total biomass or abundance  
# of species with nutrient cycling estimates and total biomass or abundance surveyed

surveys_quality <- surveys_Ntot |>
  dplyr::left_join(surveys_Btot, by="SurveyID") |>
  dplyr::left_join(surveys_Ntot_fishflux, by="SurveyID") |>
  dplyr::left_join(surveys_Btot_fishflux, by="SurveyID") |>
  tidyr::replace_na(list(Btot_fishflux = 0, Ntot_fishflux = 0) ) |>
  dplyr::mutate (p_abund_nutflux = Ntot_fishflux / Ntot) |>
  dplyr::mutate( p_biom_nutflux = Btot_fishflux / Btot )
summary(surveys_quality)
nrow(surveys_quality) # 3 631 surveys


# summary of % of biomass or abudance belonging to taxa in fishflux
summary(surveys_quality$p_biom_nutflux) # median = 0.930, Q3=0.98
summary(surveys_quality$p_abund_nutflux) # median = 0.966, Q3=0.99


# building 3 datasets: survey features, species info, fish obs ####
names(data_rls_fishflux)

# keeping only observation of species with bioenergetics models and actual abundance
data_rls_fishflux <- data_rls_fishflux |>
  dplyr::filter(number>0) |>
  dplyr::filter(fishflux=="yes")


# surveys metadata and quality metrics for nutrient recycling estimates ----
  # gravity data
  gravity <- RLS_socio_data |>
    dplyr::select( SurveyID, gravtot2, HDI, MarineEcosystemDependency) |>
    dplyr::mutate(SurveyID = as.character(SurveyID))


metadata_surveys <- data_rls_fishflux |>
  dplyr::select(SurveyID, SurveyDate, SurveyDepth, 
         SiteCode, SiteLatitude, SiteLongitude, SiteCountry, SiteEcoregion,
         SiteMeanSST, SiteMinSST, mpa_name, mpa_enforcement, protection_status,
         year_of_protection, gears_allowed, mpa_iucn_cat) |>
  dplyr::distinct(SurveyID, .keep_all = TRUE) |>
  dplyr::left_join(surveys_quality, by="SurveyID") |>
  dplyr::left_join( gravity)
  

# summary of sites surveyed
summary(metadata_surveys)
dplyr::n_distinct(metadata_surveys$SurveyID) # 3 628 surveys
dplyr::n_distinct(metadata_surveys$SiteCode) # 2 011 sites
dplyr::n_distinct(metadata_surveys$SiteCountry) # 39 countries
summary(metadata_surveys$SiteLatitude) # latitude from -32.4 to 29.5
summary(metadata_surveys$SiteLongitude) # longitude from -179 to 170


# taxonomy and traits of species ----
data_species<-data_rls_fishflux |>
  dplyr::select( species=taxa, family, 
          Size=MaxLength, Position, Activity ) |>
  dplyr::distinct(species,.keep_all = TRUE) |>
  dplyr::mutate(species=as.character(species))   |>
  dplyr::mutate(family=forcats::as_factor(family))

summary(data_species)
dplyr::n_distinct(data_species$species) # 1 024 taxa
unique(data_species$family) # from 26 families (25+Scaridae)

# recoding Position as an ordinal trait after merging the 2 types of Pelagic
# capital letter for all categories
data_species$Position<-forcats::fct_recode(data_species$Position,
                                  Benthic="benthic", 
                                  Pelagic="pelagic site attached", Pelagic="pelagic non-site attached")
data_species$Position<-factor(data_species$Position, 
                              levels=c("Benthic", "Demersal", "Pelagic"), ordered = TRUE)


# Diet categories from Parravicini et al 2020 (https://doi.org/10.1371/journal.pbio.3000702)
diet_pbiol <- species_parameters |>
  dplyr::select(species, diet_cat, a = lwa_m, b = lwb_m) |>
  dplyr::mutate(Diet=factor(diet_cat) ) |>
  dplyr::mutate(Diet= forcats::fct_recode(Diet,
                          sessile_invertivores="1",
                          herbivores_microvores_detritivores="2",
                          corallivores="3",
                          piscivores="4",
                          microinvertivores="5",
                          macroinvertivores="6",
                          crustacivores="7",
                          planktivores="8")
  ) |>
  dplyr::select(species, Diet, a, b)

summary(diet_pbiol$Diet)

# merging
data_species <- data_species |>
  dplyr::left_join(diet_pbiol)

summary(data_species)

 dplyr::filter(data_species, is.na(Position))
# => 2 NA for Naso_lopezi, replacing with value taken from FihsBase
 data_species[which(data_species$species=="Naso_lopezi"),"Position"]<-forcats::fct_explicit_na("Pelagic")
 data_species[which(data_species$species=="Naso_lopezi"),"Activity"]<-forcats::fct_explicit_na("Day")
 
 
# fish number, size and biomass in surveys ----
data_surveys<-data_rls_fishflux |>
  dplyr::select(SurveyID, species=taxa, size_class, number, biomass )
summary(data_surveys)


##----------------- Deal with different species names ----------------
#Extract spec_codes from FishBase - did in March 2023 - 
sp <- data_species$species

Species <- gsub("_", " ", sp)
Species_corrected <- rep("NA",length(sp))
SpecCode <- rep("NA",length(sp))

correctNames <- parallel::mclapply(1:length(Species), mc.cores = parallel::detectCores()-5, function(k){
  test <- rfishbase::validate_names(Species[k])
  if(length(test)==1){
    Species_corrected[k] <- test
    SpecCode[k] <- as.numeric((unique(rfishbase::species(test,fields="SpecCode"))))
  }else{
    next
  }
  c(Species_corrected[k], SpecCode[k])
})#end of k

species_spec_code <- as.data.frame(gsub(" ", "_", cbind(Species, do.call(rbind, correctNames))))
colnames(species_spec_code) <- c("species", "species_corrected", "spec_code")
species_spec_code$spec_code <- as.numeric(species_spec_code$spec_code)
head(species_spec_code)

save(species_spec_code, file = here::here("data", "species_spec_code.Rdata"))
# load(file = here::here("data", "species_spec_code.Rdata"))


##----------------- Merge dataset ----------------
data_species <- data_species |>
   dplyr::left_join(species_spec_code) |>
   dplyr::select(species, species_corrected, spec_code, family, Size, a, b, Position, Activity, Diet)

data_surveys <- data_surveys |>
  dplyr::left_join(species_spec_code) |>
  dplyr::select(SurveyID, species, spec_code, size_class, number, biomass) 

## saving as Rdata ####
save( metadata_surveys, file=here::here("data", "metadata_surveys.Rdata")  )
save( data_species, file=here::here("data", "data_species.Rdata")  )
save( data_surveys, file=here::here("data", "data_surveys.Rdata")  )

