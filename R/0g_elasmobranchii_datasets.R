################################################################################################
#' Elasmobranch informations
#'
#'This script Extract elasmobanchii observations from raw rls data and keep only 
#' surveys of metadata_surveys.Rdata
#'
#'@author Ulysse Flandrin, \email{ulysse.flandrin@@gmail.com}
#'
#'
#' @date 2022/11/09 first created
################################################################################################

## cleaning memory
rm(list=ls())

## loading data from RLS ##
load(here::here("data", "metadata_surveys.Rdata"))
load(here::here("data_raw", "rls_all_raw.Rdata"))
head(rls_all_raw)


## filtering to keep only transects from tropical sites #####
# i.e with min(SST)>=17°C

rls_trop<-rls_all_raw |> 
  dplyr::filter( SiteMinSST >= 17)

dplyr::n_distinct(rls_trop$SurveyID) # 3 794 surveys
dplyr::n_distinct(rls_trop$SiteCode) # 2 057 sites
dplyr::n_distinct(rls_trop$taxa) # 2 321 taxa
dplyr::n_distinct(rls_trop$family) # 123 families


##  keeping only elasmobranchii ##
# looking at names of NO family and NO class
rls_trop |>
  dplyr::filter(class=="" | family=="") |> 
  dplyr::distinct(taxa, taxo_level, family, class)
# 83 taxa => check, only Actinopterygii except "Dasyatid"  ###### @@@@@ some xxid spp from families not removed in next scripts (e.g. Acanthurid spp, Scarid Scarus spp)

# filtering to remove all Actinopterygii #
rls_trop_elasmo<- rls_trop |>
  dplyr::filter(class=="Elasmobranchii")

head(rls_trop_elasmo)

# diversity remaining
dplyr::n_distinct(rls_trop_elasmo$SurveyID) # 730 surveys
dplyr::n_distinct(rls_trop_elasmo$family) # 17 families
dplyr::n_distinct(rls_trop_elasmo$taxa) # 60 taxa


names <- sort(rls_trop_elasmo$taxa) 
unique(names)
names[stringr::str_detect(names, pattern="spp.")]



#### species dataset ####
data_species_elasmo<-rls_trop_elasmo |>
  dplyr::select( species=taxa, family, 
          Size=MaxLength, Position, Activity ) |>
  dplyr::distinct(species,.keep_all = TRUE) |>
  dplyr::mutate(species=as.character(species))   |>
  dplyr::mutate(family=forcats::as_factor(family))

summary(data_species_elasmo)
dplyr::n_distinct(data_species_elasmo$species) # 60 taxa
unique(data_species_elasmo$family) # 17 families

# recoding Position as an ordinal trait after merging the 2 types of Pelagic
# capital letter for all categories
data_species_elasmo$Position <- forcats::fct_recode(data_species_elasmo$Position,
                                  Benthic="benthic", 
                                  Pelagic="pelagic site attached", Pelagic="pelagic non-site attached")
data_species_elasmo$Position<-factor(data_species_elasmo$Position, 
                              levels=c("Benthic", "Demersal", "Pelagic"), ordered = TRUE)


#### surveys dataset ####
data_surveys_elasmo <- rls_trop_elasmo |>
  dplyr::select(SurveyID, species=taxa, size_class, number, biomass )

summary(data_surveys_elasmo) #1130 surveys data
length(unique(data_surveys_elasmo$SurveyID)) #730 surveys

data_surveys_elasmo <- data_surveys_elasmo |>
  dplyr::filter(SurveyID %in% metadata_surveys$SurveyID)

#elasmobranchii data remaining
summary(data_surveys_elasmo) #1038 surveys data 
length(unique(data_surveys_elasmo$SurveyID)) #686 surveys (44 surveys absent in actinopterygii dataset)



#check scientific names
names <- data_species_elasmo$species
corrected_name <- rfishbase::validate_names( stringr::str_replace(names, "_", " ") )
corrected_name <- stringr::str_replace(corrected_name, " ", "_")  

#order data
data_species_elasmo <- cbind(data_species_elasmo, species_corrected = corrected_name) |>
  dplyr::mutate(class = "elasmobranchii") |>
  dplyr::select(species, species_corrected, family, class, Size, Position, Activity)
 
## saving as Rdata ##
save( data_species_elasmo, file=here::here("data", "data_species_elasmobranchii.Rdata")  )
save( data_surveys_elasmo, file=here::here("data", "data_surveys_elasmobranchii.Rdata")  )

