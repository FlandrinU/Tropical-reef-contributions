################################################################################
##
## Assess the number of elasmobranch and IUCN species in surveys
##
## iucn_and_elasmobranch_indices.R
##
## 09/11/2022
##
## Ulysse Flandrin
##
################################################################################
## cleaning memory
rm(list=ls())

##-------------loading data-------------
#datasets
load(here::here("data","data_species.Rdata"))
load(here::here("data","data_surveys.Rdata"))

load(here::here("data","data_species_elasmobranchii.Rdata"))
load(here::here("data","data_surveys_elasmobranchii.Rdata"))

load(here::here("biodiversity", "data","iucn_category_all_species_NLoiseau.Rdata"))

#IUCN key
IUCN_KEY <- "78d79dce9a7f762ddc88584f7b2b75e81170e9e6869d9dce3b1fe8fbb3b7f8f6" #Paste your IUCN key

## -------------IUCN index-------------

## merge data
all_sp <- rbind( dplyr::select(data_species, species, species_corrected),
                 dplyr::select(data_species_elasmo, species, species_corrected))

names <- questionr::na.rm(all_sp$species_corrected) # 3 elasmobranch species identifies only at the genus level

## IUCN redlist data
iucn_name <- gsub("_", " ", names)
iucn_data_raw <- parallel::mclapply(iucn_name, mc.cores = parallel::detectCores()-5,
                          function(i){rredlist ::rl_search(name = i,  key= IUCN_KEY)}) #/!\ long time to run
iucn_data <- do.call(rbind, lapply(iucn_data_raw, "[[", "result"))

table(iucn_data$category) # 3 CR, 35 DD, 10 EN, 819 LC, 18 NT, 34 VU

save(iucn_data, file = here::here("biodiversity", "outputs", "iucn_data_all_species.Rdata"))
# load(file = here::here("biodiversity", "outputs", "iucn_data_all_species.Rdata"))

#keep the iucn category
iucn_category <- iucn_data |>
  dplyr::mutate(scientific_name = gsub(" ", "_", iucn_data$scientific_name)) |>
  dplyr::select( species_corrected = scientific_name, category) |>
  dplyr::right_join(all_sp)

table(iucn_category$category, useNA = "always") # 3 CR, 35 DD, 10 EN, 818 LC, 18 NT, 34 VU, 166 <NA>


#Complete category by random forrest and deep learning (cf N.Loiseau)
dat_network$IUCN_final <- as.character(dat_network$IUCN_final)

for(i in which(is.na(iucn_category$category))){
  name <- iucn_category$species[i]
  name_corrected <- iucn_category$species_corrected[i]
  
  if(length(which(dat_network$species == name))==1){
    iucn_category$category[i] <- dat_network$IUCN_final[which(dat_network$species == name)] 
  }else{ if(length(which(dat_network$species == name_corrected))==1){
    iucn_category$category[i] <- dat_network$IUCN_final[which(dat_network$species == name_corrected)]
  }}
}

table(iucn_category$category, useNA = "always") 
# 3 CR, 34 DD,  10 EN, 755 LC,  6 No Status, 202 Non Threatened, 18 NT, 17 Threatened, 32 VU, 7 NA



#merge survey data and iucn category
all_surveys_iucn <- rbind( dplyr::select(data_surveys, SurveyID, species, size_class, number, biomass),
                           dplyr::select(data_surveys_elasmo, SurveyID, species, size_class, number, biomass)) |>
  dplyr::left_join(iucn_category) 
  

#remove species without iucn category
all_surveys_iucn <- all_surveys_iucn |>
  dplyr::filter(is.na(all_surveys_iucn$category) == F,
                category != "No Status")

## Number of IUCN species per surveys
iucn_category_surveys <- all_surveys_iucn |>
  dplyr::mutate(category = forcats::fct_recode(category, 
                                        "0" = "LC",
                                        "0" = "NE",
                                        "0" = "DD",
                                        "0" = "NT",
                                        "0" = "Non Threatened",
                                        "1" = "VU",
                                        "1" = "EN",
                                        "1" = "CR",
                                        "1" = "Threatened")) 


iucn_by_surveys <- iucn_category_surveys|>
  dplyr::mutate(category = as.numeric(as.character(category))) |>
  dplyr::group_by(SurveyID) |>
  dplyr::distinct(species, .keep_all=TRUE) |> #keep only one size class per species
  dplyr::summarise(iucn_species = sum(category))
  
table(iucn_by_surveys$iucn_species)
      #nb of iucn species:   0    1    2    3    4    5    6 
      # nb of surveys     2244  979  304   78   18    3    1  



## -------------Elasmobranch index-------------
elasmo_by_surveys <- data_surveys_elasmo |>
  dplyr::mutate( number = 1) |>
  dplyr::group_by(SurveyID) |>
  dplyr::distinct(species, .keep_all=TRUE) |> #keep only one size class per species
  dplyr::summarise(elasmobranch_diversity = sum(number))
  
table(elasmo_by_surveys$elasmobranch_diversity)
      # nb of elasmobranch species: 1   2   3   4 
      # nb of surveys             576  92  17   1 

##merge indices
iucn_and_elasmo_by_surveys <- dplyr::left_join(iucn_by_surveys, elasmo_by_surveys) |>
  dplyr::mutate( elasmobranch_diversity = tidyr::replace_na(elasmobranch_diversity, 0) )



## -------------save iucn and chondrichtyans data-------------
save(iucn_and_elasmo_by_surveys, 
     file = here::here("biodiversity", "outputs", "surveys_iucn_and_chondrichtyans_scores.Rdata"))
