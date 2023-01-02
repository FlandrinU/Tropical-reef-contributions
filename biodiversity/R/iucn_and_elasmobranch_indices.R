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

#IUCN key
IUCN_KEY <- "78d79dce9a7f762ddc88584f7b2b75e81170e9e6869d9dce3b1fe8fbb3b7f8f6"

## -------------IUCN index-------------

## merge data
all_sp <- rbind( select(data_species, species, species_corrected),
                 select(data_species_elasmo, species, species_corrected))

names <- na.rm(all_sp$species_corrected) # 3 elasmobranch species identifies only at the genus level

## IUCN redlist data
iucn_name <- gsub("_", " ", names)
iucn_data_raw <- mclapply(iucn_name, mc.cores = 8,
                          function(i){rredlist ::rl_search(name = i,  key= IUCN_KEY)}) #long time to run
iucn_data <- do.call(rbind, lapply(iucn_data_raw, "[[", "result"))

table(iucn_data$category) # 3 CR, 35 DD, 10 EN, 819 LC, 18 NT, 34 VU

save(iucn_data, file = here::here("biodiversity", "outputs", "iucn_data_all_species.Rdata"))

#keep the iucn category
iucn_category <- iucn_data %>% 
  mutate(scientific_name = gsub(" ", "_", iucn_data$scientific_name)) %>%
  select( species_corrected = scientific_name, category) %>%
  right_join(all_sp)

table(iucn_category$category, useNA = "always") # 3 CR, 35 DD, 10 EN, 819 LC, 18 NT, 34 VU, 165 <NA>

#merge survey data and iucn category
all_surveys_iucn <- rbind( select(data_surveys, SurveyID, species, size_class, number, biomass),
                           select(data_surveys_elasmo, SurveyID, species, size_class, number, biomass)) %>%
  left_join(iucn_category) 
  
table(all_surveys_iucn$category, useNA = "always") 
    #   CR     DD     EN     LC     NT     VU   <NA> 
    #   54   1286    457 167390   1958   2094  30579 

#remove species without iucn category
all_surveys_iucn <- all_surveys_iucn %>% filter(is.na(all_surveys_iucn$category) == F)

## Number of IUCN species per surveys
iucn_category_surveys <- all_surveys_iucn %>%
  mutate(category = forcats::fct_recode(category, 
                                        "0" = "LC",
                                        "0" = "NE",
                                        "0" = "DD",
                                        "0" = "NA",
                                        "0" = "NT",
                                        "1" = "VU",
                                        "1" = "EN",
                                        "1" = "CR")) 

iucn_by_surveys <- iucn_category_surveys%>%
  mutate(category = as.numeric(as.character(category))) %>%
  group_by(SurveyID) %>%
  distinct(species, .keep_all=TRUE) %>% #keep only one size class per species
  summarise(iucn_species = sum(category))
  
table(iucn_by_surveys$iucn_species)
      #nb of iucn species:   0    1    2    3    4    5    6 
      # nb of surveys     2427  874  230   70   17    3    1 



## -------------Elasmobranch index-------------
elasmo_by_surveys <- data_surveys_elasmo %>%
  mutate( number = 1) %>%
  group_by(SurveyID) %>%
  distinct(species, .keep_all=TRUE) %>% #keep only one size class per species
  summarise(elasmobranch_diversity = sum(number))
  
table(elasmo_by_surveys$elasmobranch_diversity)
      # nb of elasmobranch species: 1   2   3   4 
      # nb of surveys             576  92  17   1 

##merge indices
iucn_and_elasmo_by_surveys <- left_join(iucn_by_surveys, elasmo_by_surveys) %>%
  mutate( elasmobranch_diversity = replace_na(elasmobranch_diversity, 0) )



## -------------save iucn and chondrichtyans data-------------
save(iucn_and_elasmo_by_surveys, 
     file = here::here("biodiversity", "outputs", "surveys_iucn_and_chondrichtyans_scores.Rdata"))
