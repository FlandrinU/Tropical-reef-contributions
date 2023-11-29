#################################################################################
#'
#'This script runs all scripts of `productivity/R/` to obtain the production and
#' productivity of biomass per day for each individuals and at the survey scale
#' 
#'
#'@author Raphael Seguin, 
#'        Nicolas Loiseau, \email{nicolas.loiseau1@@gmail.com}
#'
#'
################################################################################

#-----------------cleaning memory-------------------
rm(list=ls())
cat("Run productivity analysis... \n")

#-----------------Loading all data---------------------
load(here::here("data", "metadata_surveys.Rdata"))
load(here::here("data", "data_species.Rdata"))
load(here::here("data", "data_surveys.Rdata"))

#-----------------Loading functions---------------------
path = (here::here("productivity", "R"))
setwd(path)
sapply(list.files(path), source)

#----------------Prepping RLS data------------------------
#rename variables
surveys <- dplyr::select(data_surveys, SurveyID, Species = species, Num = number, 
                         Sizeclass = size_class, Biomass = biomass)

sp <- dplyr::select(data_species,Species = species,Family = family,
                    MaxLength = Size, lwa = a, lwb = b)

metadata <- dplyr::select(metadata_surveys, SurveyID, Temperature = SiteMeanSST)

#Merge variables
RLS_clean <- surveys |> 
  dplyr::left_join(sp) |> 
  dplyr::left_join(metadata)

#add new variables
RLS_clean$Mmax = RLS_clean$lwa*RLS_clean$MaxLength^RLS_clean$lwb
RLS_clean$logMmax  = log(RLS_clean$Mmax)
RLS_clean$logLmax  = log(RLS_clean$MaxLength)

data_final <- RLS_clean  |> dplyr::mutate(Area=50*10)  |> dplyr::arrange(SurveyID)


#----------------Filtering only fishable fish -> productivity of available biomass ------------------------
## keep only targeted families, according to Cinner et al. 2020 and expert opinion##
#families in which all species are targeted in the family: (Cinner et al. 2020 + other expert opinion )
all_sp <- c("Acanthuridae", "Caesionidae", "Carangidae", "Ephippidae", "Haemulidae", "Kyphosidae",
            "Labridae", "Lethrinidae", "Lutjanidae", "Mullidae", "Nemipteridae", "Scaridae",
            "Scombridae", "Serranidae", "Siganidae", "Sparidae", "Sphyraenidae",
            "Mugilidae", "Bothidae", "Scorpaenidae", "Sciaenidae")

#families with species targeted if larger than 20cm: (Cinner et al. 2020 + other expert opinion )
target_larger_20cm <- c("Balistidae", "Holocentridae", "Pomacanthidae", "Priacanthidae")

#non-targeted families: (Cinner et al. 2020 + other expert opinion )
untargeted_fam <- c("Chaetodontidae","Cirrhitidae", "Diodontidae", "Grammatidae", 
                    "Monacanthidae", "Pempheridae", "Pinguipedidae", "Pseudochromidae", 
                    "Synodontidae", "Tetraodontidae", "Zanclidae","Fistulariidae",
                    "Ostraciidae", "Pomacentridae")

unclassified_fam <- setdiff(unique(data_final$family), c(all_sp, target_larger_20cm, untargeted_fam))
unclassified_fam #OK

#keep only targeted fishes
data_final_fishery <- data_final |> 
  dplyr::filter(Family %in% all_sp | (Family %in% target_larger_20cm & MaxLength > 20))

unique(data_final_fishery$Family)


#----------------Run functions to predict productivity------------------------
#Calculating productivity
#--------- individual level ---------#
RLS_prod_indiv = calc_prod_rfishprod(data_final_fishery)

save(RLS_prod_indiv, file = here::here("productivity","outputs", "RLS_prod_indiv.Rdata"))

#--------- transect level ---------#
print("Aggregating productivity prediction at transect level")
RLS_prod_transect = calc_prod_transect(RLS_prod_indiv,metadata_surveys)

save(RLS_prod_transect, file = here::here("productivity/outputs/RLS_prod_transect.Rdata"))
