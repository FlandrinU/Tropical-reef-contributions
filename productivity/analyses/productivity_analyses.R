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
RLS_clean$Mmax = RLS_clean$lwa*RLS_clean$MaxLength*RLS_clean$lwb
RLS_clean$logMmax  = log(RLS_clean$Mmax)
RLS_clean$logLmax  = log(RLS_clean$MaxLength)

data_final <- RLS_clean  |> dplyr::mutate(Area=50*10)  |> dplyr::arrange(SurveyID)


#----------------Run functions to predict productivity------------------------
#Calculating productivity
#--------- individual level ---------#
RLS_prod_indiv = calc_prod_rfishprod(data_final)

save(RLS_prod_indiv, file = here::here("productivity","outputs", "RLS_prod_indiv.Rdata"))

#--------- transect level ---------#
print("Aggregating productivity prediction at transect level")
RLS_prod_transect = calc_prod_transect(RLS_prod_indiv,metadata_surveys)

save(RLS_prod_transect, file = here::here("productivity/outputs/RLS_prod_transect.Rdata"))

#--------- site level ---------#
print("Aggregating productivity prediction at site level")
RLS_prod_site = calc_prod_site(RLS_prod_indiv,metadata_surveys)

save(RLS_prod_site, file = here::here("productivity/outputs/RLS_prod_site.Rdata"))



#--------- species level ---------#
print("Aggregating productivity prediction at species level")
productivity_sp <- RLS_prod_indiv |> 
  dplyr::mutate(productivity = (Prod/Biom)*100) |> 
  dplyr::group_by(Species) |> 
  dplyr::mutate(mean_size = mean(Size),
         mean_Kmax = mean(Kmax),
         mean_productivity = mean(productivity)) |> 
  dplyr::ungroup() |> 
  dplyr::select(Species, mean_size, mean_Kmax, mean_productivity)

save(productivity_sp, file = here::here("productivity/outputs/RLS_prod_species.Rdata"))


#--------- assess sensitivity ---------#
#Sensitivity to account for variability between transects
RLS_prod_sensitivity = RLS_prod_indiv |> 
  dplyr::filter(Size < quantile(Size, 0.95)) |>  
  dplyr::filter(Num < quantile(Num, 0.95))
RLS_prod_sensitivity_transect = calc_prod_transect(RLS_prod_sensitivity, metadata_surveys)

RLS_prod_sensitivity_transect = RLS_prod_sensitivity_transect |> 
  dplyr::filter(Biom < quantile(Biom,0.99)) |> 
  dplyr::filter(Productivity < quantile(Productivity,0.99))

save(RLS_prod_sensitivity_transect, file = here::here("productivity/outputs/RLS_prod_sensitivity_transect.Rdata"))


RLS_prod_sensitivity_site = calc_prod_site(RLS_prod_sensitivity, metadata_surveys)

RLS_prod_sensitivity_site = RLS_prod_sensitivity_site |> 
  dplyr::filter(Biom < quantile(Biom,0.99)) |> 
  dplyr::filter(Productivity < quantile(Productivity,0.99))

save(RLS_prod_sensitivity_site, file = here::here("productivity/outputs/RLS_prod_sensitivity_site.Rdata"))


summary(RLS_prod_transect)
range(RLS_prod_transect$Productivity)



