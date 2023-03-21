#-----------------cleaning memory-------------------
rm(list=ls())

#-----------------Loading packages-------------------
if (!("rfishprod" %in% utils::installed.packages())) devtools::install_github("renatoamorais/rfishprod")

pkgs <- c("here", "dplyr", "stringr", "rfishprod", "stats")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

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
RLS_clean <- surveys %>%
  dplyr::left_join(sp) %>%
  dplyr::left_join(metadata)

#add new variables
RLS_clean$Mmax = RLS_clean$lwa*RLS_clean$MaxLength*RLS_clean$lwb
RLS_clean$logMmax  = log(RLS_clean$Mmax)
RLS_clean$logLmax  = log(RLS_clean$MaxLength)

data_final <- RLS_clean %>% mutate(Area=50*10) %>% arrange(SurveyID)


#----------------Run functions  to predict productivity------------------------
setwd(here::here())

#Calculating productivity
#--------- individual level ---------#
RLS_prod_indiv = calc_prod_rfishprod(data_final)

save(RLS_prod_indiv, file = "productivity/outputs/RLS_prod_indiv.Rdata")

#--------- transect level ---------#
print("Aggregating at transect level")
RLS_prod_transect = calc_prod_transect(RLS_prod_indiv,metadata_surveys)

# save(RLS_prod_transect,file = "productivity/outputs/RLS_prod_transect_with_outliers0.001.Rdata")
# 
# #Removing outliers, Biomass and productivity values superior to 99.9% of values
# RLS_prod_transect = RLS_prod_transect %>% filter(Biom < quantile(RLS_prod_transect$Biom,0.999)) %>% 
#                                           filter(Productivity < quantile(RLS_prod_transect$Productivity,0.999))

save(RLS_prod_transect, file = "productivity/outputs/RLS_prod_transect.Rdata")

#--------- site level ---------#
print("Aggregating at site level")
RLS_prod_site = calc_prod_site(RLS_prod_indiv,metadata_surveys)

# save(RLS_prod_site,file = "productivity/outputs/RLS_prod_site_with_outliers0.001.Rdata")
# 
# #Removing outliers, Biomass and productivity values superior to 99.9% of values
# RLS_prod_site = RLS_prod_site %>% filter(Biom < quantile(RLS_prod_site$Biom,0.999)) %>% 
#   filter(Productivity < quantile(RLS_prod_site$Productivity,0.999))

save(RLS_prod_site, file = "productivity/outputs/RLS_prod_site.Rdata")



#--------- species level ---------#
print("Aggregating at species level")
productivity_sp <- RLS_prod_indiv %>% 
  mutate(productivity = (Prod/Biom)*100) %>%
  group_by(Species) %>%
  mutate(mean_size = mean(Size),
         mean_Kmax = mean(Kmax),
         mean_productivity = mean(productivity)) %>%
  ungroup() %>%
  dplyr::select(Species, mean_size, mean_Kmax, mean_productivity)

save(productivity_sp, file = "productivity/outputs/RLS_prod_species.Rdata")


#--------- assess sensitivity ---------#
#Sensitivtiy to account for variability between transects
RLS_prod_sensitivity = RLS_prod_indiv %>% filter(Size < quantile(Size, 0.95)) %>% filter(Num < quantile(Num, 0.95))
RLS_prod_sensitivity_transect = calc_prod_transect(RLS_prod_sensitivity, metadata_surveys)

RLS_prod_sensitivity_transect = RLS_prod_sensitivity_transect %>%
                                  filter(Biom < quantile(Biom,0.99)) %>%
                                  filter(Productivity < quantile(Productivity,0.99))

save(RLS_prod_sensitivity_transect, file = "productivity/outputs/RLS_prod_sensitivity_transect.Rdata")


RLS_prod_sensitivity_site = calc_prod_site(RLS_prod_sensitivity, metadata_surveys)

RLS_prod_sensitivity_site = RLS_prod_sensitivity_site %>%
  filter(Biom < quantile(Biom,0.99)) %>%
  filter(Productivity < quantile(Productivity,0.99))

save(RLS_prod_sensitivity_site, file = "productivity/outputs/RLS_prod_sensitivity_site.Rdata")


summary(RLS_prod_transect)
range(RLS_prod_transect$Productivity)


