################################################################################
##
## Extract trophic levels from fishbase and assess common body size of species 
##
## common_length_and_TL_of_species.R
##
## 09/02/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse", "rfishbase", "dplyr", "DHARMa")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load(here::here("data","data_species.Rdata"))


##-------------extract data from fishbase-------------
#raw data
species_traits <- rfishbase::species( server = "fishbase") 
species_ecology <- rfishbase::ecology( server = "fishbase")

#filter data
info_species <- species_traits |>
  dplyr::select(SpecCode, Species, DemersPelag) |>
  dplyr::left_join( dplyr::select(species_ecology,SpecCode, FoodTroph)) |>
  dplyr::filter(SpecCode %in% data_species$spec_code)

unfound_sp <- data_species$species[which(
  data_species$spec_code %in% setdiff(data_species$spec_code, info_species$SpecCode))]
#  Kyphosus_analogus should be Kyphosus vaigiensis

#merge all data traits
info_species <- info_species |>
  dplyr::rename( spec_code = SpecCode, trophic_level = FoodTroph,
               water_column = DemersPelag, species_fishbase_name = Species) |>
  dplyr::right_join(data_species)


##-------------infer common length from max length-------------
### Common length ~ max length relationship
#check normality
hist(log10(species_traits$Length))
hist(log10(species_traits$CommonLength))

reg <- lm(log10(species_traits$CommonLength)~log10(species_traits$Length))
#Check the regression
plot(DHARMa::simulateResiduals(reg)) 
summary(reg)

#plot the regression
library(ggplot2)
ggplot(species_traits, aes(y=log10(CommonLength), x=log10(Length)))+
  geom_point()+
  geom_smooth(colour="red", method="lm", fill="red") +
  ylab("log(Common lenght)")+
  xlab("log(Max length)") +
  theme_classic()+
  annotate("text", x = 0.8, y = 2.5, label = "log(Common_length) = 0.901*log(Max_length) \n Common_length = Max_length^0.901 \n (r² = 0.899 , pval < 10e-16)") 
ggsave(here::here("trophic_web", "outputs", "log_linear_relationship_commonlength_maxlength.png"), width = 12, height = 8)


### Extract common length
data_species_trophic_web <- info_species |>
  dplyr::mutate( common_length = Size^reg[["coefficients"]][[2]]) #common_length = max_length^0.901
 

##-------------Complete data from the median of the gender-------------
# names_na_length <- rownames(data_species_trophic_web[which(
#   is.na(data_species_trophic_web$common_length)),])
# species_length <- list()
# 
# for(i in 1:length(names_na_length)){
#   cat("i : ",i,"\n")
#   gender_length <- data_species_trophic_web[grep(do.call(rbind,strsplit(names_na_length,"_"))[i,1],
#                                      rownames(data_species_trophic_web)),"common_length"]
#   species_length[[i]] <- median(gender_length,na.rm=TRUE)
# } 
# 
# names(species_length) <- names_na_length
# species_length <- do.call(rbind,species_length)
# 
# data_species_trophic_web[rownames(species_length),"common_length"]  <- species_length[,1]
# data_species_trophic_web <- data_species_trophic_web[-which(is.na(data_species_trophic_web[,"common_length"])),]



##-------------save final data-------------
save(data_species_trophic_web, file = here::here("trophic_web", "outputs", "data_species_trophic_web.Rdata"))

