################################################################################
##
## Assess species endemism according to their geographic range (Duhamet et al. 2023)
##  and define an endemism score in each surveys
##
## endemism_survey_index.R
##
## 21/03/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c()
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load(here::here("data","metadata_surveys.Rdata"))
load(here::here("data","data_species.Rdata"))
load( here::here("biodiversity", "outputs", "occurrence_matrix_sp_survey_01.Rdata"))

mat_PA_teleo <- readRDS(here::here("biodiversity", "data","Mat_Pa_teleo.RDS"))
mat_PA_chond <- readRDS(here::here("biodiversity", "data","Mat_Pa_chond.RDS"))

coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

##-------------compute species endemism-------------
# Deal with species names
names <- gsub(" ", "_", colnames(mat_PA_teleo[, -c(1,2)]))

for( i in 1:length(names)){
  if(names[i] %in% data_species$species){
    names[i] <- data_species$species_corrected[which(data_species$species == names[i])]
  }
}
# data_species$species_corrected[which(!(data_species$species_corrected %in% names))] #6 species without range (old names or too few points on OBIS)

colnames(mat_PA_teleo) <- c("Longitude", "Latitude", names)
mat_PA_rls <- mat_PA_teleo[, which(colnames(mat_PA_teleo) %in% data_species$species_corrected)]

## species range
range_sp <- as.data.frame(colSums(mat_PA_rls)) 
colnames(range_sp) <- "range"

# #check ranges
#     summary(range_sp$range)
#     head(range_sp)
#     test <- c( rownames(range_sp)[which(range_sp$range == max(range_sp$range))],
#                rownames(range_sp)[which(range_sp$range == min(range_sp$range))],
#                rownames(range_sp)[which(range_sp$range == median(range_sp$range))] )
#     mat_test <- mat_PA_rls |> dplyr::select(dplyr::all_of(test)) |>
#       dplyr::bind_cols(mat_PA_teleo[, c(1,2)])
#     
#     library(ggplot2)
#     ggplot(mat_test) +
#       aes(x = Longitude, y = Latitude,
#           fill =  Diproctacanthus_xanthurus ) +
#       geom_tile() +
#       scale_fill_gradient() +
#       theme_minimal()

## Endemism
endemism_sp <- range_sp |>
  dplyr::mutate(endemism = (max(range) - range)/(max(range) - min(range)))
hist(endemism_sp$endemism, breaks = 20)

##-------------survey endemic score------------- = mean of species endemism in a given survey
endemism_survey <- rep(NA, nrow(surveys_sp_occ))
i=3465
for(i in 1:nrow(surveys_sp_occ)){
  Names <- names(surveys_sp_occ[i,which(surveys_sp_occ[i,] == 1)])
  corrected_names <- dplyr::filter(data_species, species %in% Names)$species_corrected
  endemism_survey[i] <- mean(endemism_sp[corrected_names, "endemism"], na.rm=T)
}

endemism_survey_rls <- data.frame(SurveyID = rownames(surveys_sp_occ), Endemism = endemism_survey)
save(endemism_survey_rls, file = here::here("biodiversity", "outputs", "survey_endemism_score.Rdata"))

#check survey endemism
hist(endemism_survey_rls$Endemism, breaks = 20)
endemism_survey_map <- endemism_survey_rls |>
  dplyr::left_join( dplyr::select(metadata_surveys, SurveyID, SiteLongitude, SiteLatitude))

ggplot(endemism_survey_map) +
  geom_sf(data = coast, color = "grey30", fill = "lightgrey",
          aes(size=0.1)) +
  geom_point(data=endemism_survey_map,
             size = 4, shape = 20,
             aes(x = SiteLongitude, y = SiteLatitude,
                 colour= Endemism)) +
  scale_colour_gradient(low = "dodgerblue", high="darkred",
                        na.value="black") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_minimal()+
  labs(title = paste0("Endemism", " geographic distribution"),
       x="", y= "") +
  theme(legend.position = "right",
        plot.title = element_text(size=10, face="bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
  )
