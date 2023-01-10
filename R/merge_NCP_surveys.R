################################################################################
##
## Gather all NCPs in one dataframe, at the survey scale
##
## merge_NCP_surveys.R
##
## 21/10/2022
##
## Ulysse Flandrin
##
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse","dplyr", "questionr")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

##-------------loading data-------------
load(here::here("data","metadata_surveys.Rdata"))
load(here::here("recycling", "outputs","flux_final_data_surveys.Rdata"))
load(here::here("productivity", "outputs","RLS_prod_transect.Rdata"))
load(here::here("nutrients", "outputs", "nutrient_concentration_surveys.Rdata"))
load(here::here("nutrients", "outputs", "fishery_tot_biomass.Rdata"))
load(here::here("biodiversity", "outputs", "surveys_biodiversity.Rdata"))
load(here::here("biodiversity", "outputs", "phylogenetic_indices_surveys_without_outliers.Rdata"))
load(here::here("biodiversity", "outputs", "surveys_iucn_and_chondrichtyans_scores.Rdata"))
load(here::here("aesthetic", "outputs", "survey_aesth_without_outliers.Rdata"))
load( here::here("carbonates", "outputs", "caco3_per_day_without_outliers.Rdata"))

#list of coral reef sites
load(here::here("data_raw", "source", "coral_reef_allen_habitat__data.RData"))
benthic_imputed <- read.csv(here::here("data_raw", "source", "RLS_benthic_imputed.txt"), sep= " ")

#world coast shapefile
#coast_high_resolution <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))
coast <- sf::st_read(here::here("data", "coastline_shapefile", "ne_10m_coastline.shp"))

##-------------Filtering NCP-------------
recycl <- dplyr::select(task3_data_surveys,
                        SurveyID, Btot, recycling_N, recycling_P) #2399 surveys

prod <- dplyr::select(RLS_prod_transect, SurveyID, Productivity) #3619 surveys

nutrients <- RLS_nut_without_outliers[,-2] #3606 surveys

biodiv <- dplyr::select(surveys_biodiversity_without_outliers, SurveyID, taxo_richness, funct_entropy,
                        funct_distinctiveness, biom_lowTL, biom_mediumTL, biom_highTL)#3610 surveys

phylo <- dplyr::select(phylo_indices_surveys, SurveyID, phylo_entropy, ED_Mean)#3613 surveys

aesthetic <- dplyr::select(survey_aesth, SurveyID, aesthe_survey) #7005 surveys

carbonates <- dplyr::select(caco3_per_day_without_outliers, SurveyID, low_mg_calcite, high_mg_calcite,
                            aragonite, monohydrocalcite, amorphous_carbonate) #2156 surveys

coral_cover <- dplyr::select(benthic_imputed, SurveyID, coral_imputation) %>%
  mutate( SurveyID = as.character(SurveyID))

NCP <- metadata_surveys %>%
  select(SurveyID, SiteCode, SiteCountry, SiteEcoregion, SiteLatitude, SiteLongitude, 
         SurveyDepth, SiteMeanSST, gravtot2, HDI, MarineEcosystemDependency ) %>%
  left_join(coral_cover) %>%
  left_join(recycl) %>%
  left_join(prod) %>%
  left_join(biodiv) %>%
  left_join(nutrients) %>%
  left_join(surveys_fishery_biom) %>%
  left_join(phylo) %>%
  left_join(aesthetic) %>%
  left_join(iucn_and_elasmo_by_surveys)%>%
  left_join(carbonates)


#### No filter
NCP_surveys <- questionr::na.rm(NCP) #remains 1817 surveys with HDI and Marine dependency / 1792 with gravity
save(NCP_surveys, file = here::here("outputs", "all_NCP_surveys.Rdata"))


#### keep only coral reef sites (according to Allen Atlas)
# allen_list <- dplyr::mutate(hab_filt, SurveyID = as.character(SurveyID))
# NCP_coral <- dplyr::filter(NCP, SurveyID %in% allen_list$SurveyID)
# NCP_surveys <- questionr::na.rm(NCP_coral) #remains 1576

# #### test without Australia
# NCP_wo_aust <- dplyr::filter(NCP, SiteCountry != "Australia")
# NCP_surveys <- questionr::na.rm(NCP_wo_aust) #remains 731

# #### test without Spain
# NCP_wo_spain <- dplyr::filter(NCP, SiteCountry != "Spain")
# NCP_surveys <- questionr::na.rm(NCP_wo_spain) #remains 1576

# #### remove non-coral tropical reef
# ID_wo_no_coral <- dplyr::filter(benthic_imputed, coral_imputation >0) %>%
#   dplyr::mutate(SurveyID = as.character(SurveyID))
# NCP_wo_no_coral <- dplyr::filter(NCP, SurveyID %in% ID_wo_no_coral$SurveyID)
# NCP_surveys <- questionr::na.rm(NCP_wo_no_coral) #remains 1682


# #### tropical reef with more than 1% coral cover
# ID_wo_1_coral <- dplyr::filter(benthic_imputed, coral_imputation >0.01) %>%
#   dplyr::mutate(SurveyID = as.character(SurveyID))
# NCP_wo_no_coral <- dplyr::filter(NCP, SurveyID %in% ID_wo_1_coral$SurveyID)
# NCP_surveys <- questionr::na.rm(NCP_wo_no_coral) #remains 1607


# #### tropical reef with more than 5% coral cover
# ID_wo_5_coral <- dplyr::filter(benthic_imputed, coral_imputation >0.05) %>%
#   dplyr::mutate(SurveyID = as.character(SurveyID))
# NCP_wo_no_coral <- dplyr::filter(NCP, SurveyID %in% ID_wo_5_coral$SurveyID)
# NCP_surveys <- questionr::na.rm(NCP_wo_no_coral) #remains 1607


# #### tropical reef with more than 10% coral cover
# ID_wo_10_coral <- dplyr::filter(benthic_imputed, coral_imputation >0.1) %>%
#   dplyr::mutate(SurveyID = as.character(SurveyID))
# NCP_wo_no_coral <- dplyr::filter(NCP, SurveyID %in% ID_wo_10_coral$SurveyID)
# NCP_surveys <- questionr::na.rm(NCP_wo_no_coral) #remains 1313



# #### Site with minSST > 20
# survey_sst20 <- dplyr::filter(metadata_surveys, SiteMinSST > 20)
# NCP_surveys <- dplyr::filter(NCP, SurveyID %in% survey_sst20$SurveyID) %>%
#   questionr::na.rm() #remains 1475


# #### test with only Australia
# NCP_only_aust <- dplyr::filter(NCP, SiteCountry == "Australia")
# NCP_surveys <- questionr::na.rm(NCP_only_aust) #remains 1061

# #### test with only in french polynesia and galapagos
# NCP_only_aust <- dplyr::filter(NCP, SiteCountry == "French Polynesia" | SiteCountry == "Ecuador")
# NCP_surveys <- questionr::na.rm(NCP_only_aust) #remains 1061

# #### test without Australia, without spain
# NCP_wo_aust_spain <- dplyr::filter(NCP, SiteCountry != "Australia" & SiteCountry != "Spain")
# NCP_surveys <- questionr::na.rm(NCP_wo_aust_spain) #remains 672


#Mean per site
NCP_site <- NCP_surveys %>%
  group_by( SiteCode, SiteCountry, SiteEcoregion) %>%
  summarise(across(.cols = c("SurveyDepth","SiteMeanSST", "gravtot2","HDI", "MarineEcosystemDependency",
                        "SiteLatitude", "SiteLongitude", "coral_imputation",
                       "Btot","recycling_N","recycling_P","Productivity","taxo_richness", "funct_entropy",
                       "funct_distinctiveness","Selenium_C","Zinc_C","Omega_3_C","Calcium_C","Iron_C","Vitamin_A_C",
                       "phylo_entropy","ED_Mean","aesthe_survey", "iucn_species", "elasmobranch_diversity",
                       "low_mg_calcite", "high_mg_calcite", "aragonite", "monohydrocalcite",
                       "amorphous_carbonate", "biom_lowTL", "biom_mediumTL", "biom_highTL", "fishery_biomass"),
                   .fns = mean, .names = "{.col}"))

summary(NCP_site) # 1223 sites

save(NCP_site, file = here::here("outputs", "all_NCP_site.Rdata"))

##-------------plot NCPs on map-------------
#plot function
plot_NCP_on_world_map <- function(NCP = "taxo_richness"){
  map <- ggplot(NCP_site) +
    geom_sf(data = coast, color = "grey30", fill = "lightgrey",
            aes(size=0.1)) +
    
    geom_point(data=NCP_site,
               size = 2, alpha = 0.7, shape = 20,
               aes(x = SiteLongitude, y = SiteLatitude,
                   colour= NCP_site[,NCP][[1]])) +
    scale_colour_gradient(NCP,
                          low = "dodgerblue", high="darkred",
                          na.value=NA) +
    
    coord_sf(ylim = c(-36, 31), expand = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_minimal()+
    labs(title = paste0(NCP, " geographic distribution"),
         x="", y= "") +
    theme(legend.position = "right",
          plot.title = element_text(size=10, face="bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
    
  ggsave( here::here("outputs", "figures", "NCP_on_world_map", paste("world map with ", NCP, ".jpg")), plot = map, width=15, height = 7 )
  #map
}

# save maps
NCP_list <- colnames(subset(NCP_site, select = -c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                                  SiteMeanSST, SiteLatitude, SiteLongitude,
                                                  HDI, gravtot2, MarineEcosystemDependency,
                                                  coral_imputation)))
for( NCP in NCP_list){
  plot_NCP_on_world_map(NCP)
}
