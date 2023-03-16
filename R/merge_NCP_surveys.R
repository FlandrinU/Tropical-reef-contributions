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
load(here::here("cultural_contributions", "outputs", "survey_aesth_without_outliers.Rdata"))
load(here::here("cultural_contributions", "outputs", "cultural_contributions_surveys.Rdata"))
load(here::here("carbonates", "outputs", "caco3_per_day_without_outliers.Rdata"))
load(here::here("trophic_web", "outputs", "trophic_indicators_survey.Rdata"))

#list of coral reef sites
load(here::here("data_raw", "source", "coral_reef_allen_habitat__data.RData"))
benthic_imputed <- read.csv(here::here("data_raw", "source", "RLS_benthic_imputed.txt"), sep= " ")

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

troph_ind <- dplyr::select(as.data.frame(trophic_indicators_survey), SurveyID, 
                           b_power_law, weighted_mTL) |> #3627 surveys
             dplyr::rename(robustness = b_power_law, mean_TL = weighted_mTL) |>
             dplyr::mutate(across(.cols = c("SurveyID"),.fns = as.character, .names = "{.col}" ))
 

coral_cover <- dplyr::select(benthic_imputed, SurveyID, coral_imputation) |>
  dplyr::mutate( SurveyID = as.character(SurveyID))

NCP <- metadata_surveys |>
  dplyr::select(SurveyID, SiteCode, SiteCountry, SiteLatitude, SiteLongitude, 
                SiteMeanSST, SurveyDepth) |> 
  dplyr::left_join(coral_cover) |>
  dplyr::left_join(recycl) |>
  dplyr::left_join(prod) |>
  dplyr::left_join(biodiv) |>
  dplyr::left_join(nutrients) |>
  dplyr::left_join(surveys_fishery_biom) |>
  dplyr::left_join(phylo) |>
  dplyr::left_join(aesthetic) |>
  dplyr::left_join(iucn_and_elasmo_by_surveys)|>
  dplyr::left_join(carbonates) |>
  dplyr::left_join(troph_ind) |>
  dplyr::left_join(
    dplyr::select(cultural_contribution_surveys, 
                  SurveyID, academic_knowledge, public_interest)) |>
  
  dplyr::select(SurveyID, SiteCode, SiteCountry, SiteLatitude, SiteLongitude, 
                SiteMeanSST, SurveyDepth, coral_imputation,
                Biomass = Btot, 
                N_recycling = recycling_N, P_recycling = recycling_P,
                Taxonomic_Richness = taxo_richness,
                Functional_Entropy = funct_entropy, 
                Phylogenetic_Entropy = phylo_entropy,
                Functional_Distinctiveness = funct_distinctiveness,
                Evolutionary_distinctiveness = ED_Mean,
                Low_TL_Biomass = biom_lowTL,
                Medium_TL_Biomass = biom_mediumTL,
                High_TL_Biomass = biom_highTL,
                IUCN_Species = iucn_species,
                Elasmobranch_Diversity = elasmobranch_diversity,
                Low_Mg_Calcite = low_mg_calcite,
                High_Mg_Calcite = high_mg_calcite,
                Aragonite = aragonite,
                Monohydrocalcite = monohydrocalcite,
                Amorphous_Carbonate = amorphous_carbonate,
                Trophic_web_robustness = robustness,
                mean_Trophic_Level = mean_TL,
                
                Productivity = Productivity,
                Selenium = Selenium_C,
                Zinc = Zinc_C,
                Omega_3 = Omega_3_C,
                Calcium = Calcium_C,
                Iron = Iron_C,
                Vitamin_A = Vitamin_A_C,
                Fishery_Biomass = fishery_biomass,
                Aesthetic = aesthe_survey,
                Public_Interest = public_interest,
                Academic_Knowledge = academic_knowledge
                )
  


#### No filter
NCP_surveys <- questionr::na.rm(NCP) #remains 1809 surveys --> no filter on metadata

NCP_surveys <- NCP_surveys |>
  dplyr::left_join( dplyr::select( metadata_surveys,
                             SurveyID, SiteEcoregion, 
                             HDI, MarineEcosystemDependency, gravtot2,
                             mpa_name, mpa_enforcement, protection_status, mpa_iucn_cat) )
  
save(NCP_surveys, file = here::here("outputs", "all_NCP_surveys.Rdata"))


### keep only coral reef sites (according to Allen Atlas)
allen_list <- dplyr::mutate(hab_filt, SurveyID = as.character(SurveyID))
NCP_coral <- dplyr::filter(NCP, SurveyID %in% allen_list$SurveyID)
NCP_surveys_coral_reef <- questionr::na.rm(NCP_coral) #remains 1606

#### test without Australia
NCP_wo_aust <- dplyr::filter(NCP, SiteCountry != "Australia")
NCP_surveys_wo_australia <- questionr::na.rm(NCP_wo_aust) #remains 732

# #### test without Spain
# NCP_wo_spain <- dplyr::filter(NCP, SiteCountry != "Spain")
# NCP_surveys_wo_spain <- questionr::na.rm(NCP_wo_spain) #remains 1576

#### remove non-coral tropical reef
ID_wo_no_coral <- dplyr::filter(benthic_imputed, coral_imputation >0) |>
  dplyr::mutate(SurveyID = as.character(SurveyID))
NCP_wo_no_coral <- dplyr::filter(NCP, SurveyID %in% ID_wo_no_coral$SurveyID)
NCP_surveys_coral_0_imputed <- questionr::na.rm(NCP_wo_no_coral) #remains 1714


# #### tropical reef with more than 1% coral cover
# ID_wo_1_coral <- dplyr::filter(benthic_imputed, coral_imputation >0.01) |>
#   dplyr::mutate(SurveyID = as.character(SurveyID))
# NCP_wo_no_coral <- dplyr::filter(NCP, SurveyID %in% ID_wo_1_coral$SurveyID)
# NCP_surveys <- questionr::na.rm(NCP_wo_no_coral) #remains 1607


#### tropical reef with more than 5% coral cover
ID_wo_5_coral <- dplyr::filter(benthic_imputed, coral_imputation >0.05) |>
 dplyr::mutate(SurveyID = as.character(SurveyID))
NCP_wo_no_coral <- dplyr::filter(NCP, SurveyID %in% ID_wo_5_coral$SurveyID)
NCP_surveys_coral_5_imputed <- questionr::na.rm(NCP_wo_no_coral) #remains 1505


#### tropical reef with more than 10% coral cover
ID_wo_10_coral <- dplyr::filter(benthic_imputed, coral_imputation >0.1) |>
 dplyr::mutate(SurveyID = as.character(SurveyID))
NCP_wo_no_coral <- dplyr::filter(NCP, SurveyID %in% ID_wo_10_coral$SurveyID)
NCP_surveys_coral_10_imputed <- questionr::na.rm(NCP_wo_no_coral) #remains 1335



#### Site with minSST > 20
survey_sst20 <- dplyr::filter(metadata_surveys, SiteMinSST > 20)
NCP_surveys_SST20 <- dplyr::filter(NCP, SurveyID %in% survey_sst20$SurveyID) |>
  questionr::na.rm() #remains 1506


#### test with only Australia
NCP_only_aust <- dplyr::filter(NCP, SiteCountry == "Australia")
NCP_surveys_only_australia <- questionr::na.rm(NCP_only_aust) #remains 1092

# #### test with only in french polynesia and galapagos
# NCP_only_aust <- dplyr::filter(NCP, SiteCountry == "French Polynesia" | SiteCountry == "Ecuador")
# NCP_surveys <- questionr::na.rm(NCP_only_aust) #remains 1061

# #### test without Australia, without spain
# NCP_wo_aust_spain <- dplyr::filter(NCP, SiteCountry != "Australia" & SiteCountry != "Spain")
# NCP_surveys <- questionr::na.rm(NCP_wo_aust_spain) #remains 672


### NCP distribution
plot_distribution <- function(ncp, data){
  col <- fishualize::fish(n = ncol(data), option = "Ostracion_whitleyi", begin = 0, end = 0.8)
  names(col) <- colnames(data)
  
  ggplot(data) +
    aes(x = data[,ncp]) +
    geom_histogram(bins = 40L,
                   fill = col[ncp][[1]],
                   col = "black") +
    labs(title = ncp) +
    xlab("") + ylab("") +
    theme_minimal() +
    theme( legend.position = "none")
}
library(ggplot2)
library(patchwork)
plots <- lapply(colnames(NCP_surveys)[-c(1:8)], FUN = plot_distribution, data = NCP_surveys )
all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + plots[[30]] + plots[[31]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))
all_plot

## NCPs distribution with right skewed distribution
NCP_skewed_distribution <- c("Biomass","N_recycling","P_recycling","Productivity",
                             "Functional_Distinctiveness","Omega_3","Calcium","Vitamin_A",
                             "Phylogenetic_Entropy","Evolutionary_distinctiveness",
                             "IUCN_Species", "Elasmobranch_Diversity",
                             "Low_Mg_Calcite", "High_Mg_Calcite", "Aragonite", "Monohydrocalcite",
                             "Amorphous_Carbonate", "Low_TL_Biomass", "Medium_TL_Biomass", "High_TL_Biomass",
                             "Fishery_Biomass", "mean_Trophic_Level")


### Mean per site
mean_per_site_and_log <- function(var_survey = NCP_surveys){
  var_survey |> dplyr::group_by( SiteCode, SiteCountry) |>
    dplyr::summarise(across(.cols = c("SurveyDepth","SiteMeanSST", "SiteLatitude", "SiteLongitude", "coral_imputation",
                               "Biomass","N_recycling","P_recycling","Taxonomic_Richness", "Functional_Entropy",
                               "Phylogenetic_Entropy", "Functional_Distinctiveness","Evolutionary_distinctiveness",
                               "Low_TL_Biomass", "Medium_TL_Biomass", "High_TL_Biomass","IUCN_Species",
                               "Elasmobranch_Diversity", "Low_Mg_Calcite", "High_Mg_Calcite", "Aragonite", 
                               "Monohydrocalcite", "Amorphous_Carbonate", "Trophic_web_robustness",
                               "mean_Trophic_Level", 
                               "Productivity","Selenium","Zinc","Omega_3","Calcium","Iron","Vitamin_A",
                               "Fishery_Biomass", "Aesthetic","Public_Interest","Academic_Knowledge" ), #/!\ the order matter
                     .fns = mean, .names = "{.col}"))
}


NCP_site <- mean_per_site_and_log(NCP_surveys)
summary(NCP_site) # 1229 sites
NCP_site <- NCP_site |>
  dplyr::left_join( 
      dplyr::distinct(metadata_surveys,
                          SiteCode, SiteEcoregion,
                          HDI, MarineEcosystemDependency, gravtot2,
                          mpa_name, mpa_enforcement, protection_status, mpa_iucn_cat))
save(NCP_site, file = here::here("outputs", "all_NCP_site.Rdata"))


NCP_site_log_transformed <- NCP_site |>
  dplyr::mutate(across(.cols = all_of(NCP_skewed_distribution),
                       .fns = ~ .x +1 , .names = "{.col}")) |>      # Adds 1 to values to log transformed
  dplyr::mutate(across(.cols = all_of(NCP_skewed_distribution),
                       .fns = log10 , .names = "{.col}"))          # log(x+1) to avoid negative values

save(NCP_site_log_transformed, file = here::here("outputs", "all_NCP_site_log_transformed.Rdata"))



#save with different filter
for(condition in c("coral_reef", "coral_0_imputed", "wo_australia", "coral_5_imputed",  "coral_10_imputed",
"only_australia", "SST20")){
  NCP_surveys_condition <- get(paste0("NCP_surveys_", condition))
  NCP_site_condition <- mean_per_site_and_log(NCP_surveys_condition)
  NCP_site_condition <- NCP_site_condition |>
    dplyr::left_join( dplyr::distinct( metadata_surveys,
                                     SiteCode, SiteEcoregion,
                                     HDI, MarineEcosystemDependency, gravtot2,
                                     mpa_name, mpa_enforcement, protection_status, mpa_iucn_cat),
                      multiple = "all") |>
    dplyr::mutate(across(.cols = all_of(NCP_skewed_distribution),
                         .fns = ~ .x +1 , .names = "{.col}")) |>      # Adds 1 to values to log transformed
    dplyr::mutate(across(.cols = all_of(NCP_skewed_distribution),
                         .fns = log10 , .names = "{.col}"))          # log(x+1) to avoid negative values
  
  save(NCP_site_condition, file = here::here("outputs", paste0("NCP_site_log_", condition,".Rdata")))
}


