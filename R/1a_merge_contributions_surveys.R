################################################################################
##
## This script gathers all the contributions studied in this project, 
##  in one dataframe, at the survey and site scale
##
## 1a_merge_contributions_surveys.R
##
## 21/10/2022
##
## Ulysse Flandrin
##
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

# #-----------------Loading packages-------------------
# pkgs <- c("here", "tidyverse","dplyr", "questionr", "purrr", "robustbase")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

##-------------loading data-------------
#Metadata surveys
load(here::here("data","metadata_surveys.Rdata"))
#Habitat information of sites
load(here::here("data_raw", "source", "coral_reef_allen_habitat__data.RData"))
benthic_imputed <- read.csv(here::here("data_raw", "source", "RLS_benthic_imputed.txt"), sep= " ")


#Contributions
load(here::here("recycling", "outputs","flux_final_data_surveys.Rdata"))
load(here::here("productivity", "outputs","RLS_prod_transect.Rdata"))
load(here::here("nutrients", "outputs", "nutrient_concentration_surveys.Rdata"))
load(here::here("nutrients", "outputs", "fishery_tot_biomass.Rdata"))
load(here::here("biodiversity", "outputs", "surveys_biodiversity.Rdata"))
load(here::here("biodiversity", "outputs", "phylogenetic_indices_surveys.Rdata"))
load(here::here("biodiversity", "outputs", "surveys_iucn_and_chondrichtyans_scores.Rdata"))
load(here::here("biodiversity", "outputs", "survey_endemism_score.Rdata"))
load(here::here("cultural_contributions", "outputs", "survey_aesth.Rdata"))
load(here::here("cultural_contributions", "outputs", "cultural_contributions_surveys.Rdata"))
load(here::here("carbonates", "outputs", "caco3_per_day.Rdata"))
load(here::here("trophic_web", "outputs", "trophic_indicators_survey.Rdata"))


##-------------Filtering contributions-------------
recycl <- dplyr::select(task3_data_surveys,
                        SurveyID, Btot, recycling_N, recycling_P) #2402 surveys

prod <- dplyr::select(RLS_prod_transect, SurveyID, Productivity) #3628 surveys

nutrients <- dplyr::select(RLS_nut_surv,-c("biom_tot")) #3626 surveys

biodiv <- dplyr::select(surveys_biodiversity, SurveyID, taxo_richness, funct_entropy,
                        funct_distinctiveness, biom_lowTL, biom_mediumTL, biom_highTL)#3626 surveys

phylo <- dplyr::select(phylo_indices_surveys_all, SurveyID, phylo_entropy_Mean, ED_Mean)#3626 surveys

aesthetic <- dplyr::select(survey_aesth_all, SurveyID, aesthe_survey) #7013 surveys

carbonates <- dplyr::select(caco3_per_day, SurveyID, low_mg_calcite, high_mg_calcite,
                            aragonite, monohydrocalcite, amorphous_carbonate) #2169 surveys

troph_ind <- dplyr::select(as.data.frame(trophic_indicators_survey), SurveyID, 
                           b_power_law, weighted_mTL) |> 
  dplyr::rename(robustness = b_power_law, mean_TL = weighted_mTL) |>
  dplyr::mutate(across(.cols = c("SurveyID"),.fns = as.character, .names = "{.col}" )) #3628 surveys

cultural <- dplyr::select(cultural_contribution_surveys, 
                          SurveyID, public_interest) #3628 surveys

coral_cover <- dplyr::select(benthic_imputed, SurveyID, coral_imputation) |>
  dplyr::mutate( SurveyID = as.character(SurveyID))


Contrib <- metadata_surveys |>
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
  dplyr::left_join(endemism_survey_rls) |>
  dplyr::left_join(carbonates) |>
  dplyr::left_join(troph_ind) |>
  dplyr::left_join(cultural) |>
  
  dplyr::select(SurveyID, SiteCode, SiteCountry, SiteLatitude, SiteLongitude, 
                SiteMeanSST, SurveyDepth, coral_imputation,
                Biomass = Btot, 
                N_Recycling = recycling_N,
                P_Recycling = recycling_P,
                Taxonomic_Richness = taxo_richness,
                Functional_Entropy = funct_entropy, 
                Phylogenetic_Entropy = phylo_entropy_Mean,
                Traits_Distinctiveness = funct_distinctiveness,
                Evolutionary_Distinctiveness = ED_Mean,
                Herbivores_Biomass = biom_lowTL,
                Invertivores_Biomass = biom_mediumTL,
                Piscivores_Biomass = biom_highTL,
                # IUCN_Species = iucn_species,
                Endemism,
                Elasmobranch_Diversity = elasmobranch_diversity,
                Low_Mg_Calcite = low_mg_calcite,
                High_Mg_Calcite = high_mg_calcite,
                Aragonite = aragonite,
                Monohydrocalcite = monohydrocalcite,
                Amorphous_Carbonate = amorphous_carbonate,
                Trophic_Web_Robustness = robustness,
                Mean_Trophic_Level = mean_TL,
                
                Productivity = Productivity,
                Selenium = Selenium_C,
                Zinc = Zinc_C,
                Omega_3 = Omega_3_C,
                Calcium = Calcium_C,
                Iron = Iron_C,
                Vitamin_A = Vitamin_A_C,
                Fishery_Biomass = fishery_biomass,
                Aesthetic = aesthe_survey,
                Public_Interest = public_interest
                # Academic_Knowledge = academic_knowledge
  )



#### No filter
Contrib_surveys <- questionr::na.rm(Contrib) #remains 1845 surveys --> no filter on metadata

Contrib_surveys <- Contrib_surveys |>
  dplyr::left_join( dplyr::select( metadata_surveys,
                                   SurveyID, SurveyDate, SiteEcoregion, 
                                   HDI, MarineEcosystemDependency, gravtot2,
                                   mpa_name, mpa_enforcement, protection_status, mpa_iucn_cat) )
save(Contrib_surveys, file = here::here("outputs", "all_Contrib_surveys.Rdata"))


##### log right skewed distribution ###
### Contribution distribution
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
plots <- lapply(colnames(Contrib_surveys)[-c(1:8)], FUN = plot_distribution, data = Contrib_surveys )
all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + plots[[30]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))
all_plot

median <-  robustbase::colMedians(as.matrix(dplyr::select(Contrib_surveys, Biomass:Public_Interest)))
max <- apply(dplyr::select(Contrib_surveys, Biomass:Public_Interest),2,max)
magnitude <- max/median
which(magnitude > 10)

## Contributions distribution with right skewed distribution: more than one order of magnitude between median and max
Contrib_skewed_distribution <- c("Biomass","N_Recycling","P_Recycling",
                             "Low_Mg_Calcite", "High_Mg_Calcite", "Aragonite",
                             "Monohydrocalcite", "Amorphous_Carbonate", 
                             "Herbivores_Biomass", "Invertivores_Biomass", 
                             "Piscivores_Biomass", "Fishery_Biomass",
                             "gravtot2")


Contrib_survey_log_transformed <- Contrib_surveys |>
  dplyr::mutate(across(.cols = all_of(Contrib_skewed_distribution),
                       .fns = ~ .x +1 , .names = "{.col}")) |>      
  dplyr::mutate(across(.cols = all_of(Contrib_skewed_distribution),
                       .fns = log10 , .names = "{.col}"))          # log(x+1) to avoid negative values

save(Contrib_survey_log_transformed, file = here::here("outputs", "all_Contrib_survey_log_transformed.Rdata"))



### randomisation test : melt sites
set.seed(06)
Contrib_surveys_random <- questionr::na.rm(Contrib) |>
  dplyr::mutate(across(.cols = c("Biomass":"Public_Interest"), .fns = sample, .names = "{.col}"))

### keep only coral reef sites (according to Allen Atlas)
allen_list <- dplyr::mutate(hab_filt, SurveyID = as.character(SurveyID))
Contrib_coral <- dplyr::filter(Contrib, SurveyID %in% allen_list$SurveyID)
Contrib_surveys_coral_reef <- questionr::na.rm(Contrib_coral) #remains 1623

#### test without Australia
Contrib_wo_aust <- dplyr::filter(Contrib, SiteCountry != "Australia")
Contrib_surveys_wo_australia <- questionr::na.rm(Contrib_wo_aust) #remains 746


#### remove non-coral tropical reef
ID_wo_no_coral <- dplyr::filter(benthic_imputed, coral_imputation >0) |>
  dplyr::mutate(SurveyID = as.character(SurveyID))
Contrib_wo_no_coral <- dplyr::filter(Contrib, SurveyID %in% ID_wo_no_coral$SurveyID)
Contrib_surveys_coral_0_imputed <- questionr::na.rm(Contrib_wo_no_coral) #remains 1720


#### tropical reef with more than 5% coral cover
ID_wo_5_coral <- dplyr::filter(benthic_imputed, coral_imputation >0.05) |>
  dplyr::mutate(SurveyID = as.character(SurveyID))
Contrib_wo_no_coral <- dplyr::filter(Contrib, SurveyID %in% ID_wo_5_coral$SurveyID)
Contrib_surveys_coral_5_imputed <- questionr::na.rm(Contrib_wo_no_coral) #remains 1510


#### tropical reef with more than 10% coral cover
ID_wo_10_coral <- dplyr::filter(benthic_imputed, coral_imputation >0.1) |>
  dplyr::mutate(SurveyID = as.character(SurveyID))
Contrib_wo_no_coral <- dplyr::filter(Contrib, SurveyID %in% ID_wo_10_coral$SurveyID)
Contrib_surveys_coral_10_imputed <- questionr::na.rm(Contrib_wo_no_coral) #remains 1346



#### Site with minSST > 20
survey_sst20 <- dplyr::filter(metadata_surveys, SiteMinSST > 20)
Contrib_surveys_SST20 <- dplyr::filter(Contrib, SurveyID %in% survey_sst20$SurveyID) |>
  questionr::na.rm() #remains 1514


#### test with only Australia
Contrib_only_aust <- dplyr::filter(Contrib, SiteCountry == "Australia")
Contrib_surveys_only_australia <- questionr::na.rm(Contrib_only_aust) #remains 1084



### Mean per site
mean_per_site <- function(var_survey = Contrib_surveys){
  var_survey |> dplyr::group_by( SiteCode, SiteCountry) |>
    dplyr::summarise(across(.cols = c("SiteLatitude":"Public_Interest"),
                            .fns = mean, .names = "{.col}"))
}


Contrib_site <- mean_per_site(Contrib_surveys)
summary(Contrib_site) # 1237 sites
Contrib_site <- Contrib_site |>
  dplyr::left_join( 
    dplyr::distinct(metadata_surveys,
                    SiteCode, SurveyDate, SiteEcoregion,
                    HDI, MarineEcosystemDependency, gravtot2,
                    mpa_name, mpa_enforcement, protection_status, mpa_iucn_cat))
save(Contrib_site, file = here::here("outputs", "all_Contrib_site.Rdata"))


Contrib_site_log_transformed <- Contrib_site |>
  dplyr::mutate(across(.cols = all_of(Contrib_skewed_distribution),
                       .fns = ~ .x +1 , .names = "{.col}")) |>      
  dplyr::mutate(across(.cols = all_of(Contrib_skewed_distribution),
                       .fns = log10 , .names = "{.col}"))          # log(x+1) to avoid negative values

save(Contrib_site_log_transformed, file = here::here("outputs", "all_Contrib_site_log_transformed.Rdata"))



#save with different filter
for(condition in c("coral_reef", "coral_0_imputed", "wo_australia", "coral_5_imputed",  "coral_10_imputed",
                   "only_australia", "SST20", "random")){
  Contrib_surveys_condition <- get(paste0("Contrib_surveys_", condition))
  Contrib_site_condition <- mean_per_site(Contrib_surveys_condition)
  Contrib_site_condition <- Contrib_site_condition |>
    dplyr::left_join( dplyr::distinct( metadata_surveys,
                                       SiteCode, SurveyDate, SiteEcoregion,
                                       HDI, MarineEcosystemDependency, gravtot2,
                                       mpa_name, mpa_enforcement, protection_status, mpa_iucn_cat),
                      multiple = "all") |>
    dplyr::mutate(across(.cols = all_of(Contrib_skewed_distribution),
                         .fns = ~ .x +1 , .names = "{.col}")) |>      # Adds 1 to values to log transformed
    dplyr::mutate(across(.cols = all_of(Contrib_skewed_distribution),
                         .fns = log10 , .names = "{.col}"))          # log(x+1) to avoid negative values
  
  save(Contrib_site_condition, file = here::here("outputs", paste0("Contrib_site_log_", condition,".Rdata")))
}


