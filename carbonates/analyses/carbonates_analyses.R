################################################################################
##'
##' clean carbonates outputs (file: `carbonates/outputs/survey_caco3_composition.rds`) 
##'  from Ghiraldi et al. 2023
##'
##' carbonates_analyses.R
##'
##' 31/10/2022
##'
##' Ulysse Flandrin
##'
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

# #-----------------Loading packages-------------------
# pkgs <- c("here", "tidyverse", "dplyr")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

##-------------loading data-------------
caco3_composition <- readRDS(file = here::here("carbonates", "outputs", "survey_caco3_composition.rds")) 
caco3_production <- readRDS(file = here::here("carbonates", "outputs", "survey_caco3_production.rds")) 

##-------------merging data-------------
caco3_per_day <- caco3_production  |> 
  dplyr::left_join(caco3_composition)  |> 
  dplyr::select(SiteCode, SurveyID, prop_biomass, 
                carbonate_tot = caco3_umol_day, #total production of carbonate in a survey in umol/day
                low_mg_calcite = L_umol_day,
                high_mg_calcite = H_umol_day,
                aragonite = AR_umol_day,
                monohydrocalcite = M_umol_day,
                amorphous_carbonate = AC_umol_day
  )

caco3_per_day$SurveyID <- as.character(caco3_per_day$SurveyID)


##study correlations
rownames(caco3_per_day) <- caco3_per_day$SurveyID
pca <- FactoMineR::PCA(caco3_per_day[,-c(1:3)], scale.unit = T, graph=T, ncp=10)
factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 60))
var <- factoextra::get_pca_var(pca)
corrplot::corrplot(var$contrib, is.corr=FALSE)  
factoextra::fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE )
factoextra::fviz_pca_var(pca, col.var = "cos2",
             axe=c(3,4),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


##-------------saving data-------------
save(caco3_per_day, file=  here::here("carbonates", "outputs", "caco3_per_day.Rdata") )
