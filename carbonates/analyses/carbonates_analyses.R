################################################################################
##
## clean carbonates outputs
##
## carbonates_analyses.R
##
## 31/10/2022
##
## Ulysse Flandrin
##
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse", "dplyr")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

##-------------loading data-------------
caco3_composition <- readRDS(file = here::here("carbonates", "outputs", "survey_caco3_composition.rds")) 
caco3_production <- readRDS(file = here::here("carbonates", "outputs", "survey_caco3_production.rds")) 

##-------------merging data-------------
caco3_per_day <- caco3_production %>%
  left_join(caco3_composition) %>%
  select(SiteCode, SurveyID, prop_biomass, 
         carbonate_tot = caco3_umol_day, #total production of carbonate in a survey in umol/day
         low_mg_calcite = L_umol_day,
         high_mg_calcite = H_umol_day,
         aragonite = AR_umol_day,
         monohydrocalcite = M_umol_day,
         amorphous_carbonate = AC_umol_day
         )

caco3_per_day$SurveyID <- as.character(caco3_per_day$SurveyID)

# ## filtering 0.1% outliers
# caco3_per_day_without_outliers <- caco3_per_day %>%
#   filter(carbonate_tot < quantile(caco3_per_day$carbonate_tot,0.999)) %>% 
#   filter(low_mg_calcite < quantile(caco3_per_day$low_mg_calcite,0.999)) %>% 
#   filter(high_mg_calcite < quantile(caco3_per_day$high_mg_calcite,0.999)) %>% 
#   filter(aragonite < quantile(caco3_per_day$aragonite,0.999)) %>% 
#   filter(monohydrocalcite < quantile(caco3_per_day$monohydrocalcite,0.999)) %>% 
#   filter(amorphous_carbonate < quantile(caco3_per_day$amorphous_carbonate,0.999)) 

##study correlations
rownames(caco3_per_day) <- caco3_per_day$SurveyID
pca <- FactoMineR::PCA(caco3_per_day[,-c(1:3)], scale.unit = T, graph=T, ncp=10)
factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 60))
var <- get_pca_var(pca)
corrplot(var$contrib, is.corr=FALSE)  
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE )
fviz_pca_var(pca, col.var = "cos2",
             axe=c(3,4),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
site_coord_in_pca <- as.data.frame(pca$ind$coord) 
plot_ly(site_coord_in_pca, 
        x= ~Dim.1, y= ~Dim.2, z= ~Dim.3,
        size = 10)

##-------------saving data-------------
save(caco3_per_day, file=  here::here("carbonates", "outputs", "caco3_per_day.Rdata") )
# save(caco3_per_day_without_outliers, file=  here::here("carbonates", "outputs", "caco3_per_day_without_outliers.Rdata") )
