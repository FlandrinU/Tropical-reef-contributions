################################################################################
##
## Script in order to test the weighted correlation mean from Kark 2002
##
## test_weighted_mean_scores.R
##
## 06/02/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse","FactoMineR", "tibble", "questionr", "corrplot",
          "factoextra", "ggpubr", "scico", "RColorBrewer", "plotly", "fishualize", 
          "ggplot2", "patchwork", "colormap", "grDevices", "ggnewscale", "sf",
          "rdacca.hp", "faux")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load(here::here("outputs","all_NCP_site.Rdata"))
# load(here::here("outputs","NCP_site_coral_reef.Rdata"))
# load(here::here("outputs","NCP_site_SST20.Rdata"))
# load(here::here("outputs","NCP_site_coral_5_imputed.Rdata"))
# load(here::here("outputs","NCP_site_wo_australia.Rdata"))
# NCP_site <- NCP_site_condition

load(here::here("outputs","NN_NS_score_wheighted_mean.Rdata"))
load(here::here("outputs","NN_NS_score_PC1.Rdata"))

##--------------- preping data ---------------
grp <- as.factor(c(recycling_N="NN", recycling_P="NN",Productivity="NS",taxo_richness="NN", funct_entropy="NN",
                   funct_distinctiveness="NN", Selenium_C="NS", Zinc_C="NS", Omega_3_C="NS", Calcium_C="NS",
                   Iron_C="NS", Vitamin_A_C="NS", phylo_entropy="NN", ED_Mean="NN", aesthe_survey="NS", iucn_species="NN",
                   elasmobranch_diversity="NN", low_mg_calcite="NN", high_mg_calcite="NN", aragonite="NN",
                   monohydrocalcite="NN", amorphous_carbonate="NN", biom_lowTL="NN", biom_mediumTL="NN",
                   biom_highTL="NN", fishery_biomass="NS", mean_TL = "NN", robustness = "NN",
                   scientific_interest = "NS", public_interest = "NS")) # /!\ the order matter


NCP_to_transform <- c("Btot","recycling_N","recycling_P","Productivity",
                      "funct_distinctiveness","Omega_3_C","Calcium_C","Vitamin_A_C",
                      "phylo_entropy","ED_Mean", "iucn_species", "elasmobranch_diversity",
                      "low_mg_calcite", "high_mg_calcite", "aragonite", "monohydrocalcite",
                      "amorphous_carbonate", "biom_lowTL", "biom_mediumTL", "biom_highTL",
                      "fishery_biomass")

NCP_site <- NCP_site |>
  dplyr::mutate(across(.cols = all_of(NCP_to_transform),
                       .fns = ~ .x +1 , .names = "{.col}")) |>     
  dplyr::mutate(across(.cols = all_of(NCP_to_transform),
                       .fns = log10 , .names = "{.col}")) |>
  dplyr::left_join(NN_NS_scores) |>  #Add NN and NS scores calculated from Kark formula
  dplyr::left_join(NN_NS_PC1_surveys) |> 
  dplyr::rename(NN_kark = "NN_score", NS_kark = "NS_score")

NCP_site_selected <- subset(NCP_site, select = 
                              c(recycling_N, recycling_P,Productivity,taxo_richness, funct_entropy,
                                funct_distinctiveness, Selenium_C, Zinc_C, Omega_3_C, Calcium_C,
                                Iron_C, Vitamin_A_C, phylo_entropy, ED_Mean, aesthe_survey, iucn_species,
                                elasmobranch_diversity, low_mg_calcite, high_mg_calcite, aragonite,
                                monohydrocalcite, amorphous_carbonate, biom_lowTL, biom_mediumTL,
                                biom_highTL, fishery_biomass, mean_TL , robustness,
                                scientific_interest , public_interest,
                                NN_kark, NS_kark, NN_PC1, NS_PC1)) 

NCP_site_for_pca <- scale(NCP_site_selected)

##-------------Check the Kark 2002 formula-------------
set.seed(6)
dat <- faux::rnorm_multi(n = 1000, 
                   mu = c(0, 10, 20, 25,10),
                   sd = c(1, 5, 5,1,2),
                   r = c(0.9,0.9, 0.7, 0,0.9,0.7,0,0.7,0,0), 
                   varnames = c("A", "B", "C", "D","E"),
                   empirical = FALSE)
plot(dat$A~dat$E)
plot(dat$A~dat$C)
plot(dat$B~dat$C)
plot(dat$D~dat$C)

dat <- scale(dat)
corr_pearson <- stats::cor(dat, method="pearson")

## calculates weighting parameter (Kark 2002):
weighting_par <- c()
for( i in colnames(dat)){
  S<-0
  for( j in colnames(dat)){
    S <- S + (1 - abs(corr_pearson[i,j])/2)
  }
  weighting_par <- c(weighting_par, 1/2 + S)
}

## calculates mean score: Estimator in Dependant Sample (Kark 2002):
EDS <- c()
for( site in 1:nrow(dat)){
  EDS_site <- sum(weighting_par * dat[site,]) / sum(weighting_par)
  EDS <- c(EDS, EDS_site)
}

mean <- rowMeans(dat)

## add mean score:
dat <- cbind(dat, EDS, mean)
pca_test <- FactoMineR::PCA(dat, scale.unit = FALSE, graph=F, ncp=3,
                            quanti.sup = c("EDS", "mean")) 
factoextra::fviz_pca_var(pca_test, repel = TRUE )




####### Test new formula #######
# set.seed(6)
# dat <- faux::rnorm_multi(n = 1000, ## 5 variables totaly correlated
#                          mu = c(0, 10, 20, 25,10),
#                          sd = c(1, 5, 5,1,2),
#                          r = c(.999,.999,.999,.999,.999,.999,.999,.999,.999,.999), 
#                          varnames = c("A", "B", "C", "D","E"),
#                          empirical = FALSE)
# dat <- faux::rnorm_multi(n = 1000, 
#                          mu = c(0, 10, 20, 25,10),
#                          sd = c(1, 5, 5,1,2),
#                          r = c(0.9,0.9, 0.7, 0,0.9,0.7,0,0.7,0,0), 
#                          varnames = c("A", "B", "C", "D","E"),
#                          empirical = FALSE)
# dat <- faux::rnorm_multi(n = 1000, 
#                          mu = c(0, 10, 20, 25,10),
#                          sd = c(1, 5, 5,1,2),
#                          r = c(0.9,0, 0.3, 0,0,0.3,0,-0.5,0,0), 
#                          varnames = c("A", "B", "C", "D","E"),
#                          empirical = FALSE)

dat <- scale(dat)
corr_pearson <- stats::cor(dat, method="pearson")

weighting_par <- c()
for( i in colnames(dat)){
  S<-0
  for( j in colnames(dat)){
    S <- S + (1 - abs(corr_pearson[i,j]))
  }
  weighting_par <- c(weighting_par, S)
}

## calculates mean score: Estimator in Dependant Sample (Kark 2002):
EDS_adapted <- c()
for( site in 1:nrow(dat)){
  EDS_site <- sum(weighting_par * dat[site,]) / sum(weighting_par)
  EDS_adapted <- c(EDS_adapted, EDS_site)
}

mean <- rowMeans(dat)

## add mean score:
dat <- cbind(dat, EDS_adapted)
pca_test <- FactoMineR::PCA(dat, scale.unit = FALSE, graph=F, ncp=3,
                            quanti.sup = c("EDS", "mean", "EDS_adapted")) 

png(filename = here::here("outputs", "figures","test_weighted_mean_formula.png"), 
    width= 20, height = 17, units = "cm", res = 1000)
print( factoextra::fviz_pca_var(pca_test, repel = TRUE ))
dev.off()



##-------------adapted formula from kark 2002-------------
corr_pearson_NCPs <- stats::cor(NCP_site_selected, method="pearson")

## Nature to Nature (NN) score ##
NN_names <- names(grp)[ grp=="NN" ]
corr_pearson_NN <- corr_pearson_NCPs[NN_names,NN_names]
NCP_NN <- NCP_site_for_pca[,NN_names]

# calculates weighting parameter
weighting_par <- c()
for( i in NN_names){
  S<-0
  for( j in NN_names){
    S <- S + (1 - abs(corr_pearson_NN[i,j]) )
  }
  weighting_par <- c(weighting_par, S)
}

## calculates mean score
EDS_NN <- c()
for( site in 1:nrow(NCP_NN)){
  EDS_site <- sum(weighting_par * NCP_NN[site,]) / sum(weighting_par)
  EDS_NN <- c(EDS_NN, EDS_site)
}



## Nature to Society (NS) score ##
NS_names <- names(grp)[ grp=="NS" ]
corr_pearson_NS <- corr_pearson_NCPs[NS_names,NS_names]
NCP_NS <- NCP_site_for_pca[,NS_names]

# calculates weighting parameter
weighting_par <- c()
for( i in NS_names){
  S<-0
  for( j in NS_names){
    S <- S + (1 - abs(corr_pearson_NS[i,j]) )
  }
  weighting_par <- c(weighting_par,  S)
}

# calculates mean score
EDS_NS <- c()
for( site in 1:nrow(NCP_NS)){
  EDS_site <- sum(weighting_par * NCP_NS[site,]) / sum(weighting_par)
  EDS_NS <- c(EDS_NS, EDS_site)
}



NN_NS_scores_adapted_from_kark <- cbind(NCP_site[,c("SiteCode", "SiteLongitude", "SiteLatitude")],
                      data.frame(NN_score = EDS_NN, NS_score = EDS_NS) )
save(NN_NS_scores_adapted_from_kark, file = here::here("outputs", "NN_NS_score_adapted_from_kark.Rdata"))



##-------------computing PCA of nature contributions with NN and NS scores-------------
  NCP_site_for_pca <- cbind(NCP_site_for_pca,
                            NN_mean = rowMeans(NCP_site_for_pca[,names(grp)[ grp=="NN" ]]),
                            NS_mean = rowMeans(NCP_site_for_pca[,names(grp)[ grp=="NS" ]]),
                            NN_adapted_from_kark = NN_NS_scores_adapted_from_kark$NN_score ,
                            NS_adapted_from_kark = NN_NS_scores_adapted_from_kark$NS_score) 
  
  
  pca <- FactoMineR::PCA(NCP_site_for_pca, scale.unit = FALSE, graph=F, ncp=10,
                         quanti.sup = c("NN_kark", "NS_kark", "NN_PC1", "NS_PC1", 
                                        "NN_mean", "NS_mean", "NN_adapted_from_kark",
                                        "NS_adapted_from_kark")) 
  
  NCP_site_kark <- dplyr::select(as.data.frame(NCP_site_for_pca),  -c(NN_PC1, NS_PC1, 
                      NN_mean, NS_mean, NN_adapted_from_kark, NS_adapted_from_kark))
  pca_score <- FactoMineR::PCA(NCP_site_kark, scale.unit = FALSE, graph=F, ncp=10,
                         quanti.sup = c("NN_kark", "NS_kark")) 

  #### PCA in the 2 first dimensions, with representation quality ($cos^{2}$) of each variables
  png(filename = here::here("outputs", "figures","PCA_all_NCP_with_NN_NS_scores.png"), 
      width= 20, height = 17, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = "cos2",
                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                  repel = TRUE 
  ))
  dev.off()
  
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes1-2_with_NN_NS_kark.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca_score, col.var = grp, 
                                  palette = c("forestgreen", "dodgerblue3"),
                                  legend.title = "Nature Based Contributions",
                                  repel = TRUE,
                                  col.quanti.sup = "darkred"))
  dev.off()
  
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes1-2_with_NN_NS_scores.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp, 
                                  palette = c("forestgreen", "dodgerblue3"),
                                  legend.title = "Nature Based Contributions",
                                  repel = TRUE,
                                  col.quanti.sup = "darkred"))
  dev.off()
  
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes3-4_with_NN_NS_kark.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp, 
                                  axe = c(3,4),
                                  palette = c("forestgreen", "dodgerblue3"),
                                  legend.title = "Nature Based Contributions",
                                  repel = TRUE,
                                  col.quanti.sup = "darkred"))
  dev.off()
  
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes3-4_with_NN_NS_scores.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp, 
                                  axe = c(3,4),
                                  palette = c("forestgreen", "dodgerblue3"),
                                  legend.title = "Nature Based Contributions",
                                  repel = TRUE,
                                  col.quanti.sup = "darkred"))
  dev.off()
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes5-6_with_NN_NS_scores.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp, 
                                  axe = c(5,6),
                                  palette = c("forestgreen", "dodgerblue3"),
                                  legend.title = "Nature Based Contributions",
                                  repel = TRUE,
                                  col.quanti.sup = "darkred"))
  dev.off()
  
  