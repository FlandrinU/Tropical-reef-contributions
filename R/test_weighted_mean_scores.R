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

##-------------loading data------------
load(here::here("outputs","all_NCP_site_log_transformed.Rdata"))
# load(here::here("outputs","NCP_site_log_coral_reef.Rdata"))
# load(here::here("outputs","NCP_site_log_SST20.Rdata"))
# load(here::here("outputs","NCP_site_log_coral_5_imputed.Rdata"))
# load(here::here("outputs","NCP_site_log_wo_australia.Rdata"))
# NCP_site_log_transformed <- NCP_site_condition


load(here::here("outputs","NN_NS_score_wheighted_mean.Rdata"))
load(here::here("outputs","NN_NS_score_wheighted_mean_quarter_ranked.Rdata"))

load(here::here("outputs","NN_NS_score_PC1.Rdata"))
##--------------- preping data ---------------
NCP_site_log <- NCP_site_log_transformed |>
  dplyr::left_join(NN_NS_with_product) |>  #Add NN and NS scores calculated from Kark formula, and the ranks in each quarter
  dplyr::left_join(NN_NS_PC1_surveys) |> 
  dplyr::rename(NN_kark = "NN_score", NS_kark = "NS_score") 

NCP_site_selected <- subset(NCP_site_log, select = 
                              c(SiteCode,
                                N_recycling,P_recycling,Taxonomic_Richness, Functional_Entropy,
                                Phylogenetic_Entropy, Functional_Distinctiveness,Evolutionary_distinctiveness,
                                Low_TL_Biomass, Medium_TL_Biomass, High_TL_Biomass, Endemism, #IUCN_Species,
                                Elasmobranch_Diversity, Low_Mg_Calcite, High_Mg_Calcite, Aragonite, 
                                Monohydrocalcite, Amorphous_Carbonate, Trophic_web_robustness,
                                mean_Trophic_Level, 
                                Productivity,Selenium,Zinc,Omega_3,Calcium,Iron,Vitamin_A,
                                Fishery_Biomass, Aesthetic,Public_Interest,#Academic_Knowledge,
                                NN_kark, NS_kark, NN_PC1, NS_PC1)) |>
  tibble::column_to_rownames(var = "SiteCode")

NCP_site_for_pca <- scale(NCP_site_selected)


##-------------NCP in categories-------------
## Classify variables in Nature for Nature (NN) and Nature for Society (NS)
grp_NN_NS <- as.factor(c(N_recycling = "NN",
                         P_recycling = "NN",
                         Taxonomic_Richness = "NN",
                         Functional_Entropy = "NN", 
                         Phylogenetic_Entropy = "NN",
                         Functional_Distinctiveness = "NN",
                         Evolutionary_distinctiveness = "NN",
                         Low_TL_Biomass = "NN",
                         Medium_TL_Biomass = "NN",
                         High_TL_Biomass = "NN",
                         #IUCN_Species = "NN",
                         Endemism = "NN", 
                         Elasmobranch_Diversity = "NN",
                         Low_Mg_Calcite = "NN",
                         High_Mg_Calcite = "NN",
                         Aragonite = "NN",
                         Monohydrocalcite = "NN",
                         Amorphous_Carbonate = "NN",
                         Trophic_web_robustness = "NN",
                         mean_Trophic_Level = "NN",
                         
                         Productivity = "NS",
                         Selenium = "NS",
                         Zinc = "NS",
                         Omega_3 = "NS",
                         Calcium = "NS",
                         Iron = "NS",
                         Vitamin_A = "NS",
                         Fishery_Biomass = "NS",
                         Aesthetic = "NS",
                         Public_Interest = "NS")) # /!\ the order matter

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
NN_names <- names(grp_NN_NS)[ grp_NN_NS=="NN" ]
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
NS_names <- names(grp_NN_NS)[ grp_NN_NS=="NS" ]
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



NN_NS_scores_adapted_from_kark <- cbind(NCP_site_log_transformed[,c("SiteCode", "SiteLongitude", "SiteLatitude")],
                      data.frame(NN_score = EDS_NN, NS_score = EDS_NS) )
save(NN_NS_scores_adapted_from_kark, file = here::here("outputs", "NN_NS_score_adapted_from_kark.Rdata"))



##-------------computing PCA of nature contributions with NN and NS scores-------------
  NCP_site_for_pca <- cbind(NCP_site_for_pca,
                            NN_mean = rowMeans(NCP_site_for_pca[,names(grp_NN_NS)[ grp_NN_NS=="NN" ]]),
                            NS_mean = rowMeans(NCP_site_for_pca[,names(grp_NN_NS)[ grp_NN_NS=="NS" ]]),
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
  print( factoextra::fviz_pca_var(pca_score, col.var = grp_NN_NS, 
                                  palette = c("forestgreen", "dodgerblue3"),
                                  legend.title = "Nature Based Contributions",
                                  repel = TRUE,
                                  col.quanti.sup = "darkred"))
  dev.off()
  
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes1-2_with_NN_NS_scores.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp_NN_NS, 
                                  palette = c("forestgreen", "dodgerblue3"),
                                  legend.title = "Nature Based Contributions",
                                  repel = TRUE,
                                  col.quanti.sup = "darkred"))
  dev.off()
  
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes3-4_with_NN_NS_kark.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp_NN_NS, 
                                  axe = c(3,4),
                                  palette = c("forestgreen", "dodgerblue3"),
                                  legend.title = "Nature Based Contributions",
                                  repel = TRUE,
                                  col.quanti.sup = "darkred"))
  dev.off()
  
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes3-4_with_NN_NS_scores.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp_NN_NS, 
                                  axe = c(3,4),
                                  palette = c("forestgreen", "dodgerblue3"),
                                  legend.title = "Nature Based Contributions",
                                  repel = TRUE,
                                  col.quanti.sup = "darkred"))
  dev.off()
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes5-6_with_NN_NS_scores.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp_NN_NS, 
                                  axe = c(5,6),
                                  palette = c("forestgreen", "dodgerblue3"),
                                  legend.title = "Nature Based Contributions",
                                  repel = TRUE,
                                  col.quanti.sup = "darkred"))
  dev.off()
  
  
  
  ##-------------Highlighting outliers in PCA-------------
  #NN only
  ind_NN <- which(NCP_site_log$rank_d_r >= quantile(NCP_site_log$rank_d_r, probs=c(0.95), na.rm=T))
  #NS only
  ind_NS <- which(NCP_site_log$rank_u_l >= quantile(NCP_site_log$rank_u_l, probs=c(0.95), na.rm=T))
  #NN and NS
  ind_NNxNS <- which(NCP_site_log$rank_u_r >= quantile(NCP_site_log$rank_u_r, probs=c(0.95), na.rm=T))
  #neither NN nor NS
  ind_worse <- which(NCP_site_log$rank_d_l >= quantile(NCP_site_log$rank_d_l, probs=c(0.95), na.rm=T))

  ## plot on PCA with score outliers
  library(ggplot2)
  ind_coord <- as.data.frame(factoextra::get_pca_ind(pca_score)[["coord"]]) |>
    tibble::rownames_to_column(var = "label")
  
  factoextra::fviz_pca_biplot(pca_score, col.var = grp_NN_NS, 
                              palette = c("forestgreen", "dodgerblue3"),
                              legend.title = "Nature Based Contributions",
                              alpha.var = 0.3, alpha.ind = 0.5,
                              repel = TRUE,
                              geom=c("point", "text"), pointshape=21,
                              col.quanti.sup = "grey20",
                              label = c("var", "quanti.sup")) + 
    
    geom_point(data = ind_coord[ind_NN,], aes(x = Dim.1, y = Dim.2),
               size = 2.5,col = "darkgreen") +
    geom_point(data = ind_coord[ind_NS,], aes(x = Dim.1, y = Dim.2),
               size = 2.5,col = "dodgerblue4") +
    geom_point(data = ind_coord[ind_NNxNS,], aes(x = Dim.1, y = Dim.2),
               size = 2.5,col = "darkred") +
    geom_point(data = ind_coord[ind_worse,], aes(x = Dim.1, y = Dim.2),
               size = 2.5,col = "grey20")+
    geom_text(data = ind_coord[c(ind_NN, ind_NS, ind_NNxNS, ind_worse),],
              aes(x= Dim.1, y = Dim.2,label=label),
              hjust=-0.1, vjust=-0.1, size = 3)
  
  ggsave(filename = here::here("outputs", "figures", "PCA_with_quantiles0.95.jpg"))    
  
