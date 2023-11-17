################################################################################
##
## Compare different methodology to assess composites scores
##
## 1f_test_composite_scores_NP_NN.R
##
## 17/11/2023
##
## Ulysse Flandrin
##
################################################################################


#-----------------Loading packages-------------------
# pkgs <- c("here", "tidyverse", "ggplot2", "sf", "patchwork", "stats", "ggrepel",
#           "ggfun", "scatterpie", "scales", "ggpattern", "ggpubr", "lsr")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
library(ggplot2)
library(patchwork)

rm(list=ls())

##-------------loading data------------
load(here::here("outputs","all_Contrib_site.Rdata"))
load(here::here("outputs","all_Contrib_site_log_transformed.Rdata"))

##-------------Assign contributions into NN and NP categories-------------
grp_NN_NP <- as.factor(c(N_Recycling = "NN",
                         P_Recycling = "NN",
                         Taxonomic_Richness = "NN",
                         Functional_Entropy = "NN", 
                         Phylogenetic_Entropy = "NN",
                         Trait_Distinctiveness = "NN",
                         Evolutionary_Distinctiveness = "NN",
                         Herbivores_Biomass = "NN",
                         Invertivores_Biomass = "NN",
                         Piscivores_Biomass = "NN",
                         Endemism = "NN",
                         Elasmobranch_Diversity = "NN",
                         Low_Mg_Calcite = "NN",
                         High_Mg_Calcite = "NN",
                         Aragonite = "NN",
                         Monohydrocalcite = "NN",
                         Amorphous_Carbonate = "NN",
                         Trophic_Web_Robustness = "NN",
                         Mean_Trophic_Level = "NN",
                         
                         Productivity = "NP",
                         Selenium = "NP",
                         Zinc = "NP",
                         Omega_3 = "NP",
                         Calcium = "NP",
                         Iron = "NP",
                         Vitamin_A = "NP",
                         Fishery_Biomass = "NP",
                         Aesthetic = "NP",
                         Public_Interest = "NP")) # /!\ the order matter

##-------------Clean data-------------
Contrib_log_transformed <- subset(Contrib_site_log_transformed, 
                                  select = -c(SiteCode, SiteCountry, SurveyDate,
                                              SiteEcoregion, SurveyDepth, 
                                              SiteMeanSST, SiteLatitude, SiteLongitude,
                                              Biomass,
                                              HDI, MarineEcosystemDependency,
                                              coral_imputation, gravtot2, mpa_name,
                                              mpa_enforcement, protection_status, 
                                              mpa_iucn_cat))

Contrib_log_scale_clean <- scale(Contrib_log_transformed)

#-------------Compute weighted mean of Contributions-------------
#### Nature to Nature (NN) score ###
NN_names <- names(grp_NN_NP)[ grp_NN_NP=="NN" ]
Contrib_NN <- Contrib_log_scale_clean[,NN_names]

colnames(Contrib_NN)
weighting_par <- c(1/7, 1/7, 1/5, 1/5,1/5, 1/5,1/5,1/5,1/5,1/5,1/5,1/5,
                   1/7,1/7,1/7,1/7,1/7,1/2,1/2)
names(weighting_par) <- colnames(Contrib_NN)
weighting_par

weighted_mean_NN <- c()
for( site in 1:nrow(Contrib_NN)){
  EDS_site <- sum(weighting_par * Contrib_NN[site,]) / sum(weighting_par)
  weighted_mean_NN <- c(weighted_mean_NN, EDS_site) }

summary(weighted_mean_NN)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -1.8112620 -0.3350228  0.0007399  0.0000000  0.3480500  1.3931706


#### Nature to People (NP) score ###
NP_names <- names(grp_NN_NP)[ grp_NN_NP=="NP" ]
Contrib_NP <- Contrib_log_scale_clean[,NP_names]

colnames(Contrib_NP)
weighting_par <- c(1/2, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/2, 1/2,1/2)
names(weighting_par) <- colnames(Contrib_NP)
weighting_par

weighted_mean_NP <- c()
for( site in 1:nrow(Contrib_NP)){
  EDS_site <- sum(weighting_par * Contrib_NP[site,]) / sum(weighting_par)
  weighted_mean_NP <- c(weighted_mean_NP, EDS_site) }

summary(weighted_mean_NP)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.96304 -0.23629  0.05656  0.00000  0.26143  1.68336




#-------------Compute correlated weighted mean of Contributions (Kark 2002)-------------
corr_pearson_Contribs <- stats::cor(Contrib_log_scale_clean, method="pearson")

#### Nature to Nature (NN) score ###
NN_names <- names(grp_NN_NP)[ grp_NN_NP=="NN" ]
corr_pearson_NN <- corr_pearson_Contribs[NN_names,NN_names]
Contrib_NN <- Contrib_log_scale_clean[,NN_names]

## calculates weighting parameter (Kark 2002):
weighting_par <- c()
for( i in NN_names){
  S<-0
  for( j in NN_names){
    S <- S + (1 - abs(corr_pearson_NN[i,j])/2 )
  }
  weighting_par <- c(weighting_par, 1/2 + S)
}

## calculates mean score: Estimator in Dependant Sample (Kark 2002):
EDS_NN <- c()
for( site in 1:nrow(Contrib_NN)){
  EDS_site <- sum(weighting_par * Contrib_NN[site,]) / sum(weighting_par)
  EDS_NN <- c(EDS_NN, EDS_site)
}


#### Nature to People (NP) score ###
NP_names <- names(grp_NN_NP)[ grp_NN_NP=="NP" ]
corr_pearson_NP <- corr_pearson_Contribs[NP_names,NP_names]
Contrib_NP <- Contrib_log_scale_clean[,NP_names]

## calculates weighting parameter (Kark 2002):
weighting_par <- c()
for( i in NP_names){
  S<-0
  for( j in NP_names){
    S <- S + (1 - abs(corr_pearson_NP[i,j])/2 )
  }
  weighting_par <- c(weighting_par, 1/2 + S)
}

## calculates mean score: Estimator in Dependant Sample (Kark 2002):
EDS_NP <- c()
for( site in 1:nrow(Contrib_NP)){
  EDS_site <- sum(weighting_par * Contrib_NP[site,]) / sum(weighting_par)
  EDS_NP <- c(EDS_NP, EDS_site)
}



#-------------Ranking scores-------------
#### Nature to Nature (NN) score ###
NN_names <- names(grp_NN_NP)[ grp_NN_NP=="NN" ]
Contrib_NN <- as.data.frame(Contrib_log_scale_clean[,NN_names])
ranked_contrib_NN <- Contrib_NN |>  
  dplyr::mutate(across(.cols= everything(), .fns = rank))

colnames(ranked_contrib_NN)
weighting_par <- c(1/7, 1/7, 1/5, 1/5,1/5, 1/5,1/5,1/5,1/5,1/5,1/5,1/5,
                   1/7,1/7,1/7,1/7,1/7,1/2,1/2)
names(weighting_par) <- colnames(ranked_contrib_NN)
weighting_par

ranked_mean_NN <- c()
for( site in 1:nrow(ranked_contrib_NN)){
  mean_rank <- sum(weighting_par * ranked_contrib_NN[site,]) / sum(weighting_par)
  ranked_mean_NN <- c(ranked_mean_NN, mean_rank) }

summary(ranked_mean_NN)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# 144.0   484.5       616.5        619.0     757.0       1058.2 
ranked_mean_NN <- scale(ranked_mean_NN)


#### Nature to People (NP) score ###
NP_names <- names(grp_NN_NP)[ grp_NN_NP=="NP" ]
Contrib_NP <- as.data.frame(Contrib_log_scale_clean[,NP_names])

ranked_contrib_NP <- Contrib_NP |>  
  dplyr::mutate(across(.cols= everything(), .fns = rank))

colnames(ranked_contrib_NP)
weighting_par <- c(1/2, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/2, 1/2,1/2)
names(weighting_par) <- colnames(ranked_contrib_NP)
weighting_par

ranked_mean_NP <- c()
for( site in 1:nrow(ranked_contrib_NP)){
  mean_rank <- sum(weighting_par * ranked_contrib_NP[site,]) / sum(weighting_par)
  ranked_mean_NP <- c(ranked_mean_NP, mean_rank) }

summary(ranked_mean_NP)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 85.61  514.06  635.11  619.00  727.14  981.36 
ranked_mean_NP <- scale(ranked_mean_NP)



#-------------Compute unweighted mean of Contributions-------------
unweighted_mean_NN <- rowMeans(Contrib_NN)

unweighted_mean_NP <- rowMeans(Contrib_NP)


#-------------TOPSIS Scores (Jouval 2023)-------------
TOPSIS_score <- function(matrix){
  TOPSIS <- c()
  for(i in 1:nrow(matrix)){
    #distance to the worse: Di-
    square_dist <- 0
    for(j in colnames(matrix)){
      dist_j <- as.numeric((matrix[i,j] - min(matrix[,j]))^2)
      square_dist <- square_dist + dist_j
    }
    dist_worst <- sqrt(square_dist)
    
    #distance to the best: Di+
    square_dist <- 0
    for(j in colnames(matrix)){
      dist_j <- as.numeric((matrix[i,j] - max(matrix[,j]))^2)
      square_dist <- square_dist + dist_j
    }
    dist_best <- sqrt(square_dist)
    
    ## Proximity score to the best
    RI_i <- dist_worst / (dist_worst + dist_best)
    TOPSIS <- c(TOPSIS, RI_i)
  }
  return(TOPSIS)
} # end of fct TOPSIS_score

#### Nature to Nature (NN) score ###
NN_names <- names(grp_NN_NP)[ grp_NN_NP=="NN" ]
Contrib_NN <- Contrib_log_scale_clean[,NN_names]

colnames(Contrib_NN)
weighting_par <- c(1/7, 1/7, 1/5, 1/5,1/5, 1/5,1/5,1/5,1/5,1/5,1/5,1/5,
                   1/7,1/7,1/7,1/7,1/7,1/2,1/2)
names(weighting_par) <- colnames(Contrib_NN)
weighting_par

Contrib_NN_weighted <- Contrib_NN * weighting_par

TOPSIS_NN <- TOPSIS_score(Contrib_NN_weighted)


#### Nature to People (NP) score ###
NP_names <- names(grp_NN_NP)[ grp_NN_NP=="NP" ]
Contrib_NP <- Contrib_log_scale_clean[,NP_names]

colnames(Contrib_NP)
weighting_par <- c(1/2, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/2, 1/2,1/2)
names(weighting_par) <- colnames(Contrib_NP)
weighting_par

Contrib_NP_weighted <- Contrib_NP * weighting_par

TOPSIS_NP <- TOPSIS_score(Contrib_NP_weighted)




#-------------Compute geometric mean of Contributions (Nardo 2005)-------------
Contrib_log_transformed
Contrib_0_to_1 <- Contrib_log_transformed |> 
  dplyr::mutate(across(.cols = everything(), .fns = scales::rescale))

#### Nature to Nature (NN) score ###
NN_names <- names(grp_NN_NP)[ grp_NN_NP=="NN" ]
Contrib_NN <- Contrib_0_to_1[,NN_names] +1 ### all values range from 1 to 2 (Avoid strange product)

colnames(Contrib_NN)
weighting_par <- c(1/7, 1/7, 1/5, 1/5,1/5, 1/5,1/5,1/5,1/5,1/5,1/5,1/5,
                   1/7,1/7,1/7,1/7,1/7,1/2,1/2)

geometric_aggregation_NN <- c()
for( site in 1:nrow(Contrib_NN)){
  geom_site <- prod(Contrib_NN[site,]^weighting_par )
  geometric_aggregation_NN <- c(geometric_aggregation_NN, geom_site) }

hist(geometric_aggregation_NN)



#### Nature to People (NP) score ###
NP_names <- names(grp_NN_NP)[ grp_NN_NP=="NP" ]
Contrib_NP <- Contrib_0_to_1[,NP_names] +1
  
colnames(Contrib_NP)
weighting_par <- c(1/2, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/2, 1/2,1/2)

geometric_aggregation_NP <- c()
for( site in 1:nrow(Contrib_NP)){
  geom_site <- prod(Contrib_NP[site,]^weighting_par )
  geometric_aggregation_NP <- c(geometric_aggregation_NP, geom_site) }

hist(geometric_aggregation_NP)




#-------------Compare all scores-------------
scores <- cbind(weighted_mean_NN, weighted_mean_NP,
                EDS_NN, EDS_NP,
                ranked_mean_NN, ranked_mean_NP,
                unweighted_mean_NN, unweighted_mean_NP,
                TOPSIS_NN, TOPSIS_NP,
                geometric_aggregation_NN, geometric_aggregation_NP)
colnames(scores) <- c("weighted_mean_NN", "weighted_mean_NP",
                      "EDS_NN", "EDS_NP",
                      "ranked_mean_NN", "ranked_mean_NP",
                      "unweighted_mean_NN", "unweighted_mean_NP",
                      "TOPSIS_NN", "TOPSIS_NP",
                      "geometric_aggregation_NN", "geometric_aggregation_NP")


M <- cor(scores)

# png(filename = here::here("outputs", "figures","corr_matrix_log_transformed_Contributions.png"), 
#     width= 40, height = 30, units = "cm", res = 1000)
# print({
corrplot(M, order = 'hclust', hclust.method = 'ward.D2', addrect = 2, addCoef.col = 'grey')
# })
# dev.off() 





### Save NN and NP scores ###
NN_NP_scores <- cbind(Contrib_site_log_transformed[ , c("SiteCode", "SiteLongitude", "SiteLatitude", 
                                                        "SiteCountry", "SiteEcoregion","year",
                                                        "SiteMeanSST", "Biomass", "HDI", "gravtot2",
                                                        "MarineEcosystemDependency", "coral_imputation",
                                                        "mpa_name", "mpa_enforcement", "protection_status",
                                                        "mpa_iucn_cat")],
                      data.frame(NN_score = mean_NN, NP_score = mean_NP) )
save(NN_NP_scores, file = here::here("outputs", "NN_NP_score_wheighted_mean.Rdata"))
