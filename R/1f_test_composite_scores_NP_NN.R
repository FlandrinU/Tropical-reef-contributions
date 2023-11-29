################################################################################
##
## Compare different methodology to assess composites indicators NN and NP
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
#           "ggfun", "scatterpie", "scales", "ggpattern", "ggpubr", "lsr", "CompInd")
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
                         
                         Turnover_Available_Biomass = "NP",
                         Selenium = "NP",
                         Zinc = "NP",
                         Omega_3 = "NP",
                         Calcium = "NP",
                         Iron = "NP",
                         Vitamin_A = "NP",
                         Available_Biomass = "NP",
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

Contrib_0_to_1 <- Contrib_log_transformed |> 
  dplyr::mutate(across(.cols = everything(), .fns = scales::rescale))

## NN contributions
NN_names <- names(grp_NN_NP)[ grp_NN_NP=="NN" ]
Contrib_NN <- Contrib_log_scale_clean[,NN_names]

## NP contributions
NP_names <- names(grp_NN_NP)[ grp_NN_NP=="NP" ]
Contrib_NP <- Contrib_log_scale_clean[,NP_names]

#-------------Compute weighted arithmetic mean of Contributions-------------
#### Nature to Nature (NN) score ###
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
hist(weighted_mean_NN) #normal distribution


#### Nature to People (NP) score ###
colnames(Contrib_NP)
weighting_par <- c(1/2, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/2, 1/2,1/2)
names(weighting_par) <- colnames(Contrib_NP)
weighting_par

weighted_mean_NP <- c()
for( site in 1:nrow(Contrib_NP)){
  EDS_site <- sum(weighting_par * Contrib_NP[site,]) / sum(weighting_par)
  weighted_mean_NP <- c(weighted_mean_NP, EDS_site) }

summary(weighted_mean_NP)
hist(weighted_mean_NP, breaks = 20) #normal distribution



#-------------Compute correlated weighted mean of Contributions (Kark 2002)-------------
corr_pearson_Contribs <- stats::cor(Contrib_log_scale_clean, method="pearson")

#### Nature to Nature (NN) score ###
corr_pearson_NN <- corr_pearson_Contribs[NN_names,NN_names]

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
corr_pearson_NP <- corr_pearson_Contribs[NP_names,NP_names]

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

hist(EDS_NN) #normal distribution
hist(EDS_NP, breaks = 20) #normal distribution

#-------------Ranking scores-------------
#### Nature to Nature (NN) score ###
ranked_contrib_NN <- as.data.frame(Contrib_log_scale_clean[,NN_names]) |>  
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
ranked_mean_NN <- scale(ranked_mean_NN)
hist(ranked_mean_NN) #normal distribution
plot(weighted_mean_NN~ranked_mean_NN)


#### Nature to People (NP) score ###
ranked_contrib_NP <- as.data.frame(Contrib_log_scale_clean[,NP_names]) |>  
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
ranked_mean_NP <- scale(ranked_mean_NP)
hist(ranked_mean_NP) #normal distribution
plot(weighted_mean_NP~ranked_mean_NP)



#-------------Compute unweighted mean of Contributions-------------
unweighted_mean_NN <- rowMeans(Contrib_NN)

unweighted_mean_NP <- rowMeans(Contrib_NP)

hist(unweighted_mean_NN)
hist(unweighted_mean_NP, breaks = 20)


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
} # end of function TOPSIS_score

#### Nature to Nature (NN) score ###
Contrib_NN_1_2 <- Contrib_0_to_1[,NN_names] +1 ### all values range from 1 to 2 (Avoid strange product)


TOPSIS_NN <- TOPSIS_score(Contrib_NN_1_2)
hist(TOPSIS_NN)
plot(weighted_mean_NN~TOPSIS_NN)

#### Nature to People (NP) score ###
Contrib_NP_1_2 <- Contrib_0_to_1[,NP_names] +1 ### all values range from 1 to 2 (Avoid strange product)

TOPSIS_NP <- TOPSIS_score(Contrib_NP_1_2)
hist(TOPSIS_NP, breaks = 30)
plot(weighted_mean_NP~TOPSIS_NP)



#-------------Compute geometric mean of Contributions (Nardo 2005)-------------

#### Nature to Nature (NN) score ###
NN_names <- names(grp_NN_NP)[ grp_NN_NP=="NN" ]
Contrib_NN <- Contrib_0_to_1[,NN_names] +1 ### all values range from 1 to 2 (Avoid strange product)

# weighting_par <- c(1/7, 1/7, 1/5, 1/5,1/5, 1/5,1/5,1/5,1/5,1/5,1/5,1/5,
#                    1/7,1/7,1/7,1/7,1/7,1/2,1/2)
weighting_par <- rep(1,19)

geometric_aggregation_NN <- c()
for( site in 1:nrow(Contrib_NN)){
  geom_site <- prod(Contrib_NN[site,]^(weighting_par/ sum(weighting_par)) ) # sum of weights is equal to 1
  geometric_aggregation_NN <- c(geometric_aggregation_NN, geom_site) }

hist(geometric_aggregation_NN, breaks = 30)



#### Nature to People (NP) score ###
NP_names <- names(grp_NN_NP)[ grp_NN_NP=="NP" ]
Contrib_NP <- Contrib_0_to_1[,NP_names] +1
  
# weighting_par <- c(1/2, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/2, 1/2,1/2)
weighting_par <- rep(1,10)

geometric_aggregation_NP <- c()

for( site in 1:nrow(Contrib_NP)){
  geom_site <- prod(Contrib_NP[site,]^(weighting_par/ sum(weighting_par)) )
  geometric_aggregation_NP <- c(geometric_aggregation_NP, geom_site) }

hist(geometric_aggregation_NP, breaks = 30)


#-------------others Composite indicators, package Compind (Vidoli and Fusco 2023) -------------
#### Nature to Nature (NN) score ###
Contrib_NN_01 <- as.data.frame(Contrib_0_to_1[,NN_names])

## Benefit of the Doubt approach : BoD composite score depends exclusively on the frontier’s distance
BoD_NN = Compind::ci_bod(Contrib_NN_01) 
hist(BoD_NN$ci_bod_est, breaks = 30) 
#on normal distribution + indicator not very used + many weights of 0 -> do not use this indicator

## Weighting method based on Factor Analysis: the composite Indicator estimated values are equal to component scores multiplied by its proportion variance
data_norm=Compind::normalise_ci(Contrib_NN_01,
                                c(1:ncol(Contrib_NN)),
                                polarity = rep("POS", 19), 
                                method=2) 
Factor_analysis_NN = Compind::ci_factor(data_norm$ci_norm,c(1:5),method="ALL")
hist(Factor_analysis_NN$ci_factor_est)
summary(Factor_analysis_NN$ci_factor_est)
plot(Factor_analysis_NN$ci_factor_est ~ weighted_mean_NN)

## Mazziotta-Pareto Index (MPI) method: non-compensative composite index which, starting from a linear aggregation, introduces a penalty for the units with unbalanced values of the indicators (De Muro et al., 2010)
data_norm = Compind::normalise_ci(Contrib_NN_01,
                                  c(1:ncol(Contrib_NN)),
                                  rep("POS", 19), 
                                  method=1,
                                  z.mean=0, 
                                  z.std=1)
MPI_NN = Compind::ci_mpi(data_norm$ci_norm, penalty="NEG")
hist(MPI_NN[["ci_mpi_est"]])
summary(MPI_NN[["ci_mpi_est"]])



#### Nature for People (NP) score ###
Contrib_NP_01 <- as.data.frame(Contrib_0_to_1[,NP_names])

## Benefit of the Doubt approach : BoD composite score depends exclusively on the frontier’s distance
BoD_NP = Compind::ci_bod(Contrib_NP_01) 
hist(BoD_NP$ci_bod_est, breaks = 30) #non normal distribution...

## Weighting method based on Factor Analysis: the composite Indicator estimated values are equal to component scores multiplied by its proportion variance
data_norm=Compind::normalise_ci(Contrib_NP_01,
                                c(1:ncol(Contrib_NP)),
                                polarity = rep("POS", 19), 
                                method=2) 
Factor_analysis_NP = Compind::ci_factor(data_norm$ci_norm,c(1:5),method="ALL")
hist(Factor_analysis_NP$ci_factor_est, breaks = 20)
summary(Factor_analysis_NP$ci_factor_est) #long tail
plot(Factor_analysis_NP$ci_factor_est ~ weighted_mean_NP) #Many contributions are negatively correlated: we lose the directionality of contributions using PCA axis.


## Mazziotta-Pareto Index (MPI) method: non-compensative composite index which, starting from a linear aggregation, introduces a penalty for the units with unbalanced values of the indicators (De Muro et al., 2010)
data_norm = Compind::normalise_ci(Contrib_NP_01,
                                  c(1:ncol(Contrib_NP)),
                                  rep("POS", 10), 
                                  method=1,
                                  z.mean=0, 
                                  z.std=1)
MPI_NP = Compind::ci_mpi(data_norm$ci_norm, penalty="NEG")
hist(MPI_NP[["ci_mpi_est"]], breaks = 20) #long tail...
summary(MPI_NP[["ci_mpi_est"]])

plot(MPI_NP[["ci_mpi_est"]] ~ weighted_mean_NP)


#-------------Compare all scores-------------
scores <- cbind(weighted_mean_NN,
                unweighted_mean_NN,
                EDS_NN, 
                ranked_mean_NN,
                geometric_aggregation_NN,
                TOPSIS_NN, 
                MPI_NN[["ci_mpi_est"]],
                
                weighted_mean_NP,
                unweighted_mean_NP,
                EDS_NP,
                ranked_mean_NP,
                geometric_aggregation_NP,
                TOPSIS_NP,
                MPI_NP[["ci_mpi_est"]])

                # BoD_NN$ci_bod_est, BoD_NP$ci_bod_est, 
                # Factor_analysis_NN$ci_factor_est, Factor_analysis_NP$ci_factor_est,
                # -CI_wroclaw_estimated[["ci_wroclaw_est"]],
                # CI_geom_estimated[["ci_mean_geom_est"]],
                # CI_r1$ci_rbod_est)

colnames(scores) <- c("weighted_mean_NN",
                      "unweighted_mean_NN",
                      "EDS_NN",
                      "ranked_mean_NN", 
                      "geometric_aggregation_NN", 
                      "TOPSIS_NN", 
                      "MPI_NN",
                      
                      "weighted_mean_NP",
                      "unweighted_mean_NP",
                      "EDS_NP",
                      "ranked_mean_NP",
                      "geometric_aggregation_NP",
                      "TOPSIS_NP",
                      "MPI_NP")

                      # "BoD_NN","BoD_NP",
                      # "Factor_analysis_NN", "Factor_analysis_NP", 
                      # "CI_wroclaw_estimated", "CI_geom_estimated",
                      # "robust_DoB")



png(filename = here::here("outputs", "figures","corr_matrix_composite_indicators.png"),
    width= 30, height = 23, units = "cm", res = 1000)
print({
  M <- cor(scores)
  corrplot::corrplot(M, order = 'original', tl.col = "black",
                   addCoef.col = 'grey60', tl.srt = 45)
})
dev.off()



#-------------Plot scores compared to arithmetic mean-------------
scatterplot_function <- function(data,
                                 x="weighted_mean_NN",
                                 y="EDS_NN"){
  ggplot(data, 
         aes_string(x=x, y=y, colour = x)) +
    geom_point()+
    theme_bw()+
    scale_colour_gradientn(name  ="Weighted arithmetic mean",
                           colours = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu")))+
    xlab(paste( x))+ ylab(paste( y)) +
    guides(colour = guide_colourbar(title.position="top")) +
    theme(legend.position = "none",
          legend.key.size = unit(0.5, 'cm'),
          legend.direction = "vertical",
          legend.title = element_text( size = 11),
          legend.background = element_rect(fill='transparent'),
          axis.title=element_text(size=11)) +
    geom_smooth(method="lm", color="#757575",linetype = "dashed",
                linewidth = 1)+
    ggpubr::stat_cor(method = "pearson")
  # +
  #   ggrepel::geom_label_repel(
  #     data=dplyr::filter(data, abs(x - y) >
  #                          quantile(abs(x - y), 0.98) ),
  #     aes(label= paste(Var1, "-", Var2)),
  #     size=2, fill = "white", 
  #     min.segment.length = 0.1,
  #     color = "black", alpha = 0.8,
  #     direction = "both",
  #     seed = 1968)
  
}

plots <- lapply(colnames(scores)[grep("NN",colnames(scores))], function(i){
  scatterplot_function(as.data.frame(scores),
                       x="weighted_mean_NN",
                       y=i)
})

plot_NN <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + 
  plots[[6]] + plots[[7]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

plot_NN


plots <- lapply(colnames(scores)[grep("NP",colnames(scores))], function(i){
  scatterplot_function(as.data.frame(scores),x="weighted_mean_NP",y=i) })

plot_NP <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + 
  plots[[6]] + plots[[7]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

plot_NP
