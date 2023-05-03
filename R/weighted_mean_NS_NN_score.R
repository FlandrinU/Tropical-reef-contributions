################################################################################
##
## Calculate NN and NS scores of sites
##
## weighted_mean_NS_NN_score.R
##
## 17/01/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse", "ggplot2", "sf", "patchwork", "stats", "ggrepel",
          "ggfun", "scatterpie", "scales", "ggpattern", "ggpubr", "lsr")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data------------
load(here::here("outputs","all_NCP_site.Rdata"))
load(here::here("outputs","all_NCP_site_log_transformed.Rdata"))
# load(here::here("outputs","NCP_site_log_coral_reef.Rdata"))
# load(here::here("outputs","NCP_site_log_SST20.Rdata"))
# load(here::here("outputs","NCP_site_log_coral_5_imputed.Rdata"))
# load(here::here("outputs","NCP_site_log_wo_australia.Rdata"))
# load(here::here("outputs","NCP_site_log_random.Rdata"))
# load(here::here("outputs","NCP_site_log_only_australia.Rdata"))
# NCP_site_log_transformed <- NCP_site_condition


load(here::here("data","metadata_surveys.Rdata"))
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))
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
                         # IUCN_Species = "NN",
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

##-------------compute correlation coefficients-------------

## Clean data
NCP_log_transformed <- subset(NCP_site_log_transformed, 
                              select = -c(SiteCode, SiteCountry, SurveyDate,
                                          SiteEcoregion, SurveyDepth, 
                                          SiteMeanSST, SiteLatitude, SiteLongitude,
                                          Biomass,
                                          HDI, MarineEcosystemDependency,
                                          coral_imputation, gravtot2, mpa_name,
                                          mpa_enforcement, protection_status, 
                                          mpa_iucn_cat))

#### NCPs distribution
plot_distribution <- function(ncp, data){
  col <- fishualize::fish(n = ncol(data), option = "Ostracion_whitleyi", begin = 0, end = 0.8)
  names(col) <- colnames(data)
  
  ggplot(data) +
    aes(x = data[,ncp][[ncp]]) +
    geom_histogram(bins = 40L,
                   fill = col[ncp][[1]],
                   col = "black") +
    labs(title = ncp) +
    xlab("") + ylab("") +
    theme_minimal() +
    theme( legend.position = "none")
}
# 
# plots <- lapply(colnames(NCP_site_clean), FUN = plot_distribution, data = NCP_site_clean )
# 
# library(patchwork)
# all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
#   plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
#   plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
#   plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
#   plots[[28]] + plots[[29]] +  plots[[30]] + plots[[31]] + 
#   theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
#   plot_annotation(tag_levels = "a") &
#   theme(plot.tag = element_text(face = 'bold'))
# 
# # ggsave(filename = here("outputs", "figures","NCP_distribution.png"), all_plot, width = 22, height =14 )

#### NCPs distribution with log correction
plots <- lapply(colnames(NCP_log_transformed), FUN = plot_distribution, data = NCP_log_transformed )
library(patchwork)
all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + 
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

all_plot
# ggsave(filename = here("outputs", "figures","NCP_log_transformed_distribution.png"), all_plot, width = 22, height =14 )

### Correlation between NCPs
correlations_between_NCPs <- stats::cor(NCP_log_transformed, method="pearson")

#### Corr-matrix for all NCPs
corrplot::corrplot(correlations_between_NCPs, order = 'AOE', type = 'lower', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r')
corrplot::corrplot(correlations_between_NCPs, add = TRUE, type = 'upper', method = 'number', order = 'AOE', insig = 'p-value',
                   diag = FALSE, tl.pos = 'n', cl.pos = 'n')



##-------------compute weighted mean of NCP (Kark 2002): NN and NS scores-------------
corr_pearson_NCPs <- stats::cor(NCP_log_transformed, method="pearson")
NCP_log_scale_clean <- scale(NCP_log_transformed)

#### Nature to Nature (NN) score ####
NN_names <- names(grp_NN_NS)[ grp_NN_NS=="NN" ]
corr_pearson_NN <- corr_pearson_NCPs[NN_names,NN_names]
NCP_NN <- NCP_log_scale_clean[,NN_names]

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
for( site in 1:nrow(NCP_NN)){
  EDS_site <- sum(weighting_par * NCP_NN[site,]) / sum(weighting_par)
  EDS_NN <- c(EDS_NN, EDS_site)
}



#### Nature to Society (NS) score ####
NS_names <- names(grp_NN_NS)[ grp_NN_NS=="NS" ]
corr_pearson_NS <- corr_pearson_NCPs[NS_names,NS_names]
NCP_NS <- NCP_log_scale_clean[,NS_names]

## calculates weighting parameter (Kark 2002):
weighting_par <- c()
for( i in NS_names){
  S<-0
  for( j in NS_names){
    S <- S + (1 - abs(corr_pearson_NS[i,j])/2 )
  }
  weighting_par <- c(weighting_par, 1/2 + S)
}

## calculates mean score: Estimator in Dependant Sample (Kark 2002):
EDS_NS <- c()
for( site in 1:nrow(NCP_NS)){
  EDS_site <- sum(weighting_par * NCP_NS[site,]) / sum(weighting_par)
  EDS_NS <- c(EDS_NS, EDS_site)
}




NN_NS_scores <- cbind(NCP_site_log_transformed[ , c("SiteCode", "SiteLongitude", "SiteLatitude", 
                                                    "SiteCountry", "SiteEcoregion",
                                                    "SiteMeanSST", "Biomass", "HDI", "gravtot2",
                                                    "MarineEcosystemDependency", "coral_imputation",
                                                    "mpa_name", "mpa_enforcement", "protection_status",
                                                    "mpa_iucn_cat")],
  data.frame(NN_score = EDS_NN, NS_score = EDS_NS) )
save(NN_NS_scores, file = here::here("outputs", "NN_NS_score_wheighted_mean.Rdata"))

##-------------plot NN and NS scores-------------
library(ggplot2)
library(patchwork)


plot_histogram <- function(data = NN_NS_scores, 
                           x=NN_NS_scores$NN_score, col= "forestgreen",
                           title="Nature for nature score site's evaluation"){
  ggplot(data, aes(x) ) +
    geom_histogram( bins = 40, color="grey60", fill= "white", linewidth=0.2)+
    
    ggnewscale::new_scale_fill() +
    geom_histogram( bins = 40, color="grey60", aes(fill = after_stat(x)), 
                    linewidth=0.2)+
    scale_fill_gradient2(low="white", mid="white", high = col ,
                         midpoint= quantile(x,0.75), 
                         limits =quantile(x, probs=c( 0.75,1))+ c(-0.1,0.1),
                         na.value="transparent")+
    
    ggnewscale::new_scale_fill() +
    geom_histogram( bins = 40, aes(fill=after_stat(x)), color="grey60", linewidth=0.2)+
    scale_fill_gradient2(low="black", mid="white", high="white",
                         midpoint= quantile(x,0.25), 
                         limits =quantile(x, probs=c(0,0.25)) + c(-0.1,0.1),
                         na.value="transparent")+
  
    labs(title = title)+
    geom_vline(xintercept = quantile(x, 0.25), linetype = 3, col= "black")+
    geom_vline(xintercept = quantile(x, 0.75), linetype = 3, col= "black")+
    theme_minimal()+
    theme(legend.position = "none")
} #end of plot_histogram function

## NN
# histo_NN <- ggplot(NN_NS_scores, aes(x = NN_score)) +
#   geom_histogram( bins = 40, aes(fill=after_stat(x)), color="grey50", linewidth=0.3)+
#   scale_fill_gradient2(low = "black", mid = "white", high = "forestgreen", 
#                        midpoint = median(NN_NS_scores$NN_score),
#                        guide = "colourbar")+
#   labs(x = "Nature for nature score site's evaluation")+
#   geom_vline(xintercept = 0, linetype = 2, col= "black")+
#   theme_minimal()
# histo_NN
# ggsave(filename = here("outputs", "figures","hist_NN_weighted_mean.png"), histo_NN, width = 8, height =6 )
histo_NN <- plot_histogram(data = NN_NS_scores, x=NN_NS_scores$NN_score, col= "forestgreen",
                           title="Nature for nature score site's evaluation")
ggsave(filename = here::here("outputs", "figures","hist_NN_weighted_mean.png"), histo_NN, width = 8, height =6 )


mpa_NN <- ggplot(NN_NS_scores, aes(mpa_iucn_cat , NN_score , alpha = 0.5)) +
  geom_point()+
  geom_smooth() +
  labs(x = "IUCN categories", y = "NN score", title = "") +
  theme_minimal()+
  theme(legend.position = 'none')
ggsave(filename = here::here("outputs", "figures", "NN_scores according to MPA categories.png"),
       mpa_NN, width = 15, height = 12)


##NS
histo_NS <- plot_histogram(data = NN_NS_scores, x=NN_NS_scores$NS_score, col= "dodgerblue3",
                           title="Nature to People score site's evaluation")
ggsave(filename = here::here("outputs", "figures","hist_NS_weighted_mean.png"), histo_NS, width = 8, height =6 )


### plot scores according to socio-envir variables
plot_score_according_variables <- function(data = NN_NS_scores,
                                           score = "NN_score",
                                           var = "SiteLatitude"){
  ggplot(data, aes(x = data[,var][[var]],y= data[,score][[score]], alpha = 0.5)) +
    geom_point()+
    geom_smooth() +
    labs(x = var, y = score , title = "") +
    theme_minimal()+
    theme(legend.position = 'none')
}
variables <- c("SiteLatitude", "SiteLongitude", "SiteMeanSST", "Biomass",
               "HDI", "MarineEcosystemDependency", "coral_imputation",
               "gravtot2")

## NN
plots <- lapply( variables, function(var){
  plot_score_according_variables(data = NN_NS_scores, score = "NN_score", var = var)
  })

plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]]+
  plots[[6]] + plots[[7]] + plots[[8]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here::here("outputs", "figures","Socio_envir_variables_with_NN_scores.png"), plot, width = 22, height =14 )

## NS
plots <- lapply( variables, function(var){
  plot_score_according_variables(data = NN_NS_scores, score = "NS_score", var = var)
})

plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]]+
  plots[[6]] + plots[[7]] + plots[[8]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here::here("outputs", "figures","Socio_envir_variables_with_NS_scores.png"), plot, width = 22, height =14 )


## --------------- NN against NS -------------
#### Define colors by quarter ####
NN_NS_with_product <- NN_NS_scores |>
  dplyr::bind_cols( protection = ifelse(NN_NS_scores$mpa_enforcement == "High" &
                      stringr::str_detect(NN_NS_scores$protection_status, "No take"),
                          "No take",
                          ifelse(is.na(NN_NS_scores$mpa_name)==F, 
                            "Restricted", "Fished"))) |>
  # # dplyr::bind_cols(
  #   protection = ifelse(is.na(NN_NS_scores$mpa_name)==F &
  #                      NN_NS_scores$mpa_enforcement != "Low" &
  #                      stringr::str_detect(NN_NS_scores$protection_status, "No take"),
  #   # protection = ifelse(is.na(NN_NS_scores$mpa_iucn_cat)==F &
  #   #                     ( NN_NS_scores$mpa_iucn_cat == "Ia" |
  #   #                         NN_NS_scores$mpa_iucn_cat == "II") ,
  #                      "protected","not protected")) |>
  dplyr::mutate(NNxNS = abs(NN_score * NS_score)) |>
  dplyr::bind_cols(rank = rank(abs(NN_NS_scores$NN_score * NN_NS_scores$NS_score))) |>
  dplyr::bind_cols(up_right = ifelse(NN_NS_scores$NN_score > 0 & NN_NS_scores$NS_score > 0,1,NA)) |>
  dplyr::bind_cols(up_left = ifelse(NN_NS_scores$NN_score < 0 & NN_NS_scores$NS_score > 0,1,NA)) |>
  dplyr::bind_cols(down_right = ifelse(NN_NS_scores$NN_score > 0 & NN_NS_scores$NS_score < 0,1,NA)) |>
  dplyr::bind_cols(down_left = ifelse(NN_NS_scores$NN_score < 0 & NN_NS_scores$NS_score < 0,1,NA))

NN_NS_with_product <- NN_NS_with_product |>
  dplyr::bind_cols(rank_u_r = rank(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.last="keep")) |>
  dplyr::bind_cols(rank_u_l = rank(NN_NS_with_product$NNxNS * NN_NS_with_product$up_left, na.last="keep")) |>
  dplyr::bind_cols(rank_d_r = rank(NN_NS_with_product$NNxNS * NN_NS_with_product$down_right, na.last="keep")) |>
  dplyr::bind_cols(rank_d_l = rank(NN_NS_with_product$NNxNS * NN_NS_with_product$down_left, na.last="keep"))
   
summary(NN_NS_with_product)
save(NN_NS_with_product, file = here::here("outputs","NN_NS_score_wheighted_mean_quarter_ranked.Rdata"))
# median_curve_u_r <- data.frame(x_curve=
#                                  c(0,0.05,sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
#                                    na.rm=T)), 3* sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
#                                    na.rm=T)), 4* sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
#                                                              na.rm=T))), 
#                                y_curve=c(4* sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
#                                                         na.rm=T)), 
#                                          3* sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
#                                                         na.rm=T)), 
#                                          sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
#                                   na.rm=T)) , 0.01,0))



#### plot NN against NN ####
library(ggplot2)
NN_NS_plot <- ggplot(NN_NS_with_product, aes( y= NS_score, x = NN_score) ) +
  #up right quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, up_right == 1),
             size = 2, aes(colour= rank_u_r, shape = protection)) +
  scale_colour_gradient(name="up_right",
                        low = "white", high="firebrick",
                        limits =quantile(NN_NS_with_product$rank_u_r,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  geom_point(aes( y= NS_score, x = NN_score, shape = protection),
             size = 2, stroke = 1,
             color= "darkred",
             data = NN_NS_with_product[which(NN_NS_with_product$rank_u_r >=
                    quantile(NN_NS_with_product$rank_u_r, probs=c(0.95), na.rm=T)), ] )+
  ggrepel::geom_label_repel( aes(label=paste(SiteCode, ":", SiteEcoregion)), size=3, 
                             nudge_y = 0.1, nudge_x = 0.1,
             data=NN_NS_with_product[which(NN_NS_with_product$rank_u_r >=
                  quantile(NN_NS_with_product$rank_u_r, probs=c(0.95), na.rm=T)), ] )+
  
  guides(color = "none") + ggnewscale::new_scale("colour") +
  
  #up left quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, up_left == 1),
             size = 2, aes(colour= rank_u_l, shape = protection)) +
  scale_colour_gradient(name="up_left",
                        low = "white", high="dodgerblue3",
                        limits =quantile(NN_NS_with_product$rank_u_l,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  geom_point(aes( y= NS_score, x = NN_score, shape = protection),
             size = 2, stroke = 1,
             color= "dodgerblue4",
             data = NN_NS_with_product[which(NN_NS_with_product$rank_u_l >=
                        quantile(NN_NS_with_product$rank_u_l, probs=c(0.95), na.rm=T)), ] )+
  ggrepel::geom_label_repel( aes(label=paste(SiteCode, ":", SiteEcoregion)), size=3, 
                data=NN_NS_with_product[which(NN_NS_with_product$rank_u_l >=
                        quantile(NN_NS_with_product$rank_u_l, probs=c(0.95), na.rm=T)), ] )+
  guides(color = "none") + ggnewscale::new_scale("colour") +
  
  #down right quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, down_right == 1),
             size = 2, aes(colour= rank_d_r, shape = protection)) +
  scale_colour_gradient(name="down_right",
                        low = "white", high="forestgreen",
                        limits =quantile(NN_NS_with_product$rank_d_r,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  geom_point(aes( y= NS_score, x = NN_score, shape = protection),
             size = 2, stroke = 1,
             color= "darkgreen",
             data = NN_NS_with_product[which(NN_NS_with_product$rank_d_r >=
                       quantile(NN_NS_with_product$rank_d_r, probs=c(0.95), na.rm=T)), ] )+
  ggrepel::geom_label_repel( aes(label=paste(SiteCode, ":", SiteEcoregion)), size=3, 
             nudge_y = -0.2,
             data=NN_NS_with_product[which(NN_NS_with_product$rank_d_r >=
                       quantile(NN_NS_with_product$rank_d_r, probs=c(0.95), na.rm=T)), ] )+
  guides(color = "none") + ggnewscale::new_scale("colour") +
  
  #down left quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, down_left == 1),
             size = 2, aes(colour= rank_d_l, shape = protection)) +
  scale_colour_gradient(name="down_left",
                        low = "white", high="grey20",
                        limits =quantile(NN_NS_with_product$rank_d_l,
                                         probs=c( 0.5,1), na.rm=T), na.value=NA) +
  geom_point(aes( y= NS_score, x = NN_score, shape = protection),
             size = 2, stroke = 1,
             color= "black",
             data = NN_NS_with_product[which(NN_NS_with_product$rank_d_l >=
                       quantile(NN_NS_with_product$rank_d_l, probs=c(0.95), na.rm=T)), ] )+
  ggrepel::geom_label_repel( aes(label= paste(SiteCode, ":", SiteEcoregion)), size=3, 
             data=NN_NS_with_product[which(NN_NS_with_product$rank_d_l >=
                       quantile(NN_NS_with_product$rank_d_l, probs=c(0.95), na.rm=T)), ] )+
  
  # see MPAs
  scale_shape_manual(values=c(20,17,18))+
  #add lines
  geom_vline(xintercept = 0, linetype = 1, col = "black", linewidth = 0.2)+
  geom_hline(yintercept = 0, linetype = 1, col = "black", linewidth = 0.2)+
  
  #add 50% square:
  # geom_function(aes(x= NN_score, y = NS_score),
  #               fun = function(x){median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
  #                                        na.rm=T) / x},
  #               xlim=c(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
  #                             na.rm=T) / max(NN_NS_with_product$NS_score),
  #                      max(NN_NS_with_product$NN_score))) +
  # geom_polygon(data = data.frame(x=c(0,2*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
  #                                    na.rm=T)), 0, -2*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
  #                                    na.rm=T))),
  #                                y = c(2*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
  #                                    na.rm=T)),0, -2*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right,
  #                                    na.rm=T)),0)),
  #              aes(x, y), alpha= 0, color = "black",  linetype = 3,) +
  # geom_smooth(data = median_curve_u_r, aes(x_curve, y_curve), se=F, span = 1, linetype= 3,
  #             size=0.5, method = "loess", color="black" ) +
  
  geom_curve(aes(x=0, xend=2.5*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)),
                y=3*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)), yend=0),
                curvature=0.3, linetype=3, linewidth=0.1)+  #up_right
  geom_curve(aes(x=0, xend=2.5*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)),
                 y=-3*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)), yend=0),
             curvature=-0.3, linetype=3, linewidth=0.1)+   #down_right
  geom_curve(aes(x=0, xend=-2.5*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)),
                 y=3*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)), yend=0),
             curvature=-0.3, linetype=3, linewidth=0.1)+ #up_left
  geom_curve(aes(x=0, xend=-2.5*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)),
                 y=-3*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)), yend=0),
             curvature=0.3, linetype=3, linewidth=0.1)+ #down_right
  

  labs( x=  "Nature to Nature", y = "Nature to People")+
  theme_minimal()+
  theme(legend.position = c(0.9,0.1),
        legend.background = element_rect())+
  guides(color = "none", shape = guide_legend("Protection status")) 

#NN_NS_plot_legend <- ggfun::keybox( "roundrect", gp = ggfun::gpar(lty="dashed"))

NN_NS_plot
ggsave( here::here("outputs", "figures", "Sites in NN and NS scores.png"), 
        plot = NN_NS_plot, 
        width=15, height = 12 )


## with mpa proportion
plot_piechart <- function(quarter= "up_right", 
                          col="firebrick"){
  df <- as.data.frame(table(dplyr::filter(NN_NS_with_product, is.na(get(quarter))==F)$protection)) |>
    dplyr::mutate(csum = rev(cumsum(rev(Freq))),
                  pct = Freq / sum(Freq) *100,
                  pos = Freq/2 + dplyr::lead(csum, 1),
                  pos = dplyr::if_else(is.na(pos), Freq/2, pos))
  
  ggplot(df, aes(x = "", y = Freq,
                 fill = forcats::fct_inorder(Var1))) +
    # geom_col(width = 1, color = 1, position = "stack") +
    ggpattern::geom_col_pattern(pattern = c("none", "stripe","none"),
                                pattern_density = 0.001,
                                pattern_spacing = 0.05,
                                color = "grey20", position = "stack") +
    geom_text(aes(label = paste(round(pct,1), "%")),
              position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    guides(fill = guide_legend(title = "status")) +
    scale_y_continuous(breaks = df$pos, labels = df$Var1) +
    scale_fill_manual(values = c("grey", col, col))+
    theme_void()+
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          # axis.text = element_text(size = 10),
          legend.position = "none") # Removes the legend
}

plot_stack <- function(quarter= "up_right", 
                       col_restricted = "coral3",
                       col_notake = "darkred"){
  df <- as.data.frame(table(dplyr::filter(NN_NS_with_product, is.na(get(quarter))==F)$protection)) |>
    dplyr::mutate(pct = Freq / sum(Freq) *100) 
  
  ggplot(df, aes(x = "", y = pct, fill = Var1[c(3,1,2)])) +
    geom_col() +
    geom_text(aes(label = paste0(round(pct, 1), "%")),
              position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c(col_notake, col_restricted, "grey"))+
    theme_void() +  
    theme(legend.position = "none")
}

NN_NS_plot_stackchart <- NN_NS_plot + 
  annotation_custom(grob = ggplotGrob(plot_stack("up_right","coral3", "darkred")),
                    xmin = 1.3, xmax = 1.5, 
                    ymin = 0.95, ymax =1.55 )+ 
  annotation_custom(grob = ggplotGrob(plot_stack("down_right","darkseagreen3", "forestgreen")),
                    xmin = 0.7, xmax = 0.9, 
                    ymin = -2, ymax =-1.4 )+ 
  annotation_custom(grob = ggplotGrob(plot_stack("down_left","grey40", "black")),
                    xmin = -1.9, xmax = -1.7, 
                    ymin = -2, ymax =-1.4 )+ 
  annotation_custom(grob = ggplotGrob(plot_stack("up_left","deepskyblue3", "darkblue")),
                    xmin = -1.9, xmax = -1.7, 
                    ymin = 0.95, ymax =1.55 )
NN_NS_plot_stackchart
ggsave(here::here("outputs", "figures", "Sites in NN and NS scores _ with stackchart.png"),
        plot = NN_NS_plot_stackchart, width=10, height = 8 )


## on map 
function_NN_NS_on_map <- function(coord_NN_NS = NN_NS_with_product, ylim = c(-36, 31),
                                  xlim= c(-180,180), title=""){
  ggplot(coord_NN_NS) +
    geom_sf(data = coast, color = NA, fill = "lightgrey") +
    
    #down left quarter
    geom_point(data= dplyr::filter(coord_NN_NS,  down_left == 1),
               size = 2, alpha = 0.7,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="down_left",
                          low = "white", high="grey20",
                          limits =quantile(coord_NN_NS$rank * coord_NN_NS$down_left,
                                           probs=c( 0.5,1), na.rm=T), na.value=NA) +
    ggnewscale::new_scale("colour") +
    
    #up left quarter
    geom_point(data= dplyr::filter(coord_NN_NS,  up_left == 1),
               size = 2, alpha = 0.7,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="up_left",
                          low = "white", high="dodgerblue3",
                          limits =quantile(coord_NN_NS$rank * coord_NN_NS$up_left,
                                           probs=c( 0.5,1), na.rm=T), na.value=NA) +
    geom_point(aes( x= SiteLongitude, y = SiteLatitude),
               size = 2, stroke = 1, shape = 1,
               color= "darkblue",
               data = coord_NN_NS[which(coord_NN_NS$rank_u_l >=
                                          quantile(coord_NN_NS$rank_u_l, probs=c(0.99), na.rm=T)), ] )+
    
    ggnewscale::new_scale("colour") +
    
    #down right quarter
    geom_point(data= dplyr::filter(coord_NN_NS,  down_right == 1),
               size = 2, alpha = 0.7,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="down_right",
                          low = "white", high="forestgreen",
                          limits =quantile(coord_NN_NS$rank * coord_NN_NS$down_right,
                                           probs=c( 0.5,1), na.rm=T), na.value=NA) +
    geom_point(aes( x= SiteLongitude, y = SiteLatitude),
               size = 2, stroke = 1, shape = 1,
               color= "darkgreen",
               data = coord_NN_NS[which(coord_NN_NS$rank_d_r >=
                                          quantile(coord_NN_NS$rank_d_r, probs=c(0.99), na.rm=T)), ] )+
    
    ggnewscale::new_scale("colour") +
    
    #up right quarter
    geom_point(data= dplyr::filter(coord_NN_NS,  up_right == 1),
               size = 2, alpha = 0.7,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="up_right",
                          low = "white", high="darkgoldenrod3",
                          limits =quantile(coord_NN_NS$rank * coord_NN_NS$up_right,
                                           probs=c( 0.5,1), na.rm=T), na.value=NA) +
    geom_point(aes( x= SiteLongitude, y = SiteLatitude),
               size = 2, stroke = 1, shape = 1,
               color= "darkgoldenrod4",
               data = coord_NN_NS[which(coord_NN_NS$rank_u_r >=
                         quantile(coord_NN_NS$rank_u_r, probs=c(0.99), na.rm=T)), ] )+

    #     #add transparency
    # scale_alpha_continuous(range = c(-0.5, 0.6)) +
    
    
    coord_sf(xlim= xlim, ylim = ylim, expand = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_minimal()+
    labs(title = title,
         x="Longitude", y= "Latitude") +
    theme(legend.position = "none",
          plot.title = element_text(size=10, face="bold"),
          # axis.text.x = element_blank(),
          # axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
}

world_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_product, 
                                          ylim = c(-36, 31),
                                          xlim= c(-180,180),
                                          title = "Trade-offs in Nature Based Contribution Worldwide")
ggsave( here::here("outputs", "figures", "world map with NN and NS score.png"), plot = world_map_NN_NS, width=10, height = 6 )


australia_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_product, 
                                              ylim = c(-39, 0),
                                              xlim= c(100,180),
                                              title = "Trade-offs in Nature Based Contribution in Australia")

ggsave( here::here("outputs", "figures", "australia map with NN and NS score.png"), plot = australia_map_NN_NS, width=10, height = 6 )


polynesia_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_product, 
                                              ylim = c(-25, -10),
                                              xlim= c(-180,-130),
                                              title = "Trade-offs in Nature Based Contribution in polynesia")
ggsave( here::here("outputs", "figures", "polynesia map with NN and NS score.png"), plot = polynesia_map_NN_NS, width=10, height = 6 )


caraib_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_product, 
                                              ylim = c(-5, 20),
                                              xlim= c(-95,-70),
                                              title = "Trade-offs in Nature Based Contribution in caraïbe")
ggsave( here::here("outputs", "figures", "caraïbe map with NN and NS score.png"), plot = caraib_map_NN_NS, width=10, height = 6 )


#world map with zooms
world_map_zoom <- world_map_NN_NS + 
  geom_rect(aes(xmin = 140, xmax = 160, ymin = -25, ymax = -10), color = "black", fill= "transparent")+
  geom_rect(aes(xmin = -95, xmax = -70, ymin = -5, ymax = 20), color = "black", fill= "transparent")


gold_coast_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_product, 
                                              ylim = c(-25,-10),
                                              xlim= c(140,160),
                                              title = "")
caraib_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_product, 
                                           ylim = c(-5, 20),
                                           xlim= c(-95,-70),
                                           title = "")


ggpubr::ggarrange(world_map_zoom, # First row with world map
                  ggpubr::ggarrange(caraib_map_NN_NS, gold_coast_map_NN_NS, 
                                    ncol = 2, labels = c("B", "C")), # Second row with zooms
                  nrow = 2, labels = "A") 
ggsave(plot = last_plot(), width=11, height =7,
       filename = here::here("outputs", "figures", "world map with NN and NS score with zoom.png"))

##------------------------Study MPAs distribution ------------------------------------
mpa_distribution <- NN_NS_with_product |>
  dplyr::mutate(quarter = ifelse(!is.na(up_right), "up_right", ifelse(
                            !is.na(up_left), "up_left", ifelse(
                              !is.na(down_left), "down_left", "down_right")))) |>
  dplyr::select(SiteCode, protection, quarter)

contingency_table <-table(mpa_distribution$quarter, 
                                     mpa_distribution$protection)
margin.table(contingency_table,1)
margin.table(contingency_table,2)
prop.table(contingency_table,1)

#chi-squared test
summary(contingency_table)
# Test for independence of all factors:
#   Chisq = 19.067, df = 6, p-value = 0.004051

# Effect size
lsr::cramersV(contingency_table) #0.0877899

## which proportion of mpa in NN and NS quantiles:
df_plot_mpa <- NN_NS_with_product |>
  dplyr::filter(NN_score > quantile(NN_NS_with_product$NN_score, 0.75)) |>
  dplyr::mutate(score = "Nature to Nature", quantile = "best 25%") |>
  dplyr::bind_rows( dplyr::mutate(
    dplyr::filter(NN_NS_with_product,
                  NN_score < quantile(NN_NS_with_product$NN_score, 0.25)),
    score = "Nature to Nature", quantile = "worst 25%")) |>
  dplyr::bind_rows( dplyr::mutate(
    dplyr::filter(NN_NS_with_product,
                  NS_score > quantile(NN_NS_with_product$NS_score, 0.75)),
    score = "Nature to Society", quantile = "best 25%")) |>
  dplyr::bind_rows( dplyr::mutate(
    dplyr::filter(NN_NS_with_product,
                  NS_score < quantile(NN_NS_with_product$NS_score, 0.25)),
    score = "Nature to Society", quantile = "worst 25%")) |>
  dplyr::group_by(score, quantile) |>
  dplyr::summarize(No_take = sum(protection == "No take") / dplyr::n(),
                   Restricted = sum(protection == "Restricted") / dplyr::n()) |>
  tidyr::gather(`No_take`, `Restricted`, key = "protection", value = "prop")

df_plot_mpa <- dplyr::bind_cols(df_plot_mpa,
                                hatches = ifelse(df_plot_mpa$protection != "No_take",
                                                 "none", "stripe"))

  
prop_plot <- ggplot(df_plot_mpa) +
  aes(x = quantile, y = prop, pattern = hatches, fill = score) +
  ggpattern::geom_col_pattern( 
           pattern_density = 0.001,
           pattern_spacing = 0.05,
           color = "black", linewidth = 0.2)+
  geom_hline(yintercept =  sum(NN_NS_with_product$protection == "No take" |
                                 NN_NS_with_product$protection == "Restricted")/nrow(NN_NS_with_product),
             color = "black", linetype = "dashed")+
  scale_fill_manual(values = c("forestgreen", "dodgerblue3"))+
  ggpattern::scale_pattern_manual(values = c("none", "stripe"))+
  labs( x = "Externe quartiles in NN and NS scores",
        y = "Percentage of protected sites") +
  theme_light() +
  theme(legend.position = "none") +
  facet_wrap(vars(score))

prop_plot
ggsave( here::here("outputs", "figures", "Proportion of MPA in NN and NS quartiles.png"), 
        plot = prop_plot, width=6, height = 6 )


#### map of NN or NS only ####
map_NN_or_NS <- function(coord_NN_NS = NN_NS_with_product,
                         NCP = NN_NS_with_product$NN_score ,
                         col_NCP= "forestgreen",
                         ylim = c(-36, 31),
                         xlim= c(-180,180), title=""){
  ggplot(coord_NN_NS) +
  geom_sf(data = coast, color = NA, fill = "lightgrey") +
  geom_point(aes(x = SiteLongitude, y = SiteLatitude,
                 color = NCP, alpha= abs(NCP),
                 size = 2)) +
  geom_point(aes(x = SiteLongitude, y = SiteLatitude),
             shape = 1, size = 2, stroke = 0.5,
             color= "black",
             data = head(coord_NN_NS[order(NCP, decreasing = T),], 15)) +
  
  scale_colour_gradientn(colours = colorRampPalette(rev(c( col_NCP ,"white", "grey30")))(1000)) +
  scale_alpha_continuous(range = c(0, 1)) +
  
  coord_sf(xlim=xlim, ylim = ylim, expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_size_continuous(range = c(0.5, 4), guide = "none") +
  theme_minimal()+
  labs(title = title,
       x="Longitude", y= "Latitude") +
  theme(
    legend.position = "none",
    plot.title = element_text(colour = col_NCP, 
                              face = "bold", size = 12,
                              margin=margin(t = 20, b = -15),
                              hjust = 0.01))
}
## NN
NN_worldwide <- map_NN_or_NS(coord_NN_NS = NN_NS_with_product,
                            NCP = NN_NS_with_product$NN_score,
                            col_NCP= "forestgreen",
                            ylim = c(-36, 31),
                            xlim= c(-180,180), title="Nature to Nature")
ggsave( here::here("outputs", "figures", "world map with NN score.png"), plot = NN_worldwide, width=10, height = 6 )



## NS
NS_worldwide <- map_NN_or_NS(coord_NN_NS = NN_NS_with_product,
                             NCP = NN_NS_with_product$NS_score,
                             col_NCP= "dodgerblue3",
                             ylim = c(-36, 31),
                             xlim= c(-180,180), title="Nature to People")
ggsave( here::here("outputs", "figures", "world map with NS score.png"), plot = NS_worldwide, width=10, height = 6 )

#panel
ggpubr::ggarrange(NN_worldwide, NS_worldwide,
                  nrow = 2, labels = c("A", "B")) 
ggsave(plot = last_plot(), width=11, height =6,
       filename = here::here("outputs", "figures", "world map with NN and NS score SEPARATLY.png"))


#### plot outliers only ####
function_outliers_on_map <- function(coord_NN_NS = NN_NS_with_product, ylim = c(-36, 31),
                                  xlim= c(-180,180), title=""){
  ggplot(coord_NN_NS) +
    geom_sf(data = coast, color = NA, fill = "lightgrey") +
    
    #down left quarter
    geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
               size = 4, stroke = 1, 
               color= "grey10",
               data = coord_NN_NS[which(coord_NN_NS$rank_d_l >=
                                          quantile(coord_NN_NS$rank_d_l, probs=c(0.95), na.rm=T)), ] )+
    
    ggnewscale::new_scale("colour") +
    
    #up left quarter
    geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
               size = 4, stroke = 1,
               color= "dodgerblue3",
               data = coord_NN_NS[which(coord_NN_NS$rank_u_l >=
                                          quantile(coord_NN_NS$rank_u_l, probs=c(0.95), na.rm=T)), ] )+
    
    ggnewscale::new_scale("colour") +
    
    #down right quarter
    geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
               size = 4, stroke = 1, 
               color= "forestgreen",
               data = coord_NN_NS[which(coord_NN_NS$rank_d_r >=
                                          quantile(coord_NN_NS$rank_d_r, probs=c(0.95), na.rm=T)), ] )+
    
    ggnewscale::new_scale("colour") +
    
    #up right quarter
    geom_point(aes( x= SiteLongitude, y = SiteLatitude, shape = protection),
               size = 4, stroke = 1, 
               color= "darkgoldenrod3",
               data = coord_NN_NS[which(coord_NN_NS$rank_u_r >=
                                          quantile(coord_NN_NS$rank_u_r, probs=c(0.95), na.rm=T)), ] )+
    
    #add shape
    scale_shape_manual(values=c(20,17,18))+
    
    coord_sf(xlim= xlim, ylim = ylim, expand = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    # scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_minimal()+
    labs(title = title,
         x="Longitude", y= "Latitude") +
    theme(legend.position = "right",
          plot.title = element_text(size=10, face="bold"),
          # axis.text.x = element_blank(),
          # axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
}

#world
outliers_world <- function_outliers_on_map( NN_NS_with_product, 
                       ylim = c(-36, 31),
                       xlim= c(-180,180),
                       title = "Trade-offs in Nature Based Contribution Worldwide")
ggsave( here::here("outputs", "figures", "world map with NN and NS outliers.png"), 
        plot = outliers_world, width=11, height = 6 )

#Australia
function_outliers_on_map( NN_NS_with_product, 
                       ylim = c(-39, 0),
                       xlim= c(100,180),
                       title = "Trade-offs in Nature Based Contribution in Australia")

#Caraib
function_outliers_on_map( NN_NS_with_product, 
                       ylim = c(-5, 20),
                       xlim= c(-95,-70),
                       title = "Trade-offs in Nature Based Contribution in caraïbe")
#Polynesia
function_outliers_on_map( NN_NS_with_product, 
                       ylim = c(-25, -10),
                       xlim= c(-180,-130),
                       title = "Trade-offs in Nature Based Contribution in polynesia")
