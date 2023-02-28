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
          "ggfun", "scatterpie", "scales")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data------------
load(here::here("outputs","all_NCP_site.Rdata"))
# load(here::here("outputs","NCP_site_coral_reef.Rdata"))
# load(here::here("outputs","NCP_site_SST20.Rdata"))
# load(here::here("outputs","NCP_site_coral_5_imputed.Rdata"))
# load(here::here("outputs","NCP_site_wo_australia.Rdata"))
# NCP_site <- NCP_site_condition


load(here::here("data","metadata_surveys.Rdata"))
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

##-------------compute correlation coefficients-------------

## Clean data
NCP_site_clean <- subset(NCP_site, select = -c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                               SiteMeanSST, SiteLatitude, SiteLongitude,
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

plots <- lapply(colnames(NCP_site_clean), FUN = plot_distribution, data = NCP_site_clean )

library(patchwork)
all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] +  plots[[30]] + plots[[31]] + 
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

# ggsave(filename = here("outputs", "figures","NCP_distribution.png"), all_plot, width = 22, height =14 )

#### NCPs distribution with log correction
NCP_skewed_distribution <- c("Btot","recycling_N","recycling_P","Productivity",
                             "funct_distinctiveness","Omega_3_C","Calcium_C","Vitamin_A_C",
                             "phylo_entropy","ED_Mean", "iucn_species", "elasmobranch_diversity",
                             "low_mg_calcite", "high_mg_calcite", "aragonite", "monohydrocalcite",
                             "amorphous_carbonate", "biom_lowTL", "biom_mediumTL", "biom_highTL",
                             "fishery_biomass")

NCP_log_clean <- NCP_site_clean |>
  dplyr::mutate(across(.cols = all_of(NCP_skewed_distribution),
                       .fns = ~ .x +1 , .names = "{.col}")) |>      # Adds 1 to values to log transformed
  dplyr::mutate(across(.cols = all_of(NCP_skewed_distribution),
                       .fns = log10 , .names = "{.col}"))          # log(x+1) to avoid negative values
  
NCP_log_transformed <-NCP_log_clean |> dplyr::rename_with(.cols = all_of(NCP_skewed_distribution),
                                                          .fn = ~ paste0("log(", .x, ")"))

plots <- lapply(colnames(NCP_log_transformed), FUN = plot_distribution, data = NCP_log_transformed )
library(patchwork)
all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + plots[[30]] + plots[[31]] + 
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

grp <- as.factor( c(recycling_N="NN", recycling_P="NN",Productivity="NS",taxo_richness="NN", funct_entropy="NN",
                    funct_distinctiveness="NN", Selenium_C="NS", Zinc_C="NS", Omega_3_C="NS", Calcium_C="NS",
                    Iron_C="NS", Vitamin_A_C="NS", phylo_entropy="NN", ED_Mean="NN", aesthe_survey="NS", iucn_species="NN",
                    elasmobranch_diversity="NN", low_mg_calcite="NN", high_mg_calcite="NN", aragonite="NN",
                    monohydrocalcite="NN", amorphous_carbonate="NN", biom_lowTL="NN", biom_mediumTL="NN",
                    biom_highTL="NN", fishery_biomass="NS", mean_TL = "NN", robustness = "NN",
                    scientific_interest = "NS", public_interest = "NS")) # /!\ the order matter


corr_pearson_NCPs <- stats::cor(NCP_log_clean, method="pearson")
NCP_log_scale_clean <- scale(NCP_log_clean)

#### Nature to Nature (NN) score ####
NN_names <- names(grp)[ grp=="NN" ]
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
NS_names <- names(grp)[ grp=="NS" ]
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




NN_NS_scores <- cbind(NCP_site[,c("SiteCode", "SiteLongitude", "SiteLatitude", 
                                  "mpa_name", "mpa_enforcement", "protection_status",
                                  "mpa_iucn_cat")],
                      data.frame(NN_score = EDS_NN, NS_score = EDS_NS) )
save(NN_NS_scores, file = here::here("outputs", "NN_NS_score_wheighted_mean.Rdata"))

##-------------plot NN and NS scores-------------
library(ggplot2)

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


## --------------- NN against NS -------------
#### Define colors by quarter ####
NN_NS_with_product <- NN_NS_scores |>
  dplyr::bind_cols(
    protection = ifelse(is.na(NN_NS_scores$mpa_name)==F &
                       NN_NS_scores$mpa_enforcement != "Low" &
                       stringr::str_detect(NN_NS_scores$protection_status, "No take"),
    # protection = ifelse(is.na(NN_NS_scores$mpa_iucn_cat)==F &
    #                     ( NN_NS_scores$mpa_iucn_cat == "Ia" |
    #                         NN_NS_scores$mpa_iucn_cat == "II") ,
                       "protected","not protected")) |>
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
  ggrepel::geom_label_repel( aes(label=SiteCode), size=2, 
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
  ggrepel::geom_label_repel( aes(label=SiteCode), size=2, 
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
  ggrepel::geom_label_repel( aes(label=SiteCode), size=2, 
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
  ggrepel::geom_label_repel( aes(label=SiteCode), size=2, 
                             data=NN_NS_with_product[which(NN_NS_with_product$rank_d_l >=
                                                             quantile(NN_NS_with_product$rank_d_l, probs=c(0.95), na.rm=T)), ] )+
  
  # see MPAs
  scale_shape_manual(values=c(20,17))+
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
                curvature=0.4, linetype=3, size=0.1)+  #up_right
  geom_curve(aes(x=0, xend=2.5*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)),
                 y=-3*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)), yend=0),
             curvature=-0.4, linetype=3, size=0.1)+   #down_right
  geom_curve(aes(x=0, xend=-2.5*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)),
                 y=3*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)), yend=0),
             curvature=-0.4, linetype=3, size=0.1)+ #up_left
  geom_curve(aes(x=0, xend=-2.5*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)),
                 y=-3*sqrt(median(NN_NS_with_product$NNxNS * NN_NS_with_product$up_right, na.rm=T)), yend=0),
             curvature=0.3, linetype=3, size=0.1)+ #down_right
  

  labs( x=  "Nature to Nature", y = "Nature to People")+
  theme_minimal()+
  theme(legend.position = c(0.9,0.1),
        legend.background = element_rect())+
  guides(color = "none", shape = guide_legend("Protection status")) 

#NN_NS_plot_legend <- ggfun::keybox( "roundrect", gp = ggfun::gpar(lty="dashed"))

NN_NS_plot
ggsave( here::here("outputs", "figures", "Sites in NN and NS scores.png"), plot = NN_NS_plot, width=10, height = 8 )


## with piechart
plot_piechart <- function(quarter= "up_right", col="firebrick"){
  df <- as.data.frame(table(dplyr::filter(NN_NS_with_product, is.na(get(quarter))==F)$protection)) |>
    dplyr::mutate(csum = rev(cumsum(rev(Freq))),
                  pct = Freq / sum(Freq) *100,
                  pos = Freq/2 + dplyr::lead(csum, 1),
                  pos = dplyr::if_else(is.na(pos), Freq/2, pos))
  
  ggplot(df, aes(x = "", y = Freq, 
                 fill = forcats::fct_inorder(Var1))) +
    geom_col(width = 1, color = 1, position = "stack") +
    geom_text(aes(label = paste(round(pct,1), "%")),
              position = position_stack(vjust = 0.5)) +
    coord_polar(theta = "y") +
    guides(fill = guide_legend(title = "status")) +
    scale_y_continuous(breaks = df$pos, labels = df$Var1) +
    scale_fill_manual(values = c("grey", col))+
    theme_void()+
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          # axis.text = element_text(size = 10),
          legend.position = "none") # Removes the legend
}

NN_NS_plot_piechart <- NN_NS_plot + 
  annotation_custom(grob = ggplotGrob(plot_piechart("up_right","firebrick")),
                    xmin = 1.1, xmax = 1.7, 
                    ymin = 0.85, ymax =1.45 )+ 
  annotation_custom(grob = ggplotGrob(plot_piechart("down_right","forestgreen")),
                    xmin = 0.5, xmax = 1.1, 
                    ymin = -2, ymax =-1.4 )+ 
  annotation_custom(grob = ggplotGrob(plot_piechart("down_left","grey20")),
                    xmin = -2, xmax = -1.4, 
                    ymin = -2, ymax =-1.4 )+ 
  annotation_custom(grob = ggplotGrob(plot_piechart("up_left","dodgerblue3")),
                    xmin = -2, xmax = -1.4, 
                    ymin = 0.85, ymax =1.45 )
NN_NS_plot_piechart
ggsave(here::here("outputs", "figures", "Sites in NN and NS scores _ with piecharts.png"),
        plot = NN_NS_plot_piechart, width=10, height = 8 )

## on map 
function_NN_NS_on_map <- function(coord_NN_NS = NN_NS_with_product, ylim = c(-36, 31),
                                  xlim= c(-180,180), title=""){
  ggplot(coord_NN_NS) +
    geom_sf(data = coast, color = NA, fill = "lightgrey") +
    
    #down left quarter
    geom_point(data= dplyr::filter(coord_NN_NS,  down_left == 1),
               size = 2, alpha = 0.5,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="down_left",
                          low = "white", high="grey20",
                          limits =quantile(NN_NS_with_product$rank * NN_NS_with_product$down_left,
                                           probs=c( 0.5,1), na.rm=T), na.value=NA) +
    ggnewscale::new_scale("colour") +
    
    #up left quarter
    geom_point(data= dplyr::filter(coord_NN_NS,  up_left == 1),
               size = 2, alpha = 0.5,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="up_left",
                          low = "white", high="dodgerblue3",
                          limits =quantile(NN_NS_with_product$rank * NN_NS_with_product$up_left,
                                           probs=c( 0.5,1), na.rm=T), na.value=NA) +
    # geom_point(aes(x = SiteLongitude, y = SiteLatitude),
    #            shape = 1, size = 2, stroke = 1,
    #            color= "dodgerblue4",
    #            data = head(coord_NN_NS[order(coord_NN_NS$product_u_l, decreasing = T),], 15)) +
    
    ggnewscale::new_scale("colour") +
    
    #down right quarter
    geom_point(data= dplyr::filter(coord_NN_NS,  down_right == 1),
               size = 2, alpha = 0.5,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="down_right",
                          low = "white", high="forestgreen",
                          limits =quantile(NN_NS_with_product$rank * NN_NS_with_product$down_right,
                                           probs=c( 0.5,1), na.rm=T), na.value=NA) +
    # geom_point(aes(x = SiteLongitude, y = SiteLatitude),
    #            shape = 1, size = 2, stroke = 1,
    #            color= "darkgreen",
    #            data = head(coord_NN_NS[order(coord_NN_NS$product_d_r, decreasing = T),], 15)) +
    
    ggnewscale::new_scale("colour") +
    
    #up right quarter
    geom_point(data= dplyr::filter(coord_NN_NS,  up_right == 1),
               size = 2, alpha = 0.5,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="up_right",
                          low = "white", high="darkred",
                          limits =quantile(NN_NS_with_product$rank * NN_NS_with_product$up_right,
                                           probs=c( 0.5,1), na.rm=T), na.value=NA) +
    
    # geom_point(aes(x = SiteLongitude, y = SiteLatitude),
    #            shape = 1, size = 2, stroke = 1,
    #            color= "black",
    #            data = head(coord_NN_NS[order(coord_NN_NS$product_u_r, decreasing = T),], 15)) +
    
    #add transparency
    scale_alpha_continuous(range = c(-0.5, 0.6)) +
    
    
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
                                              ylim = c(-5, 30),
                                              xlim= c(-100,-55),
                                              title = "Trade-offs in Nature Based Contribution in caraïbe")
ggsave( here::here("outputs", "figures", "caraïbe map with NN and NS score.png"), plot = caraib_map_NN_NS, width=10, height = 6 )


#### which proportion of mpa in quantiles: ####
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
  dplyr::summarize(prop_protect = sum(protection == "protected") / dplyr::n())

  
prop_plot <- ggplot(df_plot_mpa) +
  aes(x = quantile, y = prop_protect, fill = score) +
  geom_col() +
  geom_hline(yintercept =  sum(NN_NS_with_product$protection == "protected")/nrow(NN_NS_with_product),
             color = "black", linetype = "dashed")+
  scale_fill_manual(values = c("forestgreen", "dodgerblue3"))+
  labs( x = "Externe quartiles in NN and NS scores",
        y = "Percentage of protected sites") +
  theme_light() +
  theme(legend.position = "none")+
  facet_wrap(vars(score))
prop_plot
ggsave( here::here("outputs", "figures", "Proportion of MPA in NN and NS quartiles.png"), 
        plot = prop_plot, width=6, height = 6 )


#### NN or NS only ####
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
    plot.title = element_text(size=10, face="bold"))
}
## NN
NN_worldwide <- map_NN_or_NS(coord_NN_NS = NN_NS_with_product,
                            NCP = NN_NS_with_product$NN_score,
                            col_NCP= "forestgreen",
                            ylim = c(-36, 31),
                            xlim= c(-180,180), title="Nature to nature contributions worldwide")
ggsave( here::here("outputs", "figures", "world map with NN score.png"), plot = NN_worldwide, width=10, height = 6 )



## NS
NS_worldwide <- map_NN_or_NS(coord_NN_NS = NN_NS_with_product,
                             NCP = NN_NS_with_product$NS_score,
                             col_NCP= "dodgerblue3",
                             ylim = c(-36, 31),
                             xlim= c(-180,180), title="Nature to People contributions worldwide")
ggsave( here::here("outputs", "figures", "world map with NS score.png"), plot = NS_worldwide, width=10, height = 6 )

## According to latitude
latitude_NN <- ggplot(NN_NS_scores, aes(SiteLatitude, NN_score , alpha = 0.5)) +
  geom_point()+
  geom_smooth() +
  labs(x = "Latitude", y = "NN score", title = "Importance of lattitude in NN score") +
  theme_minimal()+
  theme(legend.position = 'none')
ggsave(filename = here::here("outputs", "figures", "NN score according to latitude.png"), latitude_NN, width = 15, height = 12)

latitude_NS <- ggplot(NN_NS_scores, aes(SiteLatitude, NS_score , alpha = 0.5)) +
  geom_point()+
  geom_smooth() +
  labs(x = "Latitude", y = "NN score", title = "Importance of lattitude in NS score") +
  theme_minimal()+
  theme(legend.position = 'none')

ggsave(filename = here::here("outputs", "figures", "NS score according to latitude.png"), latitude_NS, width = 15, height = 12)
