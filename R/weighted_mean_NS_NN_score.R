################################################################################
##
## Calculates NN and NS scores of sites
##
## weighted_mean_NS_NN_score.R
##
## 17/01/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse", "ggplot2", "sf", "patchwork", "stats")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load(here::here("outputs","all_NCP_site.Rdata"))
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

##-------------compute correlation coefficients-------------

## Clean data
NCP_site_clean <- subset(NCP_site, select = -c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                               SiteMeanSST, SiteLatitude, SiteLongitude,
                                               HDI, gravtot2, MarineEcosystemDependency,
                                               coral_imputation))


#### NCPs distribution
plot_distribution <- function(ncp, data){
  col <- fishualize::fish(n = 27, option = "Ostracion_whitleyi", begin = 0, end = 0.8)
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



##-------------compute correlation weighted mean of NCP (Kark 2002): NN and NS scores-------------

grp <- as.factor( c(recycling_N="NN", recycling_P="NN",Productivity="NS",taxo_richness="NN", funct_entropy="NN",
                    funct_distinctiveness="NN", Selenium_C="NS", Zinc_C="NS", Omega_3_C="NS", Calcium_C="NS",
                    Iron_C="NS", Vitamin_A_C="NS", phylo_entropy="NN", ED_Mean="NN", aesthe_survey="NS", iucn_species="NN",
                    elasmobranch_diversity="NN", low_mg_calcite="NN", high_mg_calcite="NN", aragonite="NN",
                    monohydrocalcite="NN", amorphous_carbonate="NN", biom_lowTL="NN", biom_mediumTL="NN",
                    biom_highTL="NN", fishery_biomass="NS")) 

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




NN_NS_scores <- cbind(NCP_site[,c("SiteCode", "SiteLongitude", "SiteLatitude")],
                      data.frame(NN_score = EDS_NN, NS_score = EDS_NS) )
save(NN_NS_scores, file = here::here("outputs", "NN_NS_score_wheighted_mean.Rdata"))

##-------------plot NN and NS scores-------------
library(ggplot2)
## NN
histo_NN <- ggplot(NN_NS_scores, aes(x = NN_score)) +
  geom_histogram( bins = 40, aes(fill=after_stat(x)), color="grey50", linewidth=0.3)+
  scale_fill_gradient2(low = "black", mid = "white", high = "forestgreen", guide = "colourbar")+
  labs(x = "Nature for nature score site's evaluation")+
  geom_vline(xintercept = 0, linetype = 2, col= "black")+
  theme_minimal()
histo_NN
ggsave(filename = here("outputs", "figures","hist_NN_weighted_mean.png"), histo_NN, width = 8, height =6 )


##NS
histo_NS <- ggplot(NN_NS_scores, aes(x = NS_score)) +
  geom_histogram( bins = 40, aes(fill=after_stat(x)), color="grey30", linewidth=0.3)+
  scale_fill_gradient2(low = "black", mid = "white", high = "dodgerblue3", guide = "colourbar")+
  labs(x = "Nature to People score site's evaluation")+
  geom_vline(xintercept = 0, linetype = 2, col= "black")+
  theme_minimal()
histo_NS
ggsave(filename = here("outputs", "figures","hist_NS_weighted_mean.png"), histo_NS, width = 8, height =6 )


## --------------- NN against NS -------------
#### Define colors by quarter ####
library(dplyr)
#up right
quarter_up_right <- NN_NS_scores |>
  filter(NN_score > 0 & NS_score > 0) |>
  mutate(product_u_r =  (NN_score * NS_score)) |>
  arrange(product_u_r) |>
  tibble::rownames_to_column(var="scale_u_r")

#up left
quarter_up_left <- NN_NS_scores |>
  filter(NN_score < 0 & NS_score > 0) |>
  mutate(product_u_l =  (-NN_score * NS_score) ) |>
  arrange(product_u_l) |>
  tibble::rownames_to_column(var="scale_u_l")

#down right
quarter_down_right <- NN_NS_scores |>
  filter(NN_score > 0 & NS_score < 0) |>
  mutate(product_d_r =  (NN_score * -NS_score)) |>
  arrange(product_d_r) |>
  tibble::rownames_to_column(var="scale_d_r")

#down left
quarter_down_left <- NN_NS_scores |>
  filter(NN_score < 0 & NS_score < 0) |>
  mutate(product_d_l =  (NN_score * NS_score)  ) |>
  arrange(product_d_l) |>
  tibble::rownames_to_column(var="scale_d_l")

#Merge the 4 quarters
NN_NS_with_product <- quarter_down_left |>
  full_join(quarter_down_right) |>
  full_join(quarter_up_left) |>
  full_join(quarter_up_right)

#### plot NN against NN ####
NN_NS_plot <- ggplot(NN_NS_with_product, aes( y= NS_score, x = NN_score) ) +
  geom_point(data= dplyr::filter(NN_NS_with_product, is.na(scale_u_r)==F),
             size = 2,
             aes(colour= as.integer(scale_u_r), alpha = as.integer(scale_u_r))) +
  scale_colour_gradient(name="scale_u_r",
                        low = "white", high="darkred",
                        na.value=NA) +
  geom_point(aes( y= NS_score, x = NN_score),
             shape = 1, size = 2, stroke = 1,
             color= "black",
             data = head(NN_NS_with_product[order(NN_NS_with_product$product_u_r, decreasing = T),], 10)) +
  
  ggnewscale::new_scale("colour") +
  
  #up left quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, is.na(product_u_l)==F),
             size = 2,
             aes(colour= as.integer(scale_u_l), alpha = as.integer(scale_u_l))) +
  scale_colour_gradient(name="scale_u_l",
                        low = "white", high="dodgerblue3",
                        na.value=NA) +
  geom_point(aes( y= NS_score, x = NN_score),
             shape = 1, size = 2, stroke = 1,
             color= "dodgerblue4",
             data = head(NN_NS_with_product[order(NN_NS_with_product$product_u_l, decreasing = T),], 10)) +
  
  
  ggnewscale::new_scale("colour") +
  
  #down right quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, is.na(product_d_r)==F),
             size = 2,
             aes(colour= as.integer(scale_d_r), alpha = as.integer(scale_d_r))) +
  scale_colour_gradient(name="scale_d_r",
                        low = "white", high="forestgreen",
                        na.value=NA) +
  geom_point(aes( y= NS_score, x = NN_score),
             shape = 1, size = 2, stroke = 1,
             color= "darkgreen",
             data = head(NN_NS_with_product[order(NN_NS_with_product$product_d_r, decreasing = T),], 10)) +
  
  ggnewscale::new_scale("colour") +
  
  #down left quarter
  geom_point(data= dplyr::filter(NN_NS_with_product, is.na(product_d_l)==F),
             size = 2,
             aes(colour= as.integer(scale_d_l), alpha = as.integer(scale_d_l))) +
  scale_colour_gradient(name="scale_d_l",
                        low = "white", high="grey20",
                        na.value=NA) +
  
  scale_alpha_continuous(range = c(-0.5, 0.8)) +
  
  #add lines
  geom_vline(xintercept = 0, linetype = 2, col = "black")+
  geom_hline(yintercept = 0, linetype = 2, col = "black")+
  labs( x=  "Nature to Nature", y = "Nature to People")+
  theme_minimal()+
  theme(legend.position = "none")

NN_NS_plot
ggsave( here::here("outputs", "figures", "Sites in NN and NS scores.jpg"), plot = NN_NS_plot, width=10, height = 8 )



## on map 
function_NN_NS_on_map <- function(coord_NN_NS = NN_NS_with_product, ylim = c(-36, 31),
                                  xlim= c(-180,180), title=""){
  ggplot(coord_NN_NS) +
    geom_sf(data = coast, color = NA, fill = "lightgrey") +
    
    #down left quarter
    geom_point(data= dplyr::filter(coord_NN_NS, is.na(product_d_l)==F),
               size = 2, alpha = 0.5,
               aes(x = SiteLongitude, y = SiteLatitude,
                   colour= as.integer(scale_d_l), alpha = as.integer(scale_d_l))) +
    scale_colour_gradient(name="scale_d_l",
                          low = "white", high="grey20",
                          na.value=NA) +
    ggnewscale::new_scale("colour") +
    
    #up left quarter
    geom_point(data= dplyr::filter(coord_NN_NS, is.na(product_u_l)==F),
               size = 2, alpha = 0.5,
               aes(x = SiteLongitude, y = SiteLatitude,
                   colour= as.integer(scale_u_l), alpha = as.integer(scale_u_l))) +
    scale_colour_gradient(name="scale_u_l",
                          low = "white", high="dodgerblue3",
                          na.value=NA) +
    geom_point(aes(x = SiteLongitude, y = SiteLatitude),
               shape = 1, size = 2, stroke = 1,
               color= "dodgerblue4",
               data = head(coord_NN_NS[order(coord_NN_NS$product_u_l, decreasing = T),], 15)) +
    
    ggnewscale::new_scale("colour") +
    
    #down right quarter
    geom_point(data= dplyr::filter(coord_NN_NS, is.na(product_d_r)==F),
               size = 2, alpha = 0.5,
               aes(x = SiteLongitude, y = SiteLatitude,
                   colour= as.integer(scale_d_r), alpha = as.integer(scale_d_r))) +
    scale_colour_gradient(name="scale_d_r",
                          low = "white", high="forestgreen",
                          na.value=NA) +
    geom_point(aes(x = SiteLongitude, y = SiteLatitude),
               shape = 1, size = 2, stroke = 1,
               color= "darkgreen",
               data = head(coord_NN_NS[order(coord_NN_NS$product_d_r, decreasing = T),], 15)) +
    
    ggnewscale::new_scale("colour") +
    
    #up right quarter
    geom_point(data= dplyr::filter(coord_NN_NS, is.na(product_u_r)==F),
               size = 2,
               aes(x = SiteLongitude, y = SiteLatitude,
                   colour= as.integer(scale_u_r), alpha = as.integer(scale_u_r))) +
    scale_colour_gradient(name="scale_u_r",
                          low = "white", high="darkred",
                          na.value=NA) +
    geom_point(aes(x = SiteLongitude, y = SiteLatitude),
               shape = 1, size = 2, stroke = 1,
               color= "black",
               data = head(coord_NN_NS[order(coord_NN_NS$product_u_r, decreasing = T),], 15)) +
    
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
ggsave( here::here("outputs", "figures", "world map with NN and NS score.jpg"), plot = world_map_NN_NS, width=10, height = 6 )


australia_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_product, 
                                              ylim = c(-39, 0),
                                              xlim= c(100,180),
                                              title = "Trade-offs in Nature Based Contribution in Australia")

ggsave( here::here("outputs", "figures", "australia map with NN and NS score.jpg"), plot = australia_map_NN_NS, width=10, height = 6 )


polynesia_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_product, 
                                              ylim = c(-25, -10),
                                              xlim= c(-180,-130),
                                              title = "Trade-offs in Nature Based Contribution in polynesia")
ggsave( here::here("outputs", "figures", "polynesia map with NN and NS score.jpg"), plot = polynesia_map_NN_NS, width=10, height = 6 )


caraib_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_product, 
                                              ylim = c(-5, 30),
                                              xlim= c(-100,-55),
                                              title = "Trade-offs in Nature Based Contribution in caraïbe")
ggsave( here::here("outputs", "figures", "caraïbe map with NN and NS score.jpg"), plot = caraib_map_NN_NS, width=10, height = 6 )


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
ggsave( here::here("outputs", "figures", "world map with NN score.jpg"), plot = NN_worldwide, width=10, height = 6 )



## NS
NS_worldwide <- map_NN_or_NS(coord_NN_NS = NN_NS_with_product,
                             NCP = NN_NS_with_product$NS_score,
                             col_NCP= "dodgerblue3",
                             ylim = c(-36, 31),
                             xlim= c(-180,180), title="Nature to People contributions worldwide")
ggsave( here::here("outputs", "figures", "world map with NS score.jpg"), plot = NS_worldwide, width=10, height = 6 )

