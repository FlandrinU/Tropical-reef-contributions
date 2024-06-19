################################################################################
##
## This script reduces the 29 contributions into 2 synthetic dimensions, by 
##  meaning the contributions according to their category: either 
##  Nature-for-Nature (NN), either Nature-for-People (NP).
## This script gathers all variables needed for this study and save the final 
##  table in "outputs/tropical_reef_contributions_final_table.csv", and studies 
##  the correlation between both scores, spatial distributions, and their 
##  relationship to protection status.
##
## 1d_weighted_mean_NP_NN_score.R
##
## 17/01/2023
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
# load(here::here("outputs","all_Contrib_survey_log_transformed.Rdata"))
# load(here::here("outputs","Contrib_site_log_coral_reef.Rdata"))
# load(here::here("outputs","Contrib_site_log_SST20.Rdata"))
# load(here::here("outputs","Contrib_site_log_coral_5_imputed.Rdata"))
# load(here::here("outputs","Contrib_site_log_wo_australia.Rdata"))
# load(here::here("outputs","Contrib_site_log_random.Rdata"))
# load(here::here("outputs","Contrib_site_log_only_australia.Rdata"))
# Contrib_site_log_transformed <- Contrib_site_condition

#Add a year column
Contrib_site_log_transformed <- Contrib_site_log_transformed |> 
  dplyr::bind_cols(
    year = as.numeric(format(as.Date(Contrib_site_log_transformed$SurveyDate, 
                                     format = "%d/%m/%Y"),
                             format = "%Y")))


load(here::here("data","metadata_surveys.Rdata"))
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

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
                         Public_Attention = "NP")) # /!\ the order matter

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


#-------------Compute weighted mean of Contributions: NN and NS scores-------------
Contrib_log_scale_clean <- scale(Contrib_log_transformed)

#### Nature to Nature (NN) score ####
NN_names <- names(grp_NN_NP)[ grp_NN_NP=="NN" ]
Contrib_NN <- Contrib_log_scale_clean[,NN_names]

colnames(Contrib_NN)
weighting_par <- c(1/7, 1/7, 1/5, 1/5,1/5, 1/5,1/5,1/5,1/5,1/5,1/5,1/5,
                   1/7,1/7,1/7,1/7,1/7,1/2,1/2)
names(weighting_par) <- colnames(Contrib_NN)
weighting_par

mean_NN <- c()
for( site in 1:nrow(Contrib_NN)){
  EDS_site <- sum(weighting_par * Contrib_NN[site,]) / sum(weighting_par)
  mean_NN <- c(mean_NN, EDS_site) }

summary(mean_NN)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -1.8112620 -0.3350228  0.0007399  0.0000000  0.3480500  1.3931706


#### Nature to People (NP) score ####
NP_names <- names(grp_NN_NP)[ grp_NN_NP=="NP" ]
Contrib_NP <- Contrib_log_scale_clean[,NP_names]

colnames(Contrib_NP)
weighting_par <- c(1/2, 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/2, 1/2,1/2)
names(weighting_par) <- colnames(Contrib_NP)
weighting_par

mean_NP <- c()
for( site in 1:nrow(Contrib_NP)){
  EDS_site <- sum(weighting_par * Contrib_NP[site,]) / sum(weighting_par)
  mean_NP <- c(mean_NP, EDS_site) }

summary(mean_NP)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.96304 -0.23629  0.05656  0.00000  0.26143  1.68336



### Save NN and NP scores ###
NN_NP_scores <- cbind(Contrib_site_log_transformed[ , c("SiteCode", "SiteLongitude", "SiteLatitude", 
                                                    "SiteCountry", "SiteEcoregion","year",
                                                    "SiteMeanSST", "Biomass", "HDI", "gravtot2",
                                                    "MarineEcosystemDependency", "coral_imputation",
                                                    "mpa_name", "mpa_enforcement", "protection_status",
                                                    "mpa_iucn_cat")],
                      data.frame(NN_score = mean_NN, NP_score = mean_NP) )
save(NN_NP_scores, file = here::here("outputs", "NN_NP_score_wheighted_mean.Rdata"))



### Save the final dataframe with all variables ###
sites_all_variables <- cbind(Contrib_site_log_transformed,
                             data.frame(NN_score = mean_NN, NP_score = mean_NP) ) |> 
  dplyr::rename_with(.cols = c("Biomass","N_Recycling","P_Recycling",
                               "Low_Mg_Calcite", "High_Mg_Calcite", "Aragonite", 
                               "Monohydrocalcite", "Amorphous_Carbonate",
                               "Herbivores_Biomass", "Invertivores_Biomass",
                               "Piscivores_Biomass", "Available_Biomass",
                               "gravtot2"),
                     .fn = ~paste0("log(", .x, ")")) #informs which variables are log-transformed


write.csv(sites_all_variables, row.names = FALSE,
          file = here::here("outputs", "tropical_reef_contributions_final_table.csv"))
###


#--------------------------Plot NN and NP scores--------------------------------
plot_histogram <- function(data = NN_NP_scores, 
                           x=NN_NP_scores$NN_score, col= "forestgreen",
                           title="Nature for nature score site's evaluation"){
  ggplot(data, aes(x) ) +
    geom_histogram( bins = 40, color="grey60", fill= "white", linewidth=0.2)+
    
    ggnewscale::new_scale_fill() +
    geom_histogram( bins = 40, color="grey60", aes(fill = after_stat(x)), 
                    linewidth=0.2)+
    scale_fill_gradient2(low="white", mid="white", high = col ,
                         midpoint= quantile(x,0.5), 
                         limits =quantile(x, probs=c( 0.5,1))+ c(-0.1,0.1),
                         na.value="transparent")+
    
    ggnewscale::new_scale_fill() +
    geom_histogram( bins = 40, aes(fill=after_stat(x)), color="grey60", linewidth=0.2)+
    scale_fill_gradient2(low="black", mid="white", high="white",
                         midpoint= quantile(x,0.5), 
                         limits =quantile(x, probs=c(0,0.5)) + c(-0.1,0.1),
                         na.value="transparent")+
    
    labs(title = title)+
    geom_vline(xintercept = quantile(x, 0.25), linetype = 3, col= "black")+
    geom_vline(xintercept = quantile(x, 0.75), linetype = 3, col= "black")+
    theme_minimal()+
    theme(legend.position = "none")
} #end of plot_histogram function


## NN
histo_NN <- plot_histogram(data = NN_NP_scores, x=NN_NP_scores$NN_score, col= "forestgreen",
                           title="Nature for Nature score in reefs")
ggsave(filename = here::here("outputs", "figures","hist_NN_weighted_mean.png"), histo_NN, width = 8, height =6 )


##NP
histo_NP <- plot_histogram(data = NN_NP_scores, x=NN_NP_scores$NP_score, col= "dodgerblue3",
                           title="Nature to People score in reefs")
ggsave(filename = here::here("outputs", "figures","hist_NP_weighted_mean.png"), histo_NP, width = 8, height =6 )


### plot scores according to socio-envir variables
variables <- c("SiteLatitude", "SiteLongitude", "SiteMeanSST", "Biomass",
               "HDI", "MarineEcosystemDependency", "coral_imputation",
               "gravtot2", "year")

plot_score_according_variables <- function(data = NN_NP_scores,
                                           score = "NN_score",
                                           var = "SiteLatitude"){
  ggplot(data, aes(x = data[,var][[var]],y= data[,score][[score]], alpha = 0.5)) +
    geom_point()+
    geom_smooth() +
    labs(x = var, y = score , title = "") +
    theme_minimal()+
    theme(legend.position = 'none')
}

## NN
plots <- lapply( variables, function(var){
  plot_score_according_variables(data = NN_NP_scores, score = "NN_score", var = var)
})

plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]]+
  plots[[6]] + plots[[7]] + plots[[8]] + plots[[9]]+
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here::here("outputs", "figures","Socio_envir_variables_with_NN_scores.png"), plot, width = 22, height =14 )

## NP
plots <- lapply( variables, function(var){
  plot_score_according_variables(data = NN_NP_scores, score = "NP_score", var = var)
})

plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]]+
  plots[[6]] + plots[[7]] + plots[[8]] +plots[[9]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here::here("outputs", "figures","Socio_envir_variables_with_NP_scores.png"), plot, width = 22, height =14 )


##-------------NN and NP correlation -----------------
cor.test(NN_NP_scores$NP_score, NN_NP_scores$NN_score) # cor = 0.2442032 

l<- lm(NP_score ~ NN_score, NN_NP_scores)
summary(l)

# # #test NN and NP against contributions
# cor.test(sites_all_variables$Turnover_Available_Biomass, sites_all_variables$NP_score)
# cor.test(sites_all_variables$`log(Biomass)`, sites_all_variables$NP_score)
# 
# summary(lm(sites_all_variables$Turnover_Available_Biomass ~ sites_all_variables$NP_score))
# plot(sites_all_variables$Turnover_Available_Biomass ~ sites_all_variables$NP_score)
# 
# cor.test(sites_all_variables$Aesthetic, sites_all_variables$NP_score)
# cor.test(sites_all_variables$Public_Attention, sites_all_variables$NP_score)
# cor.test(sites_all_variables$Zinc, sites_all_variables$NP_score)
# cor.test(sites_all_variables$Iron, sites_all_variables$NP_score)
# cor.test(sites_all_variables$`log(Available_Biomass)`, sites_all_variables$NP_score)
# 
# cor.test(sites_all_variables$`log(Biomass)`, sites_all_variables$NN_score)
# summary(lm(sites_all_variables$`log(Biomass)`~ sites_all_variables$NN_score))
# plot(sites_all_variables$`log(Biomass)`~ sites_all_variables$NN_score)
# 
# df <- sites_all_variables[,c(8:37,48,49)]
# cor_matrix <- cor(df)
# corrplot::corrplot(cor_matrix, method = "circle")
# 
# # Compare with random variables
# set.seed(123)
# random_vars <- replicate(6, rnorm(1200))
# mean_score <- rowMeans(random_vars)
# summary(lm(random_vars[,1]~mean_score))



##-------------Study the NNxNP 2D space -----------------
#### Study the contribution distributions in the NNxNP space ####
var <- c('`log(Biomass)`', 'Turnover_Available_Biomass', 'Aesthetic', 'Iron', 
         'Omega_3','Taxonomic_Richness', '`log(Available_Biomass)`',
         'Functional_Entropy', 'Endemism', 'SiteLatitude')

# var <- c('`log(Biomass)`', 'Turnover_Available_Biomass')
# sites_all_variables$Turnover_Available_Biomass <- log(sites_all_variables$Turnover_Available_Biomass)

plot_contrib <- parallel::mclapply(var, mc.cores=5, function(contrib){
  ggplot(sites_all_variables, aes( y= NP_score, x = NN_score) ) +
    geom_point(size = 2, aes_string(colour= contrib), alpha = 0.7) +
    scale_color_gradientn(colours = rev(
      RColorBrewer::brewer.pal(10,name="RdYlBu")))+
    #add lines
    geom_vline(xintercept = 0, linetype = 1, col = "black", linewidth = 0.2)+
    geom_hline(yintercept = 0, linetype = 1, col = "black", linewidth = 0.2)+
    
    labs( x=  "Nature to Nature", y = "Nature to People")+
    theme_minimal()+
    theme(legend.position = c(0.85,0.1))
})

plot <- (plot_contrib[[1]] + plot_contrib[[2]] + plot_contrib[[3]]) / 
  (plot_contrib[[4]] + plot_contrib[[5]] + plot_contrib[[6]] ) /
  (plot_contrib[[7]] + plot_contrib[[8]] + plot_contrib[[9]] )
ggsave(here::here("outputs", "figures", "contributions_distributions_in_NNxNP_space.png"),
       plot , width=13, height = 13 )

# plot_contrib[[10]]


#### Define colors by quarter in the NNxNP space ####
NN_NP_with_product <- NN_NP_scores |>
  dplyr::bind_cols( protection = ifelse(NN_NP_scores$mpa_enforcement == "High" &
                                          stringr::str_detect(NN_NP_scores$protection_status, "No take"),
                                        "No take",
                                        ifelse(is.na(NN_NP_scores$mpa_name)==F, 
                                               "Restricted", "Fished"))) |>

    dplyr::mutate(NNxNP = abs(NN_score * NP_score)) |>
  dplyr::bind_cols(rank = rank(abs(NN_NP_scores$NN_score * NN_NP_scores$NP_score))) |>
  dplyr::bind_cols(up_right = ifelse(NN_NP_scores$NN_score > 0 & NN_NP_scores$NP_score > 0,1,NA)) |>
  dplyr::bind_cols(up_left = ifelse(NN_NP_scores$NN_score < 0 & NN_NP_scores$NP_score > 0,1,NA)) |>
  dplyr::bind_cols(down_right = ifelse(NN_NP_scores$NN_score > 0 & NN_NP_scores$NP_score < 0,1,NA)) |>
  dplyr::bind_cols(down_left = ifelse(NN_NP_scores$NN_score < 0 & NN_NP_scores$NP_score < 0,1,NA))

NN_NP_with_product <- NN_NP_with_product |>
  dplyr::bind_cols(rank_u_r = rank(NN_NP_with_product$NNxNP * NN_NP_with_product$up_right, na.last="keep")) |>
  dplyr::bind_cols(rank_u_l = rank(NN_NP_with_product$NNxNP * NN_NP_with_product$up_left, na.last="keep")) |>
  dplyr::bind_cols(rank_d_r = rank(NN_NP_with_product$NNxNP * NN_NP_with_product$down_right, na.last="keep")) |>
  dplyr::bind_cols(rank_d_l = rank(NN_NP_with_product$NNxNP * NN_NP_with_product$down_left, na.last="keep"))

save(NN_NP_with_product, file = here::here("outputs","NN_NP_score_wheighted_mean_quarter_ranked.Rdata"))


#### plot NN against NP, with color by quarter, and outlier names ####
NN_NP_plot <- ggplot(NN_NP_with_product, aes( y= NP_score, x = NN_score) ) +
  #up right quarter
  geom_point(data= dplyr::filter(NN_NP_with_product, up_right == 1),
             size = 4,  
             stroke=0,
             aes(fill= rank_u_r, shape = protection)) +
  scale_fill_gradient(name="up_right",
                      low = "grey90", high="darkgoldenrod3",
                      limits =quantile(NN_NP_with_product$rank_u_r,
                                       probs=c( 0, 1), na.rm=T),  ## Set c( 0.5, 1) to begin the color gradient from the dashed line
                      na.value=NA, breaks = seq(1,400, 10)) +
  guides(fill = "none") + ggnewscale::new_scale("fill") +
  
  #up left quarter
  geom_point(data= dplyr::filter(NN_NP_with_product, up_left == 1),
             size = 4, 
             stroke = 0,
             aes(fill= rank_u_l, shape = protection)) +
  scale_fill_gradient(name="up_left",
                      low = "grey90", high="dodgerblue3",
                      limits =quantile(NN_NP_with_product$rank_u_l,
                                       probs=c( 0,1), na.rm=T), na.value=NA) +
  guides(fill = "none") + ggnewscale::new_scale("fill") +
  
  #down right quarter
  geom_point(data= dplyr::filter(NN_NP_with_product, down_right == 1),
             size = 4, 
             stroke=0,
             aes(fill= rank_d_r, shape = protection)) +
  scale_fill_gradient(name="down_right",
                      low = "grey90", high="forestgreen",
                      limits =quantile(NN_NP_with_product$rank_d_r,
                                       probs=c( 0,1), na.rm=T), na.value=NA) +
  guides(fill = "none") + ggnewscale::new_scale("fill") +
  
  
  #down left quarter
  geom_point(data= dplyr::filter(NN_NP_with_product, down_left == 1),
             size = 4,
             stroke=0,
             aes(fill= rank_d_l, shape = protection)) +
  scale_fill_gradient(name="down_left",
                      low = "grey90", high="grey30",
                      limits =quantile(NN_NP_with_product$rank_d_l,
                                       probs=c( 0,1), na.rm=T), na.value=NA) +
  guides(fill = "none") +
  
  
  
  # Add outliers
  geom_point(data = NN_NP_with_product[ c(
    which(NN_NP_with_product$rank_u_r >= 
            quantile(NN_NP_with_product$rank_u_r, probs=c(0.95), na.rm=T)) ,
    which(NN_NP_with_product$rank_d_l >= 
            quantile(NN_NP_with_product$rank_d_l, probs=c(0.95), na.rm=T)) ,
    which(NN_NP_with_product$rank_d_r >=
            quantile(NN_NP_with_product$rank_d_r, probs=c(0.95), na.rm=T)) ,
    which(NN_NP_with_product$rank_u_l >=
            quantile(NN_NP_with_product$rank_u_l, probs=c(0.95), na.rm=T))
  ), ],
  aes(y= NP_score, x = NN_score, shape = protection),
  size = 4,  
  stroke = 1)+
  
  # Add names
  ggrepel::geom_text_repel( aes(label=paste(SiteCode, ":", SiteEcoregion)), size=3, 
                            force_pull = 0,
                            force=0.5,
                            data = NN_NP_with_product[ c(
                              which(NN_NP_with_product$rank_u_r >= 
                                      quantile(NN_NP_with_product$rank_u_r, probs=c(0.95), na.rm=T)) ,
                              which(NN_NP_with_product$rank_d_l >= 
                                      quantile(NN_NP_with_product$rank_d_l, probs=c(0.95), na.rm=T)) ,
                              which(NN_NP_with_product$rank_d_r >=
                                      quantile(NN_NP_with_product$rank_d_r, probs=c(0.95), na.rm=T)) ,
                              which(NN_NP_with_product$rank_u_l >=
                                      quantile(NN_NP_with_product$rank_u_l, probs=c(0.95), na.rm=T))
                            ), ] )+
  
  
  # see MPAs
  scale_shape_manual(values=c(21,24,23))+
  
  # add 50% square: equation: abs(x)^p + abs(y)^p = 1 -> p=0.5
  geom_function(aes(x= NN_score, y = NP_score),
                fun = function(x){(1-abs(x)^0.5)^(1/0.5) },
                xlim=c(-1,1), linetype=3, linewidth=0.5) +
  geom_function(aes(x= NN_score, y = NP_score),
                fun = function(x){-(1-abs(x)^0.5)^(1/0.5) },
                xlim=c(-1,1), linetype=3, linewidth=0.5) +
  
  #add axes
  geom_vline(xintercept = 0, linetype = 1, col = "black", linewidth = 0.2)+
  geom_hline(yintercept = 0, linetype = 1, col = "black", linewidth = 0.2)+
  
  
  xlim(c(-2,2)) +
  ylim(c(-1.8,1.8)) +
  labs( x=  "", y = "")+
  theme_bw(base_line_size = 0)+
  theme(
    axis.text = element_text(size=13),
    legend.position = c(0.91,0.12),
    legend.background = element_rect(colour="black", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1.3,1,1.5,0.7, unit='cm'))+
  guides(color = "none", shape = guide_legend("Protection status")) + #shape = guide_legend("Protection status"))
  theme(legend.position = c("bottom"))

NN_NP_plot


## Add gradient in the margins of NN_NP_plot & names the corners ##
# /!\ This plot need the figures "legend_gradient" of NN and NP which are created
# in the script 2b_make_fig_2.R 
grad_NN <- cowplot::ggdraw() +
  cowplot::draw_image(here::here("outputs", "figures","legend_gradient_NN.png"))
grad_NP <- cowplot::ggdraw() +
  cowplot::draw_image(here::here("outputs", "figures","legend_gradient_NP.png"))

label_corner <- data.frame(
  X = c(-Inf,-Inf,Inf,Inf),
  Y =  c(-Inf, Inf,-Inf,Inf),
  text = c("Dark spots","NP only",
           "NN only","Bright spots"),
  x_adjust = c(-0.2,-0.2,1.5,1.3),
  y_adjust = c(-2,3,-2,3),
  col = c("grey30", "dodgerblue3", "forestgreen", "darkgoldenrod3"))

fig_S8 <- NN_NP_plot +
  geom_text(data=label_corner, 
            aes(x=X,y=Y, hjust= x_adjust,vjust=y_adjust,label=text, 
                fontface= "bold", color = col),
            size = 4.5)+
  scale_color_identity() +
  annotation_custom( ggplotGrob(grad_NN),
                     xmin= -2.1 ,#-2.2,#-2,
                     ymin = -2.6 ,#-2.5,  #-2.4,
                     xmax=2.1,  #2.2,#1,
                     ymax = -2.2 )+
  
  annotation_custom( ggplotGrob(grad_NP),
                     xmin=-2.55, #-2.7
                     xmax=-2.32, # -2.3
                     ymin = -2, #-2
                     ymax = 1.97)+
  labs( x=  "Nature for Nature", y = "Nature for People")+
  theme( legend.position = "bottom",
         legend.box.margin = margin(30, 0, 0,0), 
         axis.title.x = element_text(hjust = 0.96,
                                     vjust = -4.5,
                                     colour = "white",
                                     face = "bold",
                                     size = 15),
         axis.title.y = element_text(hjust = 0.96,
                                     vjust = 2.1,
                                     colour = "white",
                                     face = "bold",
                                     size = 15)) 


fig_S8
ggsave(plot=fig_S8, filename=here::here("outputs", "figures","Sites and outlier names in NN and NP scores.png"),
       height = 8.5, width = 11)



#### NN and NP together on map ####
function_NN_NP_on_map <- function(coord_NN_NP = NN_NP_with_product, ylim = c(-36, 31),
                                  xlim= c(-180,180), title=""){
  ggplot(coord_NN_NP) +
    geom_sf(data = coast, color = NA, fill = "lightgrey") +
    
    #up left quarter
    geom_point(data= dplyr::filter(coord_NN_NP,  up_left == 1),
               position = position_jitter(width =width_jitter, height =height_jitter),
               size = 2, alpha = 0.7,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="up_left",
                          low = "white", high="dodgerblue3",
                          limits =quantile(coord_NN_NP$rank * coord_NN_NP$up_left,
                                           probs=c( .3,1), na.rm=T), na.value=NA) +
    geom_point(aes( x= SiteLongitude, y = SiteLatitude),
               position = position_jitter(width =width_jitter, height =height_jitter),
               size = 2, stroke = 1, shape = 1,
               color= "darkblue",
               data = coord_NN_NP[which(coord_NN_NP$rank_u_l >=
                                          quantile(coord_NN_NP$rank_u_l, probs=c(0.95), na.rm=T)), ] )+
    
    ggnewscale::new_scale("colour") +
    
    #down right quarter
    geom_point(data= dplyr::filter(coord_NN_NP,  down_right == 1),
               position = position_jitter(width =width_jitter, height =height_jitter),
               size = 2, alpha = 0.7,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="down_right",
                          low = "white", high="forestgreen",
                          limits =quantile(coord_NN_NP$rank * coord_NN_NP$down_right,
                                           probs=c( 0.3,1), na.rm=T), na.value=NA) +
    geom_point(aes( x= SiteLongitude, y = SiteLatitude),
               position = position_jitter(width =width_jitter, height =height_jitter),
               size = 2, stroke = 1, shape = 1,
               color= "darkgreen",
               data = coord_NN_NP[which(coord_NN_NP$rank_d_r >=
                                          quantile(coord_NN_NP$rank_d_r, probs=c(0.95), na.rm=T)), ] )+
    
    ggnewscale::new_scale("colour") +
    
    #up right quarter
    geom_point(data= dplyr::filter(coord_NN_NP,  up_right == 1),
               position = position_jitter(width =width_jitter, height =height_jitter),
               size = 2, alpha = 0.7,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="up_right",
                          low = "white", high="darkgoldenrod3",
                          limits =quantile(coord_NN_NP$rank * coord_NN_NP$up_right,
                                           probs=c( 0.3,1), na.rm=T), na.value=NA) +
    geom_point(aes( x= SiteLongitude, y = SiteLatitude),
               position = position_jitter(width =width_jitter, height =height_jitter),
               size = 2, stroke = 1, shape = 1,
               color= "darkgoldenrod4",
               data = coord_NN_NP[which(coord_NN_NP$rank_u_r >=
                                          quantile(coord_NN_NP$rank_u_r, probs=c(0.95), na.rm=T)), ] )+
    ggnewscale::new_scale("colour") +
    
    
    #down left quarter
    geom_point(data= dplyr::filter(coord_NN_NP,  down_left == 1),
               position = position_jitter(width =width_jitter, height =height_jitter),
               size = 2, alpha = 0.7,
               aes(x = SiteLongitude, y = SiteLatitude, colour= rank)) +
    scale_colour_gradient(name="down_left",
                          low = "white", high="grey20",
                          limits =quantile(coord_NN_NP$rank * coord_NN_NP$down_left,
                                           probs=c( 0.3,1), na.rm=T), na.value=NA) +
    geom_point(aes( x= SiteLongitude, y = SiteLatitude),
               position = position_jitter(width =width_jitter, height =height_jitter),
               size = 2, stroke = 1, shape = 1,
               color= "grey5",
               data = coord_NN_NP[which(coord_NN_NP$rank_d_l >=
                                          quantile(coord_NN_NP$rank_d_l, probs=c(0.95), na.rm=T)), ] )+

  coord_sf(xlim= xlim, ylim = ylim, expand = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_minimal()+
    labs(title = title,
         x="", y= "") +
    theme(legend.position = "none",
          plot.title = element_text(size=10, face="bold"),
          # axis.text.x = element_blank(),
          # axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
}


width_jitter = 0.4
height_jitter = 0.4

world_map_zoom <- function_NN_NP_on_map( NN_NP_with_product, 
                                         ylim = c(-36, 31),xlim= c(-180,180),title = "") + 
  geom_rect(aes(xmin = 110, xmax = 160, ymin = -32, ymax = -7), color = "black", fill= "transparent")+
  geom_rect(aes(xmin = -95, xmax = -67, ymin = -3, ymax = 18), color = "black", fill= "transparent")


# gold_coast_map_NN_NP <- function_NN_NP_on_map( NN_NP_with_product, 
#                                                ylim = c(-25,-10),
#                                                xlim= c(140,160),
#                                                title = "")
caraib_map_NN_NP <- function_NN_NP_on_map( NN_NP_with_product, 
                                           ylim = c(-3, 18),
                                           xlim= c(-95,-67),
                                           title = "")
australia_map_NN_NP <- function_NN_NP_on_map( NN_NP_with_product, 
                                              ylim = c(-32, -7),
                                              xlim= c(110,160),
                                              title = "")


ggpubr::ggarrange(world_map_zoom, # First row with world map
                  ggpubr::ggarrange(caraib_map_NN_NP, australia_map_NN_NP, 
                                    ncol = 2, labels = c("B", "C")), # Second row with zooms
                  nrow = 2, labels = "A") 
ggsave(plot = last_plot(), width=11, height =7,
       filename = here::here("outputs", "figures", "world map of NN and NP score with zoom.png"))


##-------------Spatial distribution of NN and NP scores, independently -----------------
map_NN_or_NP <- function(coord_NN_NP = NN_NP_with_product,
                         Contrib = NN_NP_with_product$NP_score ,
                         col_Contrib= "dodgerblue3",
                         ylim = c(-36, 31),
                         xlim= c(-180,180), title=""){
  
  coord_NN_NP_sorted <- coord_NN_NP#[order(abs(Contrib)), ]
  # Contrib <- Contrib[order(abs(Contrib))]
  
  ggplot(coord_NN_NP_sorted) +
    geom_sf(data = coast, color = NA, fill = "lightgrey") +
    geom_point(aes(x = SiteLongitude, y = SiteLatitude,
                   color = rank(Contrib),
                   alpha= rank(abs(Contrib)),
                   size = 2)) +
    geom_point(aes(x = SiteLongitude, y = SiteLatitude),
               shape = 1, size = 2, stroke = 0.5,
               color= "black",
               data = coord_NN_NP_sorted[which(Contrib > quantile(Contrib, .98)),]) +
    
    scale_colour_gradientn(colours = colorRampPalette(rev(c( col_Contrib ,"white", "grey50")))(1000)) +
    scale_alpha_continuous(range = c(0, 0.6)) +
    
    coord_sf(xlim=xlim, ylim = ylim, expand = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    scale_size_continuous(range = c(0.5, 4), guide = "none") +
    theme_minimal()+
    labs(title = title,
         x="Longitude", y= "Latitude") +
    theme(
      legend.position = "none",
      plot.title = element_text(colour = col_Contrib, 
                                face = "bold", size = 12,
                                margin=margin(t = 20, b = -15),
                                hjust = 0.01))
}

## NN
NN_worldwide <- map_NN_or_NP(coord_NN_NP = NN_NP_with_product,
                             Contrib = NN_NP_with_product$NN_score,
                             col_Contrib= "forestgreen",
                             ylim = c(-36, 31),
                             xlim= c(-180,180), title="Nature to Nature")

## NP
NP_worldwide <- map_NN_or_NP(coord_NN_NP = NN_NP_with_product,
                             Contrib = NN_NP_with_product$NP_score,
                             col_Contrib= "dodgerblue3",
                             ylim = c(-36, 31),
                             xlim= c(-180,180), title="Nature to People")

#panel
ggpubr::ggarrange(NN_worldwide, NP_worldwide,
                  nrow = 2, labels = c("A", "B")) 
ggsave(plot = last_plot(), width=11, height =6,
       filename = here::here("outputs", "figures", "world map with NN and NP score SEPARATLY.png"))



##------------------------Study MPAs distribution in the NNxNP space------------------------------------
mpa_distribution <- NN_NP_with_product |>
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
#   Chisq = 28.164, df = 6, p-value = 8.753e-05

# Effect size
lsr::cramersV(contingency_table) #0.106695
