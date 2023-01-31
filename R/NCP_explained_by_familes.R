################################################################################
##
## Looking for correlation between families and NN or NS score
##
## NCP_explained_by_familes.R
##
## 23/01/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse", "ggplot2", "sf", "patchwork")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load(here::here("outputs","NN_NS_score_wheighted_mean.Rdata"))
load(here::here("biodiversity", "outputs", "occurrence_matrix_family_survey_relative_biomass.Rdata"))
load(here::here("biodiversity", "outputs", "occurrence_matrix_family_survey_01.Rdata"))
load(here::here("biodiversity", "outputs", "occurrence_matrix_nb_sp_per_family_survey.Rdata"))
load(here::here("data","metadata_surveys.Rdata"))

coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

##-------------preping data-------------
#Sum up relative biomass of each family at the site scale
surveys_family_pbiom <- as.data.frame(surveys_family_pbiom) %>%
  tibble::rownames_to_column("SurveyID") %>%
  dplyr::left_join(metadata_surveys[,c("SurveyID", "SiteCode")]) %>%
  dplyr::select(- SurveyID) 

site_family_pbiom <- surveys_family_pbiom %>%
  group_by( SiteCode) %>%
  summarise(across(.cols = everything(),
                   .fns = mean, 
                   .names = "{.col}")) %>%
  dplyr::filter(SiteCode %in% NN_NS_scores$SiteCode) %>%
  dplyr::select( - SiteCode)

#Sum up richness of each family at the site scale
surveys_family_richness <- as.data.frame(surveys_nb_sp_per_family) %>%
  tibble::rownames_to_column("SurveyID") %>%
  dplyr::left_join(metadata_surveys[,c("SurveyID", "SiteCode")]) %>%
  dplyr::select(- SurveyID) 

site_family_richness <- surveys_family_richness %>%
  group_by( SiteCode) %>%
  summarise(across(.cols = everything(),
                   .fns = max, 
                   .names = "{.col}")) %>%
  dplyr::filter(SiteCode %in% NN_NS_scores$SiteCode) %>%
  dplyr::select( - SiteCode)

#Sum up occurrences of each family at the site scale
surveys_family_occ <- as.data.frame(surveys_family_occ) %>%
  tibble::rownames_to_column("SurveyID") %>%
  dplyr::left_join(metadata_surveys[,c("SurveyID", "SiteCode")]) %>%
  dplyr::select(- SurveyID) 

site_family_occ <- surveys_family_occ %>%
  group_by( SiteCode) %>%
  summarise(across(.cols = everything(),
                   .fns = max, 
                   .names = "{.col}")) %>%
  dplyr::filter(SiteCode %in% NN_NS_scores$SiteCode) %>%
  dplyr::select( - SiteCode)


##------------- families distribution-------------
plot_distribution <- function(family, data){
  col <- fishualize::fish(n = ncol(data), option = "Ostracion_whitleyi", begin = 0, end = 0.8)
  names(col) <- colnames(data)
  
  ggplot(data) +
    aes(x = data[,family][[family]]) +
    geom_histogram(bins = 40L,
                   fill = col[family][[1]],
                   col = "black") +
    labs(title = family) +
    xlab("") + ylab("") +
    theme_minimal() +
    theme( legend.position = "none")
}

plots <- lapply(colnames(site_family_richness), FUN = plot_distribution, data = site_family_richness )

all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))
all_plot

ggsave(filename = here("outputs", "figures","families_distribution.png"), all_plot, width = 22, height =14 )

# 
# ###log distribution
# families_log_transformed <- site_family_pbiom |>
#   dplyr::mutate(across(.cols = everything(),
#                        .fns = ~ .x +1 , .names = "{.col}")) |>      # Adds 1 to values to log transformed
#   dplyr::mutate(across(.cols = everything(),
#                        .fns = log10 , .names = "{.col}")) |>       # log(x+1) to avoid negative values
#   dplyr::rename_with(.cols = everything(),
#                      .fn = ~ paste0("log(", .x, ")"))


##------------- plot families according to NN and NS scores-------------
plot_correlation <- function(data,x,y,i, y_title){  
  ggplot() +
    geom_point(aes(x = data[,x][[x]], y = y),
               color = col[i], alpha = 0.6, size = 1) +
    stat_smooth(data = data, aes(x = data[,x][[x]], y = y),
                method="lm", se=FALSE, linetype="dashed", 
                size=0.8, color = "black") +
    ggpubr::stat_cor(data = data, aes(x = data[,x][[x]], y = y))+
    theme_bw() +
    labs(x = x, y = y_title) +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 17))
}
## Families according to NN
col <- fishualize::fish(n = ncol(site_family_richness), option = "Coryphaena_hippurus", begin = 0.5, end = 0.8)
### family richness
plots <- lapply( 1:ncol(site_family_richness), function(i){
  plot_correlation(data= site_family_richness ,
                         x = colnames(site_family_richness)[i],
                         y = NN_NS_scores$NN_score,
                         i, y_title = "")
})

plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] + 
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] +
  
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))
plot
ggsave(filename = here::here("outputs", "figures","families_richness_according_to_NN_score.png"), plot, width = 22, height =14 )

### family biomass
plots <- lapply( 1:ncol(site_family_pbiom), function(i){
  plot_correlation(data= site_family_pbiom ,
                   x = colnames(site_family_richness)[i],
                   y = NN_NS_scores$NN_score,
                   i, y_title = "")
})
  
plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] + 
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] +
  
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))
plot
ggsave(filename = here::here("outputs", "figures","families_pbiom_according_to_NN_score.png"), plot, width = 22, height =14 )


## Families according to NS
col <- fishualize::fish(n = ncol(site_family_richness), option = "Ostracion_whitleyi", begin = 0.5, end = 0.8)
### family richness
plots <- lapply( 1:ncol(site_family_richness), function(i){
  plot_correlation(data= site_family_richness ,
                   x = colnames(site_family_richness)[i],
                   y = NN_NS_scores$NN_score,
                   i, y_title = "")
})
plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] + 
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] +
  
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))
plot
ggsave(filename = here::here("outputs", "figures","families_richness_according_to_NS_score.png"), plot, width = 22, height =14 )

### family pbiom
plots <- lapply( 1:ncol(site_family_pbiom), function(i){
  plot_correlation(data= site_family_pbiom ,
                   x = colnames(site_family_richness)[i],
                   y = NN_NS_scores$NN_score,
                   i, y_title = "")
})

plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] + 
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] +
  
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))
plot
ggsave(filename = here::here("outputs", "figures","families_pbiom_according_to_NS_score.png"), plot, width = 22, height =14 )


## Families according to Btot 
col <- rep("grey60", ncol(site_family_richness))
plots <- lapply( 1:ncol(site_family_richness), function(i){
  plot_correlation(data= site_family_richness ,
                   x = colnames(site_family_richness)[i],
                   y = NN_NS_scores$NN_score,
                   i, y_title = "")
})
plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] + 
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] +
  
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))
plot
ggsave(filename = here::here("outputs", "figures","families_richness_according_to_Btot.png"), plot, width = 22, height =14 )


##------------- test multi regression-------------
plot_hist <- function(data, title = ""){
  ggplot() +
    aes(x = as.numeric(data)) +
    geom_histogram(bins = 30L,
                   fill= "gray",col = "black") +
    labs(title = title) +
    xlab("") + ylab("") +
    theme_minimal() +
    theme( legend.position = "none")
}
#Families occurrences
site_family <- site_family_occ
multiple_regression_occ_NN <- lm(NN_NS_scores$NN_score ~ 
                                   site_family$Pomacentridae + site_family$Labridae + site_family$Lethrinidae +  
                                   site_family$Mullidae + site_family$Scaridae + site_family$Haemulidae +   
                                   site_family$Kyphosidae + site_family$Acanthuridae + site_family$Monacanthidae +
                                   site_family$Scorpaenidae + site_family$Chaetodontidae + site_family$Tetraodontidae+
                                   site_family$Pomacanthidae + site_family$Serranidae  + site_family$Ostraciidae +  
                                   site_family$Holocentridae + site_family$Lutjanidae + site_family$Sciaenidae +   
                                   site_family$Balistidae + site_family$Fistulariidae + site_family$Pempheridae +  
                                   site_family$Cirrhitidae + site_family$Siganidae + site_family$Zanclidae +    
                                   site_family$Mugilidae + site_family$Bothidae )
plot(multiple_regression_occ_NN)
summary(multiple_regression_occ_NN)
plot_hist(multiple_regression_occ_NN$residuals, title = "lm on occurences")


multiple_regression_occ_NS <- lm(NN_NS_scores$NS_score ~ 
                            site_family$Pomacentridae + site_family$Labridae + site_family$Lethrinidae +  
                            site_family$Mullidae + site_family$Scaridae + site_family$Haemulidae +   
                            site_family$Kyphosidae + site_family$Acanthuridae + site_family$Monacanthidae +
                            site_family$Scorpaenidae + site_family$Chaetodontidae + site_family$Tetraodontidae+
                            site_family$Pomacanthidae + site_family$Serranidae  + site_family$Ostraciidae +  
                            site_family$Holocentridae + site_family$Lutjanidae + site_family$Sciaenidae +   
                            site_family$Balistidae + site_family$Fistulariidae + site_family$Pempheridae +  
                            site_family$Cirrhitidae + site_family$Siganidae + site_family$Zanclidae +    
                            site_family$Mugilidae + site_family$Bothidae )
plot(multiple_regression_occ_NS)
summary(multiple_regression_occ_NS)
plot_hist(multiple_regression_occ_NS$residuals, title = "lm on occurences")

#Families richness
site_family <- site_family_richness

multiple_regression_richness_NN <- lm(NN_NS_scores$NN_score ~ 
                                     site_family$Pomacentridae + site_family$Labridae + site_family$Lethrinidae +  
                                     site_family$Mullidae + site_family$Scaridae + site_family$Haemulidae +   
                                     site_family$Kyphosidae + site_family$Acanthuridae + site_family$Monacanthidae +
                                     site_family$Scorpaenidae + site_family$Chaetodontidae + site_family$Tetraodontidae+
                                     site_family$Pomacanthidae + site_family$Serranidae  + site_family$Ostraciidae +  
                                     site_family$Holocentridae + site_family$Lutjanidae + site_family$Sciaenidae +   
                                     site_family$Balistidae + site_family$Fistulariidae + site_family$Pempheridae +  
                                     site_family$Cirrhitidae + site_family$Siganidae + site_family$Zanclidae +    
                                     site_family$Mugilidae + site_family$Bothidae )
plot(multiple_regression_richness_NN)
summary(multiple_regression_richness_NN)
plot_hist(multiple_regression_richness_NN$residuals, title = "lm on richness")


multiple_regression_richness_NS <- lm(NN_NS_scores$NS_score ~ 
                            site_family$Pomacentridae + site_family$Labridae + site_family$Lethrinidae +  
                            site_family$Mullidae + site_family$Scaridae + site_family$Haemulidae +   
                            site_family$Kyphosidae + site_family$Acanthuridae + site_family$Monacanthidae +
                            site_family$Scorpaenidae + site_family$Chaetodontidae + site_family$Tetraodontidae+
                            site_family$Pomacanthidae + site_family$Serranidae  + site_family$Ostraciidae +  
                            site_family$Holocentridae + site_family$Lutjanidae + site_family$Sciaenidae +   
                            site_family$Balistidae + site_family$Fistulariidae + site_family$Pempheridae +  
                            site_family$Cirrhitidae + site_family$Siganidae + site_family$Zanclidae +    
                            site_family$Mugilidae + site_family$Bothidae )
plot(multiple_regression_richness_NS)
summary(multiple_regression_richness_NS)
plot_hist(multiple_regression_richness_NS$residuals, title = "lm on richness")

#Families relative biomass
site_family <- site_family_pbiom

multiple_regression_pbiom_NN <- lm(NN_NS_scores$NN_score ~ 
                                  site_family$Pomacentridae + site_family$Labridae + site_family$Lethrinidae +  
                                  site_family$Mullidae + site_family$Scaridae + site_family$Haemulidae +   
                                  site_family$Kyphosidae + site_family$Acanthuridae + site_family$Monacanthidae +
                                  site_family$Scorpaenidae + site_family$Chaetodontidae + site_family$Tetraodontidae+
                                  site_family$Pomacanthidae + site_family$Serranidae  + site_family$Ostraciidae +  
                                  site_family$Holocentridae + site_family$Lutjanidae + site_family$Sciaenidae +   
                                  site_family$Balistidae + site_family$Fistulariidae + site_family$Pempheridae +  
                                  site_family$Cirrhitidae + site_family$Siganidae + site_family$Zanclidae +    
                                  site_family$Mugilidae + site_family$Bothidae )
plot(multiple_regression_pbiom_NN)
summary(multiple_regression_pbiom_NN)
plot_hist(multiple_regression_pbiom_NN$residuals, title = "lm on pbiom")


multiple_regression_pbiom_NS <- lm(NN_NS_scores$NN_score ~ 
                            site_family$Pomacentridae + site_family$Labridae + site_family$Lethrinidae +  
                            site_family$Mullidae + site_family$Scaridae + site_family$Haemulidae +   
                            site_family$Kyphosidae + site_family$Acanthuridae + site_family$Monacanthidae +
                            site_family$Scorpaenidae + site_family$Chaetodontidae + site_family$Tetraodontidae+
                            site_family$Pomacanthidae + site_family$Serranidae  + site_family$Ostraciidae +  
                            site_family$Holocentridae + site_family$Lutjanidae + site_family$Sciaenidae +   
                            site_family$Balistidae + site_family$Fistulariidae + site_family$Pempheridae +  
                            site_family$Cirrhitidae + site_family$Siganidae + site_family$Zanclidae +    
                            site_family$Mugilidae + site_family$Bothidae )
plot(multiple_regression_pbiom_NS)
summary(multiple_regression_pbiom_NS)
plot_hist(multiple_regression_pbiom_NS$residuals, title = "lm on pbiom")

