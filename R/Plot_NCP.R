################################################################################
##
## Display and record different ways to vizualize Nature based contributions
##
## Plot_NCP.R
##
## 24/01/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse", "ggplot2", "sf", "patchwork", "fishualize", "corrplot")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load(here::here("data","metadata_surveys.Rdata"))
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

load(here::here("outputs","all_NCP_site.Rdata"))
load(here::here("outputs","all_NCP_site_log_transformed.Rdata"))
# load(here::here("outputs","NCP_site_log_coral_reef.Rdata"))
# load(here::here("outputs","NCP_site_log_SST20.Rdata"))
# load(here::here("outputs","NCP_site_log_coral_5_imputed.Rdata"))
# load(here::here("outputs","NCP_site_log_wo_australia.Rdata"))
# NCP_site <- NCP_site_condition

# Preping data  
## Clean data  
NCP_site_clean_before_log <- subset(NCP_site, 
                         select = -c(SiteCode, SiteCountry, SurveyDate,
                                     SiteEcoregion, SurveyDepth, 
                                     SiteMeanSST, SiteLatitude, SiteLongitude,
                                     HDI, MarineEcosystemDependency,
                                     coral_imputation, gravtot2, mpa_name,
                                     mpa_enforcement, protection_status, 
                                     mpa_iucn_cat))
NCP_site_clean <- subset(NCP_site_log_transformed, 
                         select = -c(SiteCode, SiteCountry, SurveyDate,
                                     SiteEcoregion, SurveyDepth, 
                                     SiteMeanSST, SiteLatitude, SiteLongitude,
                                     HDI, MarineEcosystemDependency,
                                     coral_imputation, gravtot2, mpa_name,
                                     mpa_enforcement, protection_status, 
                                     mpa_iucn_cat))


##-------------Histograms of NCPs-------------
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

library(ggplot2)
library(patchwork)
plots <- lapply(colnames(NCP_site_clean_before_log), FUN = plot_distribution, data = NCP_site_clean_before_log )
all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + plots[[30]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here::here("outputs", "figures","NCP_distribution.png"), all_plot, width = 22, height =14 )


#### NCPs distribution with log correction
plots <- lapply(colnames(NCP_site_clean), FUN = plot_distribution, data = NCP_site_clean )

all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + plots[[30]]+
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here::here("outputs", "figures","NCP_log_transformed_distribution.png"), all_plot, width = 22, height =14 )



##-------------Correlations between NCPs-------------
  plot_correlation <- function(x,y,i){  
    ggplot() +
      geom_point(aes(y = NCP_site[,y][[y]], x = x),
                 color = col[i], alpha = 0.6, size = 1) +
      theme_bw() +
      labs(x = x_title, y = y) +
      theme(panel.grid.minor = element_blank(),
            axis.text = element_text(color = "black"), 
            axis.title = element_text(size = 17))
  }
  
  
  ## With Biomass
  x<- NCP_site$Biomass
  x_title = "total Biomass"
  col <- fishualize::fish(n = ncol(NCP_site_clean), option = "Ostracion_whitleyi", begin = 0, end = 0.8)
  
  plots <- lapply( 1:ncol(NCP_site_clean), function(i){
    plot_correlation(x,colnames(NCP_site_clean)[i],i)
  })
  
  plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
    plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
    plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
    plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
    plots[[28]] + plots[[29]] + plots[[30]] +
    
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  ggsave(filename = here::here("outputs", "figures","NCP_correlation_with_biomass.png"), plot, width = 22, height =14 )
  
  
  
  ## With biodiversity
  x<- NCP_site$Taxonomic_Richness
  x_title = "taxonomic richness"
  col <- fishualize::fish(n = ncol(NCP_site_clean), option = "Ostracion_whitleyi", begin = 0, end = 0.8)
  
  plots <- lapply( 1:ncol(NCP_site_clean), function(i){
    plot_correlation(x,colnames(NCP_site_clean)[i],i)
  })
  
  plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
    plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
    plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
    plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
    plots[[28]] + plots[[29]] + plots[[30]] +
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  ggsave(filename = here::here("outputs", "figures","NCP_correlation_with_biodiversity.png"), plot, width = 22, height =14 )
  
  
  
  #### Corr-matrix for all NCPs
  png(filename = here::here("outputs", "figures","corr_matrix.png"), 
      width= 40, height = 30, units = "cm", res = 1000)
  print({
    M <- cor(NCP_site_clean_before_log)
    # ## circle + black number
    corrplot::corrplot(M, order = 'AOE', type = 'lower', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r')
    corrplot::corrplot(M, add = TRUE, type = 'upper', method = 'number', order = 'AOE', insig = 'p-value',
             diag = FALSE, tl.pos = 'n', cl.pos = 'n')
    
    # corrplot(M, p.mat = testRes$p, insig = 'p-value')
    # corrplot(M, order = 'hclust', type = 'lower', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r')
    # corrplot(M, order = 'hclust', add= T, type= 'upper', p.mat = testRes$p, insig = 'p-value', 
    #          tl.pos= 'n', cl.pos = 'n')
  })
  dev.off() 
  
  #### Corr-matrix for log transformed NCPs
  png(filename = here::here("outputs", "figures","corr_matrix_log_transformed_NCP.png"), 
      width= 40, height = 30, units = "cm", res = 1000)
  print({
    M <- cor(NCP_site_clean)
    corrplot::corrplot(M, order = 'AOE', type = 'lower', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r')
    corrplot::corrplot(M, add = TRUE, type = 'upper', method = 'number', order = 'AOE', insig = 'p-value',
             diag = FALSE, tl.pos = 'n', cl.pos = 'n')
  })
  dev.off() 
  
##-------------Links of socio-envir with biomass -------------
  socio_envir <- c("HDI", "MarineEcosystemDependency", "coral_imputation",
                   "SiteLatitude", "mpa_iucn_cat")
  
  x <- as.list(NCP_site_clean[,"Biomass"])[["Biomass"]]
  x_title = "log(total Biomass)"
  col <- fishualize::fish(n = length(socio_envir), option = "Centropyge_loricula", begin = 0, end = 0.8)
  
  plots <- lapply( 1:length(socio_envir), function(i){
    plot_correlation(x,socio_envir[i],i) +
    coord_flip() 
  })
  
  plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]]+
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  ggsave(filename = here::here("outputs", "figures","Socio_envir_variables_with_biomass.png"), plot, width = 22, height =14 )
  
  
  
##-------------plot NCPs on map-------------
  #plot function
  plot_NCP_on_world_map <- function(NCP = "Taxonomic_Richness", xlim=c(-180,180), ylim = c(-36, 31),
                                    title="world map with "){
    library(ggplot2)
    map <- ggplot(NCP_site_log_transformed) +
      geom_sf(data = coast, color = "grey30", fill = "lightgrey",
              aes(size=0.1)) +
      
      geom_point(data=NCP_site_log_transformed,
                 size = 4, shape = 20,
                 aes(x = SiteLongitude, y = SiteLatitude,
                     colour= NCP_site_log_transformed[,NCP][[1]],
                     alpha = scale(NCP_site_log_transformed[,NCP][[1]]))) +
      scale_colour_gradient(NCP,
                            low = "dodgerblue", high="darkred",
                            na.value=NA) +
      
      coord_sf(xlim, ylim, expand = FALSE) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      scale_size_continuous(range = c(0.5, 4), guide = "none") +
      theme_minimal()+
      labs(title = paste0(NCP, " geographic distribution"),
           x="", y= "") +
      theme(legend.position = "right",
            plot.title = element_text(size=10, face="bold"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
      )
    
    ggsave( here::here("outputs", "figures", "NCP_on_world_map", paste0( title , NCP, ".jpg")), plot = map, width=15, height = 7 )
    #map
  }
  
  # save maps
parallel::mclapply(colnames(NCP_site_clean), mc.cores=15, function(NCP){
    plot_NCP_on_world_map(NCP, xlim=c(-180,180), ylim = c(-36, 31), title="world_map_with_")
})
  
  #focus on Biomass
  plot_NCP_on_world_map(NCP= "Biomass", ylim = c(-39, 0), xlim= c(100,180), title= "Australian_map_with_")
  plot_NCP_on_world_map(NCP= "Biomass", ylim = c(-5, 30), xlim= c(-100,-55), title= "Caraib_map_with_")
  
  #focus on robustness
  plot_NCP_on_world_map(NCP= "Trophic_web_robustness", ylim = c(-39, 0), xlim= c(100,180), title= "Australian_map_with_")
  plot_NCP_on_world_map(NCP= "Trophic_web_robustness", ylim = c(-5, 30), xlim= c(-100,-55), title= "Caraib_map_with_")
  
  #focus on mean TL
  plot_NCP_on_world_map(NCP= "mean_Trophic_Level", ylim = c(-39, 0), xlim= c(100,180), title= "Australian_map_with_")
  plot_NCP_on_world_map(NCP= "mean_Trophic_Level", ylim = c(-5, 30), xlim= c(-100,-55), title= "Caraib_map_with_")
  
  # #focus on public interest
  plot_NCP_on_world_map(NCP= "Public_Interest", ylim = c(-39, 0), xlim= c(100,180), title= "Australian_map_with_")
  plot_NCP_on_world_map(NCP= "Public_Interest", ylim = c(-5, 30), xlim= c(-100,-55), title= "Caraib_map_with_")
  
##-------------plot MPAs on map-------------
NCP_site$mpa_iucn_cat[which(NCP_site$mpa_iucn_cat == "Not Applicable" |
                              NCP_site$mpa_iucn_cat == "Not Reported")] <- NA
NCP_site <- NCP_site |>
  dplyr::bind_cols(
    protection = ifelse(NCP_site$mpa_enforcement == "High" &
                        stringr::str_detect(NCP_site$protection_status, "No take"),
                        "No take",
                        ifelse(is.na(NCP_site$mpa_name)==F, 
                               "Restricted", "Fished")) )
  # protection_med_high = ifelse(NCP_site$mpa_enforcement != "Low" &
  #                                stringr::str_detect(NCP_site$protection_status, "No take"),
  #                              "No take",
  #                              ifelse(is.na(NCP_site$mpa_name)==F, 
  #                                     "Restricted", "Fished")) )


plot_mpa <-function(NCP_site, xlim=c(-180,180), ylim = c(-36, 31)){
  ggplot(NCP_site) +
      geom_sf(data = coast, color = "grey30", fill = "lightgrey",
              aes(size=0.1)) +
      
      geom_point(size = 2, na.rm = T,
                 colour= "black",
                 alpha = 1,
                 aes(x = SiteLongitude, y = SiteLatitude,
                    shape = protection)) +
      
      coord_sf(xlim, ylim , expand = FALSE) +
      guides(alpha = "none", size = "none", colour = "none") +
      scale_shape_manual(values=c(19,17,18))+
    
      theme_minimal()+
      labs(title = paste0("MPA geographic distribution"),
           x="", y= "") +
      theme(legend.position = "right",
            plot.title = element_text(size=10, face="bold"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
      )
}

mpa <- plot_mpa(NCP_site , xlim=c(-180,180), ylim = c(-36, 31))
ggsave(filename = here::here("outputs", "figures", "RLS_sites_with_protection.jpg"),
       plot = mpa, width=15, height = 7 )

mpa <- mpa + 
  geom_rect(aes(xmin = 140, xmax = 160, ymin = -25, ymax = -10), color = "black", fill= "transparent")+
  geom_rect(aes(xmin = -95, xmax = -70, ymin = -5, ymax = 20), color = "black", fill= "transparent")

gold_coast_mpa<- plot_mpa( NCP_site,ylim = c(-25,-10), xlim= c(140,160))
caraib_mpa <- plot_mpa( NCP_site, ylim = c(-5, 20), xlim= c(-95,-70))

ggpubr::ggarrange(mpa, # First row with world map
                  ggpubr::ggarrange(caraib_mpa, gold_coast_mpa,  
                                    ncol = 2, labels = c("B", "C")), # Second row with zooms
                  nrow = 2, labels = "A") 
ggsave(plot = last_plot(), filename = here::here("outputs", "figures",
                                                 "RLS sites with protection and zoom.png"))
