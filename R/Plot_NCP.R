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
# load(here::here("outputs","NCP_site_coral_reef.Rdata"))
# load(here::here("outputs","NCP_site_SST20.Rdata"))
# load(here::here("outputs","NCP_site_coral_5_imputed.Rdata"))
# load(here::here("outputs","NCP_site_wo_australia.Rdata"))
# NCP_site <- NCP_site_condition

# Preping data  
## Clean data  
NCP_site_clean <- subset(NCP_site, select = -c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                               SiteMeanSST, SiteLatitude, SiteLongitude,
                                               HDI, MarineEcosystemDependency,
                                               coral_imputation, gravtot2, mpa_name,
                                               mpa_enforcement, protection_status, 
                                               mpa_iucn_cat))

## NCPs distribution with right skewed distribution
NCP_skewed_distribution <- c("Btot","recycling_N","recycling_P","Productivity",
                             "funct_distinctiveness","Omega_3_C","Calcium_C","Vitamin_A_C",
                             "phylo_entropy","ED_Mean", "iucn_species", "elasmobranch_diversity",
                             "low_mg_calcite", "high_mg_calcite", "aragonite", "monohydrocalcite",
                             "amorphous_carbonate", "biom_lowTL", "biom_mediumTL", "biom_highTL",
                             "fishery_biomass")

NCP_log_transformed <- NCP_site_clean |>
  dplyr::mutate(across(.cols = all_of(NCP_skewed_distribution),
                       .fns = ~ .x +1 , .names = "{.col}")) |>      # Adds 1 to values to log transformed
  dplyr::mutate(across(.cols = all_of(NCP_skewed_distribution),
                       .fns = log10 , .names = "{.col}")) |>       # log(x+1) to avoid negative values
  dplyr::rename_with(.cols = all_of(NCP_skewed_distribution),
                     .fn = ~ paste0("log(", .x, ")"))

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
plots <- lapply(colnames(NCP_site_clean), FUN = plot_distribution, data = NCP_site_clean )
all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + plots[[30]] + plots[[31]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here::here("outputs", "figures","NCP_distribution.png"), all_plot, width = 22, height =14 )


#### NCPs distribution with log correction
plots <- lapply(colnames(NCP_log_transformed), FUN = plot_distribution, data = NCP_log_transformed )

all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + plots[[30]] + plots[[31]] +
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
  x<- NCP_site$Btot
  x_title = "total Biomass"
  col <- fishualize::fish(n = ncol(NCP_site_clean), option = "Ostracion_whitleyi", begin = 0, end = 0.8)
  
  plots <- lapply( 1:ncol(NCP_site_clean), function(i){
    plot_correlation(x,colnames(NCP_site_clean)[i],i)
  })
  
  plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
    plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
    plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
    plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
    plots[[28]] + plots[[29]] + plots[[30]] + plots[[31]] +
    
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  ggsave(filename = here::here("outputs", "figures","NCP_correlation_with_biomass.png"), plot, width = 22, height =14 )
  
  
  
  ## With biodiversity
  x<- NCP_site$taxo_richness
  x_title = "taxonomic richness"
  col <- fishualize::fish(n = ncol(NCP_site_clean), option = "Ostracion_whitleyi", begin = 0, end = 0.8)
  
  plots <- lapply( 1:ncol(NCP_site_clean), function(i){
    plot_correlation(x,colnames(NCP_site_clean)[i],i)
  })
  
  plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
    plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
    plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
    plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
    plots[[28]] + plots[[29]] + plots[[30]] + plots[[31]] +
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  ggsave(filename = here::here("outputs", "figures","NCP_correlation_with_biodiversity.png"), plot, width = 22, height =14 )
  
  
  
  #### Corr-matrix for all NCPs
  png(filename = here::here("outputs", "figures","corr_matrix.png"), 
      width= 40, height = 30, units = "cm", res = 1000)
  print({
    M <- cor(NCP_site_clean)
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
    M <- cor(NCP_log_transformed)
    corrplot::corrplot(M, order = 'AOE', type = 'lower', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r')
    corrplot::corrplot(M, add = TRUE, type = 'upper', method = 'number', order = 'AOE', insig = 'p-value',
             diag = FALSE, tl.pos = 'n', cl.pos = 'n')
  })
  dev.off() 
  
##-------------Links of socio-envir with biomass -------------
  socio_envir <- c("HDI", "MarineEcosystemDependency", "coral_imputation",
                   "SiteLatitude", "mpa_iucn_cat")
  
  x <- as.list(NCP_log_transformed[,"log(Btot)"])[["log(Btot)"]]
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
  plot_NCP_on_world_map <- function(NCP = "taxo_richness", xlim=c(-180,180), ylim = c(-36, 31),
                                    title="world map with "){
    map <- ggplot(NCP_site) +
      geom_sf(data = coast, color = "grey30", fill = "lightgrey",
              aes(size=0.1)) +
      
      geom_point(data=NCP_site,
                 size = 4, shape = 20,
                 aes(x = SiteLongitude, y = SiteLatitude,
                     colour= NCP_site[,NCP][[1]],
                     alpha = scale(NCP_site[,NCP][[1]]))) +
      scale_colour_gradient(NCP,
                            low = "dodgerblue", high="darkred",
                            na.value=NA) +
      
      coord_sf(xlim, ylim, expand = FALSE) +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
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
  plot_NCP_on_world_map(NCP= "Btot", ylim = c(-39, 0), xlim= c(100,180), title= "Australian_map_with_")
  plot_NCP_on_world_map(NCP= "Btot", ylim = c(-5, 30), xlim= c(-100,-55), title= "Caraib_map_with_")
  
  #focus on robustness
  plot_NCP_on_world_map(NCP= "robustness", ylim = c(-39, 0), xlim= c(100,180), title= "Australian_map_with_")
  plot_NCP_on_world_map(NCP= "robustness", ylim = c(-5, 30), xlim= c(-100,-55), title= "Caraib_map_with_")
  
  #focus on mean TL
  plot_NCP_on_world_map(NCP= "mean_TL", ylim = c(-39, 0), xlim= c(100,180), title= "Australian_map_with_")
  plot_NCP_on_world_map(NCP= "mean_TL", ylim = c(-5, 30), xlim= c(-100,-55), title= "Caraib_map_with_")
  
  # #focus on Ntop
  # plot_NCP_on_world_map(NCP= "top_predator_proportion", ylim = c(-39, 0), xlim= c(100,180), title= "Australian_map_with_")
  # plot_NCP_on_world_map(NCP= "top_predator_proportion", ylim = c(-5, 30), xlim= c(-100,-55), title= "Caraib_map_with_")
  
##-------------plot MPAs on map-------------
NCP_site$mpa_iucn_cat[which(NCP_site$mpa_iucn_cat == "Not Applicable" |
                              NCP_site$mpa_iucn_cat == "Not Reported")] <- NA
NCP_site <- NCP_site |>
  dplyr::bind_cols(
    protection_not_low = ifelse(is.na(NCP_site$mpa_name)==F &
                            NCP_site$mpa_enforcement != "Low" &
                            stringr::str_detect(NCP_site$protection_status, "No take"),
                          "protected","not protected"),
      protection_only_high = ifelse(is.na(NCP_site$mpa_name)==F &
                                NCP_site$mpa_enforcement == "High" &
                                stringr::str_detect(NCP_site$protection_status, "No take"),
                              "protected","not protected"))


plot_mpa <-function(NCP_site){
  ggplot(NCP_site) +
      geom_sf(data = coast, color = "grey30", fill = "lightgrey",
              aes(size=0.1)) +
      
      geom_point(data=NCP_site,
                 size = 3, shape = 20,
                 aes(x = SiteLongitude, y = SiteLatitude,
                     colour= mpa_iucn_cat,
                     alpha = 0.3), na.rm = T) +
      
      coord_sf(xlim=c(-180,180), ylim = c(-36, 31), expand = FALSE) +
      guides(alpha = "none", size = "none") +
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

mpa <- plot_mpa(NCP_site = dplyr::filter(NCP_site, is.na(mpa_iucn_cat) == F))
ggsave(filename = here::here("outputs", "figures", "MPAs_iucn_cat_on_world_map.jpg"),
       plot = mpa, width=15, height = 7 )

mpa_not_low <- plot_mpa(NCP_site = dplyr::filter(NCP_site, protection_not_low == "protected"))
ggsave(filename = here::here("outputs", "figures", "MPAs_high_and_medium_enforcement_iucn_cat_on_world_map.jpg"),
       plot = mpa, width=15, height = 7 )

mpa_high <- plot_mpa(NCP_site = dplyr::filter(NCP_site, protection_only_high == "protected"))
ggsave(filename = here::here("outputs", "figures", "MPAs_high_enforcement_iucn_cat_on_world_map.jpg"),
       plot = mpa, width=15, height = 7 )
