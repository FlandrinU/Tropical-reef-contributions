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
eez <- sf::st_read(here::here("data", "ShapeFiles coast", "eez_v11.shp"))

load(here::here("outputs","all_NCP_site.Rdata"))
load(here::here("outputs","all_NCP_site_log_transformed.Rdata"))
# load(here::here("outputs","NCP_site_log_coral_reef.Rdata"))
# load(here::here("outputs","NCP_site_log_SST20.Rdata"))
# load(here::here("outputs","NCP_site_log_coral_5_imputed.Rdata"))
# load(here::here("outputs","NCP_site_log_wo_australia.Rdata"))
# load(here::here("outputs","NCP_site_log_random.Rdata"))
# load(here::here("outputs","NCP_site_log_only_australia.Rdata"))
# NCP_site_log_transformed <- NCP_site_condition

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

library(ggplot2)

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
      geom_point(aes(y = NCP_site_log_transformed[,y][[y]], x = x),
                 color = col[i], alpha = 0.6, size = 1) +
      theme_bw() +
      labs(x = x_title, y = y) +
      theme(panel.grid.minor = element_blank(),
            axis.text = element_text(color = "black"), 
            axis.title = element_text(size = 17))
  }
  
  
  ## With Biomass
  x<- NCP_site_log_transformed$Biomass
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
  x<- NCP_site_log_transformed$Taxonomic_Richness
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
    M <- cor(dplyr::select(NCP_site_clean_before_log, -Biomass))
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
    M <- cor(dplyr::select(NCP_site_clean, -Biomass))
    corrplot::corrplot(M, order = 'AOE', type = 'lower', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r')
    corrplot::corrplot(M, add = TRUE, type = 'upper', method = 'number', order = 'AOE', insig = 'p-value',
             diag = FALSE, tl.pos = 'n', cl.pos = 'n', number.digits = 1)
  })
  dev.off() 
  
  # Study correlogramm
  pairwise_corr <- M[upper.tri(M)]
  summary(pairwise_corr)
      # Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
      # -0.69889 -0.08767  0.10600  0.12593  0.33824  0.97542 
  hist(pairwise_corr)
  length(which(pairwise_corr > 0.2)) #152
  length(which(pairwise_corr < -0.2)) #54
  length(which( pairwise_corr < 0.2 & -0.2 < pairwise_corr)) #200
  length(which( pairwise_corr < 0.5 & -0.5 < pairwise_corr)) #200
  
  RdBu = c("#67001F", "#B2182B", 
           "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", 
           "#92C5DE", "#4393C3", "#2166AC", "#053061")
  
  plot_histogram <- function(data = as.data.frame(pairwise_corr), 
                             x=pairwise_corr,
                             col= RdBu,
                             title=""){
    ggplot(data, aes(x)) +
      geom_histogram( bins = 30, color="grey60", aes(fill = after_stat(x)), 
                      linewidth=0.2)+
      scale_fill_gradientn(colors = colorRampPalette(col, bias = 1.3)(30))+
    
      labs(title = title, x = "Pairwise correlation", fill= "Pearson \ncorrelation")+
      geom_vline(xintercept = 0.2, linetype = 3, col= "black")+
      geom_vline(xintercept = -0.2, linetype = 3, col= "black")+
      theme_bw()+
      theme(legend.position = "right",    
            panel.grid.major = element_blank())
  } #end of plot_histogram function
  
  histo_r_corr <- plot_histogram(data = as.data.frame(pairwise_corr), 
                                 x=pairwise_corr,
                                 col= RdBu,
                                 title="")
  ggsave(filename = here::here("outputs", "figures","hist_pearson_correlation_NCP.png"), 
         histo_r_corr, width = 8, height =6 )
  
  
##------------- Check spatial robustness with Mantel test ---------------
  corr_matrix_all_data <- cor(dplyr::select(NCP_site_clean, -Biomass))
  
  load(here::here("outputs","NCP_site_log_wo_australia.Rdata"))
  NCP_site_wo <- subset(NCP_site_condition, 
                           select = -c(SiteCode, SiteCountry, SurveyDate,
                                       SiteEcoregion, SurveyDepth, 
                                       SiteMeanSST, SiteLatitude, SiteLongitude,
                                       HDI, MarineEcosystemDependency,
                                       coral_imputation, gravtot2, mpa_name,
                                       mpa_enforcement, protection_status, 
                                       mpa_iucn_cat,
                                       Biomass))
  corr_matrix_wo_aust <- cor(NCP_site_wo)
  
  load(here::here("outputs","NCP_site_log_only_australia.Rdata"))
  NCP_site_only <- subset(NCP_site_condition, 
                           select = -c(SiteCode, SiteCountry, SurveyDate,
                                       SiteEcoregion, SurveyDepth, 
                                       SiteMeanSST, SiteLatitude, SiteLongitude,
                                       HDI, MarineEcosystemDependency,
                                       coral_imputation, gravtot2, mpa_name,
                                       mpa_enforcement, protection_status, 
                                       mpa_iucn_cat,
                                       Biomass))
  corr_matrix_only_aust <- cor(NCP_site_only)
  
  mantel_test_allVSaust <- vegan::mantel(corr_matrix_all_data, corr_matrix_only_aust)
  mantel_test_aust <- vegan::mantel(corr_matrix_wo_aust, corr_matrix_only_aust)
  
  
  load(here::here("outputs","NCP_site_log_SST20.Rdata"))
  NCP_site_SST20 <- subset(NCP_site_condition, 
                          select = -c(SiteCode, SiteCountry, SurveyDate,
                                      SiteEcoregion, SurveyDepth, 
                                      SiteMeanSST, SiteLatitude, SiteLongitude,
                                      HDI, MarineEcosystemDependency,
                                      coral_imputation, gravtot2, mpa_name,
                                      mpa_enforcement, protection_status, 
                                      mpa_iucn_cat,
                                      Biomass))
  corr_matrix_SST20 <- cor(NCP_site_SST20)
  mantel_test_SST20 <- vegan::mantel(corr_matrix_all_data, corr_matrix_SST20)
  
  
  
  ## Plot correlation between corr_matrix
  corr_matrix_wo_aust[upper.tri(corr_matrix_wo_aust)] <- NA
  corr_matrix_only_aust[upper.tri(corr_matrix_only_aust)] <- NA
  corr_matrix_all_data[upper.tri(corr_matrix_all_data)] <- NA
  
  mantel_data <- dplyr::left_join(reshape2::melt(corr_matrix_wo_aust),
                                  reshape2::melt(corr_matrix_only_aust),
                                  by=c("Var1", "Var2")) |> 
    dplyr::left_join(reshape2::melt(corr_matrix_all_data)) |> 
    dplyr::rename(wo_aust= "value.x", only_aust = "value.y", all_data="value")|>
    na.omit()
    
  
  plot_mantel <- ggplot(mantel_data, aes(x=wo_aust, y=only_aust, colour = all_data)) +
    geom_point()+theme_bw()+
    scale_colour_gradientn(name  ="correlation with \n all data",
                           colours = RColorBrewer::brewer.pal(n = 8, name = "RdBu"))+
    xlab("Correlations without australia")+
    ylab("Correlations in australia only") +
    guides(colour = guide_colourbar(title.position="top")) +
    theme(legend.position = "right",
          legend.key.size = unit(0.5, 'cm'),
          legend.direction = "vertical",
          legend.title = element_text( size = 11),
          legend.background = element_rect(fill='transparent'),
          axis.title=element_text(size=14)) +
    geom_abline(intercept = 0, slope = 1,color="#757575",linetype = "dashed",
                linewidth = 1)+
    geom_label(aes(label = paste0("Mantel statistic r: ", 
                                  round(mantel_test_aust[["statistic"]], 3),
                                  ", p = ", mantel_test_aust[["signif"]], 
                                  "\n 95% quantiles in null model: ", 
                                  round(quantile(mantel_test_aust[["perm"]], 0.95),2)),
                   y = 0.9, x = -0.8), size = 3,
               color = "black", hjust = 0)+
    ggrepel::geom_label_repel(
      data=dplyr::filter(mantel_data, abs(wo_aust - only_aust) >
                           quantile(abs(wo_aust - only_aust), 0.98) ),
      aes(label= paste(Var1, "-", Var2)),
      size=2, fill = "white", 
      min.segment.length = 0.1,
      color = "black", alpha = 0.8,
      direction = "both",
      seed = 1968)
  
  plot_mantel
  ggsave(filename = here::here("outputs", "figures",
                               "Mantel_test_correlation_between_corrmatrix.png"),
         plot_mantel, width = 8, height =5)
  
  
##-------------Links of socio-envir with biomass -------------
  socio_envir <- c("HDI", "MarineEcosystemDependency", "coral_imputation",
                   "SiteLatitude", "mpa_iucn_cat", "gravtot2")
  
  x <- as.list(NCP_site_clean[,"Biomass"])[["Biomass"]]
  x_title = "log(total Biomass)"
  col <- fishualize::fish(n = length(socio_envir), option = "Centropyge_loricula", begin = 0, end = 0.8)
  
  plots <- lapply( 1:length(socio_envir), function(i){
    plot_correlation(x,socio_envir[i],i) +
    coord_flip() 
  })
  
  plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]]+ plots[[6]]+
    theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = 'bold'))
  
  ggsave(filename = here::here("outputs", "figures","Socio_envir_variables_with_biomass.png"), plot, width = 22, height =14 )
  
  
  
##-------------plot NCPs on map-------------
  #plot function
  plot_NCP_on_world_map <- function(NCP = "Taxonomic_Richness", xlim=c(-180,180), ylim = c(-36, 31),
                                    title="world map with ", jitter=1.5, pt_size=2){
    data <- NCP_site_log_transformed[order(NCP_site_log_transformed[,NCP][[1]]),]
    library(ggplot2)
    map <- ggplot(data) +
      geom_sf(data = coast, color = "grey30", fill = "lightgrey",
              size=0.1) +
      
      geom_point(data=data,
                 size = pt_size, shape = 20,
                 position=position_jitter(width=jitter, height = jitter),
                 aes(x = SiteLongitude, y = SiteLatitude,
                     colour= data[,NCP][[1]],
                     alpha = 0.7)) +
      scale_colour_gradientn(name  = NCP,
                             colours = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu")))+
      
      
      coord_sf(xlim, ylim, expand = FALSE) +
      guides(alpha= "none" ) +
      # scale_size_continuous(range = c(0.5, 4), guide = "none") +
      theme_minimal()+
      labs(#title = paste0(NCP, " geographic distribution"),
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
  
  # save world maps
parallel::mclapply(colnames(NCP_site_clean), mc.cores=15, function(NCP){
    plot_NCP_on_world_map(NCP, xlim=c(-180,180), ylim = c(-36, 31), 
                          title="world_map_with_", jitter=1.5, pt_size=2)
})
  
# save central america maps
parallel::mclapply(colnames(NCP_site_clean), mc.cores=15, function(NCP){
  plot_NCP_on_world_map(NCP, ylim = c(-5, 30), xlim= c(-100,-55),
                        title= "Caraib_map_with_", jitter=0.5, pt_size=4)
})

# save Pacific maps
parallel::mclapply(colnames(NCP_site_clean), mc.cores=15, function(NCP){
  plot_NCP_on_world_map(NCP,ylim = c(-27, -10), xlim= c(-180,-110),
                        title= "Pacific_map_with_", jitter=0.5, pt_size=4)
})

  #focus on Biomass
  plot_NCP_on_world_map(NCP= "Biomass", ylim = c(-39, 0), xlim= c(100,180),
                        title= "Australian_map_with_", jitter=0.5, pt_size=4)
  plot_NCP_on_world_map(NCP= "Biomass", ylim = c(-5, 30), xlim= c(-100,-55),
                        title= "Caraib_map_with_", jitter=0.5, pt_size=4)
  plot_NCP_on_world_map(NCP= "Biomass", ylim = c(-27, -10), xlim= c(-180,-110),
                        title= "Pacific_map_with_", jitter=0.5, pt_size=4)
  
  #focus on robustness
  plot_NCP_on_world_map(NCP= "Trophic_web_robustness", ylim = c(-39, 0), xlim= c(100,180),
                        title= "Australian_map_with_", jitter=0.5, pt_size=4)
  plot_NCP_on_world_map(NCP= "Trophic_web_robustness", ylim = c(-5, 30), xlim= c(-100,-55),
                        title= "Caraib_map_with_", jitter=0.5, pt_size=4)
  
  #focus on mean TL
  plot_NCP_on_world_map(NCP= "mean_Trophic_Level", ylim = c(-39, 0), xlim= c(100,180),
                        title= "Australian_map_with_", jitter=0.5, pt_size=4)
  plot_NCP_on_world_map(NCP= "mean_Trophic_Level", ylim = c(-5, 30), xlim= c(-100,-55),
                        title= "Caraib_map_with_", jitter=0.5, pt_size=4)
  
  # #focus on public interest
  plot_NCP_on_world_map(NCP= "Public_Interest", ylim = c(-39, 0), xlim= c(100,180),
                        title= "Australian_map_with_", jitter=0.5, pt_size=4)
  plot_NCP_on_world_map(NCP= "Public_Interest", ylim = c(-5, 30), xlim= c(-100,-55),
                        title= "Caraib_map_with_", jitter=0.5, pt_size=4)
  
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
                 alpha = 1,
                 colour = "black",
                 stroke=0.1,
                 aes(x = SiteLongitude, y = SiteLatitude,
                    shape = protection,
                    fill=protection)) +
      
      coord_sf(xlim, ylim , expand = FALSE) +
      guides(alpha = "none", size = "none", colour = "none") +
      scale_shape_manual(values=c(21,24,23))+
    
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
                                                 "RLS sites with protection and zoom.png"),
       width = 13, height = 7)




##-------------Number of sites per EEZ-------------
sf_use_s2(FALSE)
sites_rls <- sf::st_as_sf(NN_NS_with_product,
                          coords=c("SiteLongitude", "SiteLatitude"),
                          crs=4326)
eez$count_pt <- lengths(sf::st_intersects(eez, sites_rls))

sampling_effort <- ggplot() +
  geom_sf(data = eez, color = "grey30", aes(fill=count_pt) ) +
  scale_fill_gradient(name = "count", trans = "log10",
                      na.value = "grey70", breaks = c(1,10,100,500)) +
  geom_sf(data = coast, color = NA, fill = "grey90") +
  coord_sf(ylim= c(-45,45),expand = FALSE) +
  theme_bw()+
  labs(x="Longitude", y= "Latitude", title = "Number of site by EEZ") +
  theme(axis.title.x = element_text(face = "bold",
                                    size = 13),
        axis.title.y = element_text(face = "bold",
                                    size = 13),
        axis.text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        plot.title = element_text(size=15, face="bold")
  )

# sampling_effort
ggsave(plot = sampling_effort, filename = here::here("outputs", "figures",
                                                 "Number of RLS sites per EEZ.png"),
       width = 13, height = 5)
