################################################################################
##
## Display and record different ways to visualize Nature based contributions. 
##  All plots made in this script are saved in the file 'outputs/figures/'. 
##  All contributions are mapped and saved in 'outputs/figures/contributions_on_map'
##
## 1b_Plot_contributions.R
##
## 24/01/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
# pkgs <- c("here", "tidyverse", "ggplot2", "sf", "patchwork", "fishualize", "corrplot")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))
library(ggplot2)
library(patchwork)

rm(list=ls())

##-------------loading data-------------
load(here::here("data","metadata_surveys.Rdata"))
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

load(here::here("outputs","all_Contrib_site.Rdata"))
load(here::here("outputs","all_Contrib_site_log_transformed.Rdata"))

### To make test on other filtering condition, uncomment the wanted lines
  # load(here::here("outputs","Contrib_site_log_coral_reef.Rdata"))
  # load(here::here("outputs","Contrib_site_log_SST20.Rdata"))
  # load(here::here("outputs","Contrib_site_log_coral_5_imputed.Rdata"))
  # load(here::here("outputs","Contrib_site_log_wo_australia.Rdata"))
  # load(here::here("outputs","Contrib_site_log_random.Rdata"))
  # load(here::here("outputs","Contrib_site_log_only_australia.Rdata"))
  
  # Contrib_site_log_transformed <- Contrib_site_condition
###


# Preping data  
## Clean data  
Contrib_site_clean_before_log <- subset(Contrib_site, 
                                    select = -c(SiteCode, SiteCountry, SurveyDate,
                                                SiteEcoregion, SurveyDepth, 
                                                SiteMeanSST, SiteLatitude, SiteLongitude,
                                                HDI, MarineEcosystemDependency,
                                                coral_imputation, gravtot2, mpa_name,
                                                mpa_enforcement, protection_status, 
                                                mpa_iucn_cat))

Contrib_site_clean <- subset(Contrib_site_log_transformed, 
                         select = -c(SiteCode, SiteCountry, SurveyDate,
                                     SiteEcoregion, SurveyDepth, 
                                     SiteMeanSST, SiteLatitude, SiteLongitude,
                                     HDI, MarineEcosystemDependency,
                                     coral_imputation, gravtot2, mpa_name,
                                     mpa_enforcement, protection_status, 
                                     mpa_iucn_cat))


##-------------plot studied reefs and their protection status on map-------------
Contrib_site$mpa_iucn_cat[which(Contrib_site$mpa_iucn_cat == "Not Applicable" |
                                  Contrib_site$mpa_iucn_cat == "Not Reported")] <- NA
Contrib_site <- Contrib_site |>
  dplyr::bind_cols(
    protection = ifelse(Contrib_site$mpa_enforcement == "High" &
                          stringr::str_detect(Contrib_site$protection_status, "No take"),
                        "No take",
                        ifelse(is.na(Contrib_site$mpa_name)==F, 
                               "Restricted", "Fished")) )


plot_mpa <-function(Contrib_site, xlim=c(-180,180), ylim = c(-36, 31), legend_pos = "none"){
  ggplot(Contrib_site) +
    geom_sf(data = coast, color = "grey30", fill = "lightgrey",
            aes(size=0.1)) +
    
    geom_point(size = 2, na.rm = T,
               position = position_jitter(width =0.2, height =0.2),
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
    labs(title = "",
         x="", y= "") +
    theme(legend.position = legend_pos,
          plot.title = element_text(size=10, face="bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
}

mpa <-  plot_mpa(Contrib_site , xlim=c(-180,180), ylim = c(-36, 31), legend_pos = "none")+
  geom_rect(aes(xmin = 110, xmax = 160, ymin = -32, ymax = -7), color = "black", fill= "transparent")+
  geom_rect(aes(xmin = -95, xmax = -67, ymin = -3, ymax = 18), color = "black", fill= "transparent")

gold_coast_mpa<- plot_mpa( Contrib_site,ylim = c(-32, -7),
                           xlim= c(110,160), legend_pos = "right")
caraib_mpa <- plot_mpa( Contrib_site, ylim = c(-3, 18),
                        xlim= c(-95,-67), legend_pos = "none")

ggpubr::ggarrange(mpa, # First row with world map
                  ggpubr::ggarrange(caraib_mpa, gold_coast_mpa,  
                                    ncol = 2, labels = c("B", "C"), widths = c(1, 1.3)), # Second row with zooms
                  nrow = 2, labels = "A") 
ggsave(plot = last_plot(), filename = here::here("outputs", "figures",
                                                 "RLS sites with protection and zoom.png"),
       width = 13, height = 7)


##-------------Histograms of Contributions-------------
#### Contributions distribution
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

plots <- lapply(colnames(Contrib_site_clean_before_log), FUN = plot_distribution, data = Contrib_site_clean_before_log )
all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + plots[[30]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here::here("outputs", "figures","Contributions_raw_distribution.png"), all_plot, width = 22, height =14 )


#### Contributions distribution with log correction
plots <- lapply(colnames(Contrib_site_clean), FUN = plot_distribution, data = Contrib_site_clean )

all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + plots[[30]]+
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here::here("outputs", "figures","Contributions_log_transformed_distribution.png"), all_plot, width = 22, height =14 )



##-------------Correlations among Contributions-------------
plot_correlation <- function(x,y,i){  
  r_pearson <- cor(x, Contrib_site_log_transformed[,y][[y]])
  
  ggplot() +
    geom_point(aes(y = Contrib_site_log_transformed[,y][[y]], x = x),
               color = col[i], alpha = 0.2, size = 1) +
    theme_bw() +
    labs(x = x_title, y = y) +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 15))+
    annotate("text", x = -Inf, y = Inf, label = paste("r =", round(r_pearson, 2)), 
             hjust = -.3, vjust = 2, size = 5, color = "black")
}


## Relationship with the total biomass
x<- Contrib_site_log_transformed$Biomass
x_title = "Total Biomass"
# col <- fishualize::fish(n = ncol(Contrib_site_clean), option = "Ostracion_whitleyi", begin = 0, end = 0.8)
col <- c(rep("forestgreen", 20), rep("dodgerblue3", 10))
         
plots <- lapply( 2:ncol(Contrib_site_clean), function(i){
  plot_correlation(x,colnames(Contrib_site_clean)[i],i)
})

plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  plots[[28]] + plots[[29]] + #plots[[30]] +
  plot_layout(ncol=5)+
  
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here::here("outputs", "figures",
                             "Contrib_correlation_with_biomass.png"),
       plot, width = 20, height =18 )


corr_biom_tot <- cor(Contrib_site_clean)[,"Biomass"]
NN_cor_to_biom <- corr_biom_tot[
  c("N_Recycling","P_Recycling","Taxonomic_Richness","Functional_Entropy",
    "Phylogenetic_Entropy", "Trait_Distinctiveness","Evolutionary_Distinctiveness",
    "Herbivores_Biomass","Invertivores_Biomass","Piscivores_Biomass","Endemism",
    "Elasmobranch_Diversity","Low_Mg_Calcite","High_Mg_Calcite", "Aragonite",                   
    "Monohydrocalcite", "Amorphous_Carbonate", "Trophic_Web_Robustness", 
    "Mean_Trophic_Level" )]
summary(NN_cor_to_biom)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.1749  0.1342  0.4330  0.4282  0.6669  0.9736 

NP_cor_to_biom <- corr_biom_tot[
  c("Turnover_Available_Biomass", "Selenium","Zinc","Omega_3", "Calcium",
    "Iron", "Vitamin_A", "Available_Biomass", "Aesthetic", "Public_Attention")]
summary(NP_cor_to_biom)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.39988 -0.17313  0.02623  0.06492  0.22374  0.93905


#### Corr-matrix for all Contributions (log transformed)
M <- cor(dplyr::select(Contrib_site_clean, -Biomass))

png(filename = here::here("outputs", "figures","corr_matrix_log_transformed_Contributions.png"), 
    width= 40, height = 30, units = "cm", res = 1000)
print({
  corrplot::corrplot(M, order = 'AOE', type = 'lower', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r')
  corrplot::corrplot(M, add = TRUE, type = 'upper', method = 'number', order = 'AOE', insig = 'p-value',
                     diag = FALSE, tl.pos = 'n', cl.pos = 'n', number.digits = 1)
})
dev.off() 

## Extract p-value of this correlogram
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y)[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

p_val <- cor.test.p(as.data.frame(dplyr::select(Contrib_site_clean, -Biomass)))
length(which(p_val > 0.01))/2 #71 pairs are insignificantly correlated
summary(M[which(p_val > 0.01)]) # correlations lower than -0.07 and upper 0.07 are significant


## Study correlogram
pairwise_corr <- M[upper.tri(M)]
summary(pairwise_corr) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.70118 -0.08658  0.10082  0.12798  0.33006  0.93022 
summary(abs(pairwise_corr))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001535 0.097598 0.208872 0.259157 0.374540 0.930220 
pairwise_corr[order(pairwise_corr)]
length(which(pairwise_corr > 0.2)) #156
length(which(pairwise_corr < -0.2)) #53
length(which( pairwise_corr < 0.2 & -0.2 < pairwise_corr)) #197
length(which( pairwise_corr < 0.5 & -0.5 < pairwise_corr)) #356

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
ggsave(filename = here::here("outputs", "figures","hist_pearson_correlation_Contributions.png"), 
       histo_r_corr, width = 8, height =6 )


## Compare with random associations ##
load(here::here("outputs","Contrib_site_log_random.Rdata"))
Contrib_random <- subset(Contrib_site_condition, 
                             select = -c(SiteCode, SiteCountry, SurveyDate,
                                         SiteEcoregion, SurveyDepth, 
                                         SiteMeanSST, SiteLatitude, SiteLongitude,
                                         HDI, MarineEcosystemDependency,
                                         coral_imputation, gravtot2, mpa_name,
                                         mpa_enforcement, protection_status, 
                                         mpa_iucn_cat, Biomass))
M_random <- cor(Contrib_random)
pairwise_corr_random <- M_random[upper.tri(M_random)]


# Extract both histogramms
hist1 <- hist(pairwise_corr, breaks = 30, plot = FALSE)
df1 <- data.frame(mid = hist1$mids, density = hist1$density)

hist2 <- hist(pairwise_corr_random, plot = FALSE)
df2 <- data.frame(mid = hist2$mids, density = hist2$density)


histo_r_corr_random <- ggplot() +
  geom_bar(data = df1, aes(x = mid, y = density, fill = mid), 
           stat = "identity", color = "grey60", linewidth = 0.2,
           width = diff(hist1$breaks)[1]) +
  scale_fill_gradientn(colors = colorRampPalette(RdBu, bias = 1.3)(30))+
  
  labs(title = "Comparison between the correlation distributions of 
       calculated and randomized contributions", 
       x = "Pairwise correlation", fill = "Pearson \ncorrelation",
       y = "Density of contribution correlations") +
  ggnewscale::new_scale("fill") +
  
  geom_bar(data = df2,
           aes(x = mid, y = density / max(density) * (1.3*max(df1$density)),
               fill = "Random correlation"),
           stat = "identity", color = "grey60",linewidth = 0.2, alpha = 0.5,
           width = diff(hist2$breaks)[1]) +
  scale_fill_manual(values = c("Random correlation" = "grey60")) +
  labs(fill = "") +
  
  # Add the y axis for random correlations
  scale_y_continuous(sec.axis = sec_axis(~ . / max(df1$density) * max(df2$density),
                                         name = "Density of RANDOM contribution correlations")) +
  
  geom_vline(xintercept = 0.2, linetype = 3, col = "black") +
  geom_vline(xintercept = -0.2, linetype = 3, col = "black") +
  theme_bw() +
  theme(legend.position = "right", panel.grid.major = element_blank())

histo_r_corr_random
ggsave(filename = here::here("outputs", "figures",
       "hist_pearson_correlation_Contributions_with_random_associations.png"), 
       histo_r_corr_random, width = 8, height =6 )

##------------- Check spatial robustness with Mantel test ---------------
## all sites
corr_matrix_all_data <- cor(dplyr::select(Contrib_site_clean, -Biomass))

## all site, without Australia
load(here::here("outputs","Contrib_site_log_wo_australia.Rdata"))
Contrib_site_wo <- subset(Contrib_site_condition, 
                      select = -c(SiteCode, SiteCountry, SurveyDate,
                                  SiteEcoregion, SurveyDepth, 
                                  SiteMeanSST, SiteLatitude, SiteLongitude,
                                  HDI, MarineEcosystemDependency,
                                  coral_imputation, gravtot2, mpa_name,
                                  mpa_enforcement, protection_status, 
                                  mpa_iucn_cat,
                                  Biomass))
corr_matrix_wo_aust <- cor(Contrib_site_wo)

## australians sites only
load(here::here("outputs","Contrib_site_log_only_australia.Rdata"))
Contrib_site_only <- subset(Contrib_site_condition, 
                        select = -c(SiteCode, SiteCountry, SurveyDate,
                                    SiteEcoregion, SurveyDepth, 
                                    SiteMeanSST, SiteLatitude, SiteLongitude,
                                    HDI, MarineEcosystemDependency,
                                    coral_imputation, gravtot2, mpa_name,
                                    mpa_enforcement, protection_status, 
                                    mpa_iucn_cat,
                                    Biomass))
corr_matrix_only_aust <- cor(Contrib_site_only)

mantel_test_allVSaust <- vegan::mantel(corr_matrix_all_data, corr_matrix_only_aust)
mantel_test_aust <- vegan::mantel(corr_matrix_wo_aust, corr_matrix_only_aust)

## all sites xith SST > 20°C
load(here::here("outputs","Contrib_site_log_SST20.Rdata"))
Contrib_site_SST20 <- subset(Contrib_site_condition, 
                         select = -c(SiteCode, SiteCountry, SurveyDate,
                                     SiteEcoregion, SurveyDepth, 
                                     SiteMeanSST, SiteLatitude, SiteLongitude,
                                     HDI, MarineEcosystemDependency,
                                     coral_imputation, gravtot2, mpa_name,
                                     mpa_enforcement, protection_status, 
                                     mpa_iucn_cat,
                                     Biomass))
corr_matrix_SST20 <- cor(Contrib_site_SST20)
mantel_test_SST20 <- vegan::mantel(corr_matrix_all_data, corr_matrix_SST20) # r mantel stat = 0.91



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
  scale_colour_gradientn(name  ="Correlation with \n all data",
                         colours = RColorBrewer::brewer.pal(n = 8, name = "RdBu"))+
  xlab("Correlations without Australia")+
  ylab("Correlations in Australia only") +
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

ggsave(filename = here::here("outputs", "figures",
                             "Mantel_test_correlation_between_corrmatrix.png"),
       plot_mantel, width = 8, height =5)



##-------------plot Contributions on map-------------
#plot function
plot_Contrib_on_world_map <- function(NCP = "Taxonomic_Richness", xlim=c(-180,180), ylim = c(-36, 31),
                                  title="world map with ", jitter=1.5, pt_size=2){
  data <- Contrib_site_log_transformed[order(Contrib_site_log_transformed[,NCP][[1]]),]
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
  
  ggsave( here::here("outputs", "figures", "contributions_on_map", paste0( title , NCP, ".jpg")), plot = map, width=15, height = 7 )
  #map
}

# save world maps
parallel::mclapply(colnames(Contrib_site_clean), function(NCP){
  plot_Contrib_on_world_map(NCP, xlim=c(-180,180), ylim = c(-36, 31), 
                        title="world_map_with_", jitter=1.5, pt_size=2)
},mc.cores=parallel::detectCores()-5)

# save central america maps
parallel::mclapply(colnames(Contrib_site_clean), mc.cores=15, function(NCP){
  plot_Contrib_on_world_map(NCP, ylim = c(-5, 30), xlim= c(-100,-55),
                        title= "Caraib_map_with_", jitter=0.5, pt_size=4)
})

# save Pacific maps
parallel::mclapply(colnames(Contrib_site_clean), mc.cores=15, function(NCP){
  plot_Contrib_on_world_map(NCP,ylim = c(-27, -10), xlim= c(-180,-110),
                        title= "Pacific_map_with_", jitter=0.5, pt_size=4)
})