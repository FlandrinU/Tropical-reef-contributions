################################################################################
##
## Makes Reduction dimensions analyses on all sites, with all NCPs. 
##   /!\ need python and miniconnda on the computer
##    then download tensorflow with tensorflow::install_tensorflow()
##    and UMAP module on python (ubuntu: $ pip install umap-learn)
##
## dimensions_reduction_on_NCP.R
##
## 09/01/2023
##
## Ulysse Flandrin
##
################################################################################

#----------------- Loading packages -------------------
pkgs <- c("here", "tidyverse", "tibble", "questionr", "dplyr", "Rtsne",
          "dimRed", "tensorflow", "ggplot2", "scico")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##------------- loading data -------------
load(here::here("outputs","all_NCP_site_log_transformed.Rdata"))
load(here::here("data","metadata_surveys.Rdata"))

##------------- preping data -------------
#NCP_site <- dplyr::filter(NCP_site, SiteCountry != "Australia")
NCP_site_clean <- subset(NCP_site_log_transformed, 
                         select = -c(Biomass, SiteCode, SurveyDate,
                                     SiteCountry, SiteEcoregion, SurveyDepth, 
                                     SiteMeanSST, SiteLatitude, SiteLongitude,
                                     HDI, MarineEcosystemDependency,
                                     coral_imputation, gravtot2, mpa_name,
                                     mpa_enforcement, protection_status, 
                                     mpa_iucn_cat))
NCP_scaled <- scale(NCP_site_clean)


##------------- all analysis -------------

embed_methods <- c("PCA", "tSNE" , "UMAP", "AutoEncoder")
## apply dimensionality reduction 
data_emb <- lapply(embed_methods, function(x){ 
  dimRed::embed(NCP_scaled,x, ndim=2)}) 
names(data_emb) <- embed_methods

lapply(data_emb, plot, type = "2vars")
dimRed::plot_R_NX(data_emb)


embed_methods <- c("PCA", "tSNE" , "UMAP")
## apply dimensionality reduction 
data_emb <- lapply(embed_methods, function(x){ 
  dimRed::embed(NCP_scaled, x, ndim = 3 )}) 
names(data_emb) <- embed_methods

lapply(data_emb, plot, type = "3vars")
dimRed::plot_R_NX(data_emb)

# Tuto pkg
# embed_methods <- c("Isomap", "PCA") 
# data_set <- loadDataSet("3D S Curve", n = 1000)
# data_emb <- lapply(embed_methods, function(x) embed(data_set, x)) 
# names(data_emb) <- embed_methods 
# plot(data_set, type = "3vars")
# lapply(data_emb, plot, type = "2vars")
# plot_R_NX(data_emb)

##------------- See analysis separately -------------
## PCA
pca <- dimRed::embed( NCP_scaled, "PCA", ndim = 3)
plot(pca, type = "3vars")
dimRed::AUC_lnK_R_NX(pca)


## encoding

# /!\ need python and miniconnda on the computer
# then download tensorflow with tensorflow::install_tensorflow()
encoding <- dimRed::embed( NCP_scaled,  "AutoEncoder", ndim=3)
plot(encoding, type = "2vars")
dimRed::AUC_lnK_R_NX(encoding)


## UMAP
umap <- dimRed::embed( NCP_scaled,  "UMAP", ndim=3)
plot(umap, type = "2vars")
dimRed::AUC_lnK_R_NX(umap)
umap_coordinates <- cbind( umap@data@data , 
                           subset(NCP_site, select = c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                                       SiteMeanSST, SiteLatitude, SiteLongitude, Btot,
                                                       HDI, gravtot2, MarineEcosystemDependency,
                                                       coral_imputation)))

cluster <- ggplot(umap_coordinates) +
  aes(x = UMAP1, y = UMAP2, colour = SiteCountry) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_hue(direction = 1) +
  theme_minimal()
ggsave(filename = here::here("outputs", "figures", "umap_countries_distribution.png"), width=17, height= 10)

##-------------focus on t-SNE ------------- 
tsne <- dimRed::embed( NCP_scaled, "tSNE", ndim = 3)
plot(tsne, type = "3vars")
dimRed::AUC_lnK_R_NX(tsne)
tsne_coordinates <- cbind( tsne@data@data , 
                           subset(NCP_site_log_transformed, select = c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                                                         SiteMeanSST, SiteLatitude, SiteLongitude,
                                                                         HDI, gravtot2, MarineEcosystemDependency,
                                                                         coral_imputation)))

# NCP_scaled_tsne <- Rtsne::normalize_input( as.matrix(NCP_unique) )
# tsne <- Rtsne::Rtsne(NCP_scaled_tsne, dims = 3)
# plot(tsne$Y[,c(1,2)], col=NCP_site$SiteMeanSST, asp=1)
# dimRed::AUC_lnK_R_NX(tsne)

plotly::plot_ly(tsne_coordinates,
                x= ~tSNE1, y= ~tSNE2, #z= ~tSNE3,
                size = 10 ) |> 
  plotly::add_markers(color= ~SiteCountry)


#temperature pattern
library(ggplot2)
ggplot(data =tsne_coordinates) +
 aes(x = tSNE1, y = tSNE2, 
     colour = SiteMeanSST) +
 geom_point(shape = "circle", 
            size = 1.5) +
 geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.3)+
 geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.3)+
 scale_color_distiller(palette = "RdYlBu", direction = -1) +
 labs(title = "Sites distribution in t-SNE dimension") +
 theme_minimal()

ggplot(tsne_coordinates, aes(tSNE2 , SiteMeanSST , alpha = 0.5)) +
  geom_point()+
  geom_smooth() +
  labs(x = "Mean SST", y = "Dimension  in global PCA", title = "Importance of SST in tSNE's axe 2") +
  theme_minimal()+
  theme(legend.position = 'none')

# MED pattern
ggplot(data =tsne_coordinates) +
  aes(x = tSNE1, y = tSNE2, 
      colour = MarineEcosystemDependency) +
  geom_point(shape = "circle", 
             size = 1.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.3)+
  scale_color_distiller(palette = "RdYlBu", direction = -1) +
  labs(title = "Sites distribution in t-SNE dimension") +
  theme_minimal()

#gravtot2 pattern
ggplot(data =tsne_coordinates) +
  aes(x = tSNE1, y = tSNE2, 
      colour = log(gravtot2)) +
  geom_point(shape = "circle", 
             size = 1.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.3)+
  scale_color_distiller(palette = "RdYlBu", direction = -1) +
  labs(title = "Sites distribution in t-SNE dimension") +
  theme_minimal()

#countries distribution
tsne_coordinates |> 
  dplyr::filter(SiteCountry != '') |>
  ggplot() +
  aes(x = tSNE1, y = tSNE2, 
      colour = SiteCountry) +
  geom_point(shape = "circle", 
             size = 1.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.3)+
  scale_color_hue() +
  labs(title = "Countries distribution in t-SNE dimension") +
  theme_minimal()

##-------------Study NN and NS separately with t-SNE ------------- 
grp <- as.factor( c(recycling_N="NN", recycling_P="NN",Productivity="NS",taxo_richness="NN", funct_entropy="NN",
         funct_distinctiveness="NN", Selenium_C="NS", Zinc_C="NS", Omega_3_C="NS", Calcium_C="NS",
         Iron_C="NS", Vitamin_A_C="NS", phylo_entropy="NN", ED_Mean="NN", aesthe_survey="NS", iucn_species="NN",
         elasmobranch_diversity="NN", low_mg_calcite="NN", high_mg_calcite="NN", aragonite="NN",
         monohydrocalcite="NN", amorphous_carbonate="NN", biom_lowTL="NN", biom_mediumTL="NN",
         biom_highTL="NN", fishery_biomass="NS") )

### NN
NN <- names(grp)[ grp=="NN" ]
NCP_NN <- NCP_scaled[,NN]


# embed_methods <- c("PCA", "tSNE" , "UMAP", "AutoEncoder")
# ## apply dimensionality reduction 
# data_emb <- lapply(embed_methods, function(x){ 
#   dimRed::embed(NCP_NN,x, ndim=2)}) 
# names(data_emb) <- embed_methods
# 
# lapply(data_emb, plot, type = "2vars")
# plot_R_NX(data_emb)

tsne_NN <- dimRed::embed( NCP_NN, "tSNE", ndim = 2)
plot(tsne_NN, type = "2vars")
dimRed::AUC_lnK_R_NX(tsne_NN)
tsne_NN_coordinates <- cbind( tsne_NN@data@data , 
                           subset(NCP_site, select = c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                                       SiteMeanSST, SiteLatitude, SiteLongitude, Btot,
                                                       HDI, gravtot2, MarineEcosystemDependency,
                                                       coral_imputation)))


#countries distribution
tsne_NN_coordinates |> 
  filter(SiteCountry != '') |>
  ggplot() +
  aes(x = tSNE1, y = tSNE2, 
      colour = SiteCountry) +
  geom_point(shape = "circle", 
             size = 1.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.3)+
  scale_color_hue() +
  labs(title = "Countries distribution in t-SNE dimension") +
  theme_minimal()

# NN score
ggplot(data =tsne_NN_coordinates) +
  aes(x = tSNE1, y = tSNE2, 
      colour = tSNE1+tSNE2) +
  geom_point(shape = "circle", 
             size = 1.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.3)+
  scale_colour_gradient2(low = "black", mid = "white", high = "forestgreen", guide = "colourbar")+
  labs(title = "Sites distribution in t-SNE dimension") +
  theme_minimal()

# plot on map
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

map_NN <- ggplot(tsne_NN_coordinates) +
  geom_sf(data = coast, color = NA, fill = "lightgrey") +
  geom_point(aes(x = SiteLongitude, y = SiteLatitude,
                 color = tSNE1+tSNE2, alpha= abs(tSNE1+tSNE2),
                 size = 2)) +
  # geom_point(aes(x = SiteLongitude, y = SiteLatitude),
  #            shape = 1, size = 2, stroke = 0.5,
  #            color= "black",
  #            data = head(NN_NS_with_col[order(NN_NS_with_col$product_d_r, decreasing = T),], 15)) +
  
  scale_colour_gradientn(colours = colorRampPalette(rev(c( "forestgreen","white", "grey30")))(1000)) +
  scale_alpha_continuous(range = c(0, 1)) +
  
  coord_sf(ylim = c(-36, 31), expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
  theme_minimal()+
  labs(title = "Nature for Nature",
       x="Longitude", y= "Latitude") +
  theme(
    plot.title = element_text(size=10, face="bold"),
    legend.position = 'none'
  )
map_NN






### NS
NS <- names(grp)[ grp=="NS" ]
NCP_NS <- NCP_scaled[,NS]


# embed_methods <- c("PCA", "tSNE" , "UMAP", "AutoEncoder")
# ## apply dimensionality reduction 
# data_emb <- lapply(embed_methods, function(x){ 
#   dimRed::embed(NCP_NS,x, ndim=2)}) 
# names(data_emb) <- embed_methods
# 
# lapply(data_emb, plot, type = "2vars")
# plot_R_NX(data_emb)

tsne_NS <- dimRed::embed( NCP_NS, "tSNE", ndim = 2)
plot(tsne_NS, type = "2vars")
dimRed::AUC_lnK_R_NX(tsne_NS)
tsne_NS_coordinates <- cbind( tsne_NS@data@data , 
                              subset(NCP_site, select = c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                                          SiteMeanSST, SiteLatitude, SiteLongitude, Btot,
                                                          HDI, gravtot2, MarineEcosystemDependency,
                                                          coral_imputation)))


#countries distribution
tsne_NS_coordinates |> 
  filter(SiteCountry != '') |>
  ggplot() +
  aes(x = tSNE1, y = tSNE2, 
      colour = SiteCountry) +
  geom_point(shape = "circle", 
             size = 1.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.3)+
  scale_color_hue() +
  labs(title = "Countries distribution in t-SNE dimension") +
  theme_minimal()

# NS score
ggplot(data =tsne_NS_coordinates) +
  aes(x = tSNE1, y = tSNE2, 
      colour = tSNE1+tSNE2) +
  geom_point(shape = "circle", 
             size = 1.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size = 0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.3)+
  scale_colour_gradient2(low = "black", mid = "white", high = "dodgerblue3", guide = "colourbar")+
  labs(title = "Sites distribution in t-SNE dimension") +
  theme_minimal()

# plot on map
map_NS <- ggplot(tsne_NS_coordinates) +
  geom_sf(data = coast, color = NA, fill = "lightgrey") +
  geom_point(aes(x = SiteLongitude, y = SiteLatitude,
                 color = tSNE1+tSNE2, alpha= abs(tSNE1+tSNE2),
                 size = 2)) +
  # geom_point(aes(x = SiteLongitude, y = SiteLatitude),
  #            shape = 1, size = 2, stroke = 0.5,
  #            color= "black",
  #            data = head(NS_NS_with_col[order(NS_NS_with_col$product_d_r, decreasing = T),], 15)) +
  
  scale_colour_gradientn(colours = colorRampPalette(rev(c( "dodgerblue3","white", "grey30")))(1000)) +
  scale_alpha_continuous(range = c(0, 1)) +
  
  coord_sf(ylim = c(-36, 31), expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
  theme_minimal()+
  labs(title = "Nature for Society",
       x="Longitude", y= "Latitude") +
  theme(
    plot.title = element_text(size=10, face="bold"),
    legend.position = 'none'
  )
map_NS
