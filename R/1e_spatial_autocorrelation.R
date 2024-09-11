################################################################################
##
## This script make the spatial autocorrelation analyse of the NN and NP scores
##
## 1e_spatial_autocorrelation.R
##
## 09/03/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
# pkgs <- c("spdep", "geodist", "ape", "ncf", "nlme", "ggplot2")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)

rm(list=ls())

##------------- loading data -------------
load( file = here::here("outputs", "NN_NP_score_wheighted_mean.Rdata"))
load(here::here("outputs","all_Contrib_site_log_transformed.Rdata"))

##------------- preping data -------------
all_data <- dplyr::left_join(Contrib_site_log_transformed, NN_NP_scores)

unique_points <- all_data[which(duplicated(paste0(all_data$SiteLongitude, all_data$SiteLatitude)) == F),]
coord <- dplyr::select(unique_points, longitude = SiteLongitude, latitude = SiteLatitude)

# Calculate the spatial distance between reefs
site_dist <- as.matrix( geodist::geodist(coord[,c("longitude", "latitude")]
                                         , measure = "geodesic"))
site_dist_inv <- 1/site_dist
diag(site_dist_inv) <- 0


## -------------Assess the Moran index per distance class -------------
correlog_plot <- function(score = "NN_score",
                          title = "Nature for Nature scores",
                          increment = 200, 
                          resamp= 99){
  cat("compute correlog ... \n")
  spatial_cor <- ncf::correlog(x= all_data$SiteLongitude, 
                               y= all_data$SiteLatitude, 
                               z = all_data[,score][[score]], increment = increment, 
                               resamp= resamp, latlon = TRUE) #distance and increment are in km
  library(ggplot2)
  ggplot() +
    geom_point(aes(y = spatial_cor$correlation[ which(spatial_cor$mean.of.class < 15000) ], 
                   x = spatial_cor$mean.of.class[ which(spatial_cor$mean.of.class < 15000) ]),
               alpha = 0.6, size = 1) +
    geom_smooth(aes(y = spatial_cor$correlation[ which(spatial_cor$mean.of.class < 15000) ], 
                    x = spatial_cor$mean.of.class[ which(spatial_cor$mean.of.class < 15000) ]))+
    
    geom_point(aes(y = spatial_cor$correlation[ which(spatial_cor$p <= 0.01 &
                                                        spatial_cor$mean.of.class < 15000) ] ,
                   x = spatial_cor$mean.of.class[ which(spatial_cor$p <= 0.01 &
                                                          spatial_cor$mean.of.class < 15000) ]),
               alpha = 0.8, size = 1, col ="red") +
    geom_vline(xintercept = spatial_cor[["x.intercept"]], linetype="dashed", 
               color = "black", linewidth=1)+
    
    xlab("Distance class (in km)") + ylab("spatial correlation (Moran I)")+
    labs(title =  title,
         subtitle = paste("Increment = ", increment, "km, ", "permutations = ", resamp,
                          ", x intercept = ", round(spatial_cor[["x.intercept"]],1), "km"))+
    # ylim(-0.7,0.9)+
    theme_bw()
}

#NN score
spatial_cor_NN_100 <- correlog_plot(score = "NN_score",title = "Nature for Nature", increment = 100, resamp= 99)
spatial_cor_NN_100

#NP score
spatial_cor_NP_100 <- correlog_plot(score = "NP_score",title = "Nature for People", increment = 100, resamp= 99)
spatial_cor_NP_100

#panel
ggpubr::ggarrange(spatial_cor_NN_100, spatial_cor_NP_100,
                  ncol = 2, labels = c("A", "B")) 
ggsave(plot = last_plot(), width=12, height =6,
       filename = here::here("outputs", "figures", "Spatial correlation per distance _ PANEL NN and NP.png"))



##-------------Variogram-------------
dist_matrix <- geodist::geodist(x= coord[,c("longitude", "latitude")],
                                measure = "geodesic") #distance in meters
range(dist_matrix) #max ~20 000 km => OK

site_dist <- dist_matrix[lower.tri(dist_matrix, diag = FALSE)]
vario <- nlme::Variogram(unique_points$NN_score, distance = site_dist)

variogram <- ggplot(data= vario) +
  geom_point(aes(y = variog, x = dist),
             col = "grey", alpha = 0.2, size = 1) +
  geom_smooth(aes(y = variog, x = dist))+
  theme_bw()
ggsave(plot = variogram, filename = here::here("outputs", "figures",
                                               "Spatial correlation variogram.png"))

