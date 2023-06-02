################################################################################
##
## Makes spatial autocorrelation analyses on NN and NS scores
##
## spatial_autocorrelation.R
##
## 09/03/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("spdep", "geodist", "ape", "ncf", "nlme", "ggplot2")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load( file = here::here("outputs", "NN_NS_score_wheighted_mean.Rdata"))
load(here::here("outputs","all_NCP_site_log_transformed.Rdata"))

##-------------preping data-------------
all_data <- dplyr::left_join(NCP_site_log_transformed, NN_NS_scores)

unique_points <- all_data[which(duplicated(paste0(all_data$SiteLongitude, all_data$SiteLatitude)) == F),]
coord <- dplyr::select(unique_points, longitude = SiteLongitude, latitude = SiteLatitude)

site_dist <- as.matrix( geodist::geodist(coord[,c("longitude", "latitude")]
                                         , measure = "geodesic"))
site_dist_inv <- 1/site_dist
diag(site_dist_inv) <- 0

## -------------Moran index-------------
ape::Moran.I(unique_points$NN_score, site_dist_inv)

listw <- spdep::mat2listw(site_dist_inv, style = "M")
spdep::moran.test(x= unique_points$NN_score, listw = listw, randomisation = T )

#check Spatial autocorrelation 
coords<-as.matrix(cbind(NN_NS_scores[,"SiteLongitude"],NN_NS_scores[,"SiteLatitude"])) 
nb <- spdep::knn2nb(spdep::knearneigh(coords, 1, longlat = TRUE)) 
lstw <- spdep::nb2listw((spdep::knn2nb(spdep::knearneigh(coords, k=1, longlat = TRUE))))
spdep::moran.test(NN_NS_scores$NN_score, lstw) # there is SAC for model residuals

#check Spatial autocorrelation 
lstw <- spdep::nb2listw((spdep::knn2nb(spdep::knearneigh(coords, k=2, longlat = TRUE))))
spdep::moran.test(NN_NS_scores$NN_score, lstw) # there is SAC 
spdep::moran.mc(NN_NS_scores$NN_score, lstw,nsim=999)
spdep::moran.plot(NN_NS_scores$NN_score, lstw)

## -------------Moran index per distance class-------------
correlog_plot <- function(score = "NN_score",
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
    labs(title = paste("Spatial correlation of", score),
         subtitle = paste("Increment = ", increment, "km, ", "resample = ", resamp,
                          ", x intercept = ", round(spatial_cor[["x.intercept"]],1), "km"))+
    # ylim(-0.7,0.9)+
    theme_bw()
}

#NN score
spatial_cor_NN <- correlog_plot(score = "NN_score", increment = 500, resamp= 30)
spatial_cor_NN
ggsave(plot = spatial_cor_NN, filename = here::here("outputs", "figures",
        "Spatial correlation per distance NN, increment 500.png"))

spatial_cor_NN <- correlog_plot(score = "NN_score", increment = 100, resamp= 99)
ggsave(plot = spatial_cor_NN, filename = here::here("outputs", "figures",
        "Spatial correlation per distance NN, increment 100.png"))

spatial_cor_NN <- correlog_plot(score = "NN_score", increment = 200, resamp= 99)
ggsave(plot = spatial_cor_NN, width=6, height =6,
       filename = here::here("outputs", "figures",
        "Spatial correlation per distance NN, increment 200.png"))


spatial_cor_NN <- correlog_plot(score = "NN_score", increment = 5000, resamp= 30)
ggsave(plot = spatial_cor_NN, filename = here::here("outputs", "figures",
         "Spatial correlation per distance NN, increment 5000.png"))

#NS score
spatial_cor_NS <- correlog_plot(score = "NS_score", increment = 500, resamp= 30)
spatial_cor_NS
ggsave(plot = spatial_cor_NS, filename = here::here("outputs", "figures",
          "Spatial correlation per distance NS, increment 500.png"))

spatial_cor_NS <- correlog_plot(score = "NS_score", increment = 100, resamp= 30)
spatial_cor_NS
ggsave(plot = spatial_cor_NS, filename = here::here("outputs", "figures",
            "Spatial correlation per distance NS, increment 100.png"))

spatial_cor_NS <- correlog_plot(score = "NS_score", increment = 200, resamp= 99)
spatial_cor_NS
ggsave(plot = spatial_cor_NS, width=6, height =6,
       filename = here::here("outputs", "figures","Spatial correlation per distance NS, increment 200.png"))

#All contributions individually
parallel::mclapply(colnames(all_data)[c(8:37)], mc.cores = 15, FUN=function(i){
  spatial_cor_NN <- correlog_plot(score = i, increment = 200, resamp= 99)
  ggsave(plot = spatial_cor_NN, width=6, height =6,
         filename = here::here("outputs", "figures", "Spatial correlation contributions",
                               paste0("Spatial correlation in ", i, ".png")))
})



##-------------Variogram-------------
dist_matrix <- geodist::geodist(x= coord[,c("longitude", "latitude")],
                              measure = "geodesic") #distance in meters
range(dist_matrix) #max ~20 000 km => OK

site_dist <- dist_matrix[lower.tri(dist_matrix, diag = FALSE)]
vario <- nlme::Variogram(unique_points$NN_score, distance = site_dist)
# plot(vario, smooth = TRUE)
variogram <- ggplot(data= vario) +
  geom_point(aes(y = variog, x = dist),
             col = "grey", alpha = 0.2, size = 1) +
  geom_smooth(aes(y = variog, x = dist))+
  theme_bw()
ggsave(plot = variogram, filename = here::here("outputs", "figures",
           "Spatial correlation variogram.png"))


# ##-------------check  very distant points >18 000 km -------------
# long_dist <- which(dist_matrix > 19000000)
# names <- c()
# for( i in long_dist){
#   r <- i %% nrow(dist_matrix)
#   c <- i %/% nrow(dist_matrix) + 1
#   names <- c(names, paste0(coord$SiteCode[max(r,c)], "/", coord$SiteCode[min(r,c)]))
# }
# names_u <- unique(names)
# names_u
