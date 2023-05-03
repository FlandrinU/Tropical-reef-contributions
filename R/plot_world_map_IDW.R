################################################################################
##
## Plot Scores on world map with buffer
##
## plot_world_map_IDW.R
##
## 29/03/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "sf", "raster", "geosphere")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data------------
load(here::here("outputs","NN_NS_score_wheighted_mean.Rdata"))
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))
##-------------set parameters------------
data <- NN_NS_scores
score <- "NN_score"
res <- c(1000, 1000) #resolution of map x, y
my_buffer <- 200*1000 # The size of the buffers around the points in metres
p <- 0.02 #power for IDW
#1.6 for ella Clausius
#0.5 was the default in the old code 
#lower values mean its more influenced by distances further away
# higher values = less smoothing 


##-------------clean data------------
sites <- data |>
  dplyr::select(SiteCode, SiteLongitude, SiteLatitude, all_of(score)) |>
  na.omit()

# sites <- sites |> dplyr::filter(SiteLatitude >-39 & SiteLatitude < 0 &
#                                   SiteLongitude > 100 & SiteLongitude <180) ###########################"
sites_xy <- cbind(sites$SiteLongitude, sites$SiteLatitude)

##-------------Use raster to make buffers arounf sites------------
#create raster with resolution defined above 
r <- raster::raster(nrow = res[1], ncol = res[2]) 

#pulls raster cells that match coordinates from sites_xy
icell <- raster::cellFromXY(r, sites_xy) 

#assigns grid cells present in icell a value of 1 
r[icell] <- 1 

#computes the distance for all cells that are NA to the nearest cell that is not NA (i.e., to the sites in sites_xy)
# rdist <- raster::distance(r) ##this may take a while to compute if working with a large number of sites 
# save(rdist, file = here::here("outputs", "raster_with_distance_from_RLS_points.Rdata"))
load(file = here::here("outputs", "raster_with_distance_from_RLS_points.Rdata"))


#list of grid cells in rdist that fall within the defined buffer 
ifilter <- which(rdist[] < my_buffer)  

#generate a list of coordinates for each of the grid cells that fall within the buffer of each site
xyinterp <- raster::xyFromCell(r, ifilter) 

#calculate the distance of each grid cell inside the buffer from each site
xydist <- geosphere::distm(sites_xy, xyinterp)/1000 

#pull values of interest: "score"
values <- sites |> 
  dplyr::pull(score)

#weights each grid cell value by it's distance from the site according to the power you assign. 
w <- 1/(xydist^p) 

indic <- matrix(values, nrow = 1)
isreal <- which(!is.na(indic))
nvals <- length(isreal)

val <- (indic[isreal] %*% w[isreal,])/colSums(w[isreal,]) # matrix multiplication(score of sites x weight_distance) / sum of weights for each 
val <- as.numeric(val)

rval <- r
rval[ifilter] <- as.numeric(val)
value_rast <- rval 

#create dataframe of interpreted points
dat_global <- data.frame(xyinterp, val)
names(dat_global) <- c("SiteLongitude", "SiteLatitude", score)

##-------------plot map------------
#Raster
plot(value_rast)

#Ggplot2
library(ggplot2)
map_NN_or_NS <- function(coord_NN_NS = NN_NS_with_product,
                         NCP = NN_NS_with_product$NN_score ,
                         sites= NN_NS_scores,
                         col_NCP= "forestgreen",
                         ylim = c(-36, 31),
                         xlim= c(-180,180), title="",
                         title_legend = ""){
  ggplot(coord_NN_NS) +
    geom_tile(aes(x = SiteLongitude, y = SiteLatitude,
                   fill = NCP)) +
    # scale_fill_gradientn(colours = colorRampPalette(rev(c( col_NCP ,"white", "grey30")))(1000), midpoint=0) +
    scale_fill_gradient2(low="grey30", mid="white", high=col_NCP, midpoint = 0)+
    geom_sf(data = coast, color = NA, fill = "lightgrey") +
    coord_sf(xlim=xlim, ylim = ylim, expand = FALSE) +
    
    geom_point( aes(x = SiteLongitude, y = SiteLatitude),
                shape = 1, size = 2, stroke = 0.5,
                color= "black",
                data = head(sites[order(NN_NS_scores$NS_score, decreasing = T),], 10)) +##########################change NN_score
    # geom_point( aes(x = SiteLongitude, y = SiteLatitude),
    #             shape = 4, size = 2, stroke = 0.5,
    #             color= "black",
    #             data =sites) +
    
    theme_minimal()+
    labs(title = title, fill = title_legend,
         x="Longitude", y= "Latitude") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme(
      # legend.position = "none",
      legend.title = element_text(size = 10),
      plot.title = element_text(size=10, face="bold"))
}

## NN
map_NN_or_NS(coord_NN_NS = dat_global,
             NCP = dat_global$NN_score,
             sites=NN_NS_scores,
             col_NCP= "forestgreen",
             ylim = c(-36, 31),
             xlim= c(-180,180) ,
             title="",
             title_legend = "Nature to \n Nature")

map_NN_or_NS(coord_NN_NS = dat_global,
             NCP = dat_global$NN_score,
             sites= NN_NS_scores,
             col_NCP= "forestgreen",
             ylim = c(-36, 0),
             xlim= c(130,180) ,
             title="",
             title_legend = "Nature to \n Nature")

# ggsave( here::here("outputs", "figures", "world map with NN score.png"), plot = NN_worldwide, width=10, height = 6 )
# 
# ## NS
# map_NN_or_NS(coord_NN_NS = dat_global,
#              NCP = dat_global$NS_score,
#              sites=NN_NS_scores,
#              col_NCP= "dodgerblue3",
#              ylim = c(-36, 31),
#              xlim= c(-180,180) ,
#              title="",
#              title_legend = "Nature to \n People")
# 
# map_NN_or_NS(coord_NN_NS = dat_global,
#              NCP = dat_global$NS_score,
#              col_NCP= "dodgerblue3",
#              ylim = c(-36, 0),
#              xlim= c(130,180) ,
#              title="",
#              title_legend = "Nature to \n People")
# 
# 
# map_NN_or_NS(coord_NN_NS = dat_global,
#              NCP = dat_global$NS_score,
#              col_NCP= "dodgerblue3",
#              ylim = c(-6, 20),
#              xlim= c(-90,-60) ,
#              title="",
#              title_legend = "Nature to \n People")
# 
# 
