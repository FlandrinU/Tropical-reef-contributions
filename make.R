#' Run the Entire Project
#'
#' This script runs the entire project and produces all figures present in the
#' Flandrin U _et al._ 2023.
#'
#' @author Ulysse Flandrin, \email{ulysse.flandrin@@gmail.com}
#'
#' @date 2023/06/13
#' 
#' @note
#' **WARNING:** Running this script calls all the scripts and takes time so if 
#' you want to work on one or a few scripts, you should run lines 17-42 of this 
#' script and then go to the other one.


## Clean environment ----

rm(list = ls())


#-----------------Install required dependencies---------------------

if (!("remotes" %in% utils::installed.packages())) 
  install.packages("remotes")
renv::install("r-lib/devtools") #to install packages from github
renv::install() #install all packages noted in file DESCRIPTION
renv::snapshot()


## Load packages & functions + Setup project ----

devtools::load_all(here::here())



#-----------------Preping data---------------------

## save a world coastline shapefile
behrmann <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs" # Behrmann projection - Pacific-centered ----
world <- robinmap::robinmap(center = 150, crs = behrmann)
sf::st_write(obj=world, here::here("data", "ShapeFiles coast", "shapefile_coast_pacific_centered.shp"))

## Extract, filter and clean  data to run analysis
source(here::here("analyses", "preping_data_rls.R"))


#-----------------Assess all survey contributions---------------------

## Assess the recycling of N and P from tropical fishes in each communities
source(here::here("recycling", "analyses", "recycling_analysis.R"))

## Assess the biodiversity metrics in each community
source(here::here("biodiversity", "analyses", "productivity_analyses.R"))

## Assess the productivity in each community
source(here::here("productivity", "analyses", "productivity_analyses.R"))

## Assess the nutrient concentration and the fishery biomass in each community
source(here::here("nutrients", "analyses", "nutrients_analyses.R"))

## Assess the cultural contributions of each community
source(here::here("cultural_contributions", "analyses", "cultural_analyses.R"))

## Assess the carbonate excretion in each community
source(here::here("carbonates", "analyses", "carbonates_analyses.R"))

## Assess the trophic structure in each community
source(here::here("trophic_web", "trophic_web_analyses.R"))



#-----------------Study contributions at the reef (=site) level---------------------

## Merge assessed contributions at the reef level
source(here::here('R', 'merge_NCP_surveys.R'))

## Explore NCP
source(here::here('R', 'Plot_NCP.R'))

## Study the dimensionnality of contributions
source(here::here('R', 'PCA_analyses_on_NCP.R'))



#----------------- 'Nature-for-Nature' and 'Nature-for-People' framework ---------------------

## Sum up contributions into 2 dimensions: NN and NP
source(here::here('R', 'weighted_mean_NP_NN_score.R'))

## Test the spatial correlation
source(here::here('R', 'spatial_autocorrelation.R'))



#----------------- Consttruct the figures of the paper ---------------------

source(here::here('R', 'make_fig_1.R'))
source(here::here('R', 'make_fig_2.R'))
