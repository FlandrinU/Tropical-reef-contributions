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
# renv::desactivate()
# renv::activate()

if (!("remotes" %in% utils::installed.packages())) 
  install.packages("remotes")
if (!("renv" %in% utils::installed.packages())) 
  install.packages("renv")
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
source(here::here("biodiversity", "analyses", "biodiversity_analyses.R")) 

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

## Merge all the assessed contributions at the reef level. Filter 
source(here::here('R', '1a_merge_contributions_surveys.R')) 

## Explore the correlations among contributions, and the spatial distribution of
##  studied reef and their contributions
source(here::here('R', '1b_Plot_contributions.R')) 

## Perform a Principal Component Analyses (PCA)on contributions and study the 
##  dimensionnality of contributions
source(here::here('R', '1c_PCA_analyses_on_contributions.R'))



#----------------- 'Nature-for-Nature' and 'Nature-for-People' framework ---------------------

## Synthesize all contributions into 2 dimensions: Nature-for-Nature (NN) and
##  Nature-for-People (NP). Study their distribution and the relationship between them.
source(here::here('R', '1d_weighted_mean_NP_NN_score.R')) 

## Test the spatial correlation of both NN and NP scores
source(here::here('R', '1e_spatial_autocorrelation.R')) 

## Test the the senseibility of composite indicators to the aggregation sensibility
source(here::here('R', '1f_test_composite_scores_NP_NN.R')) 


#----------------- Construct the figures of the paper ---------------------

source(here::here('R', '2a_make_fig_1.R')) 
source(here::here('R', '2b_make_fig_2.R')) 
