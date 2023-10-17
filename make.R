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


#-----------------Assess all reef contributions---------------------

## Assess the recycling of N and P from tropical fishes in each communities
source(here::here("recycling", "analyses", "recycling_analysis.R"))

## Assess the productivity in each community
source(here::here("productivity", "analyses", "productivity_analyses.R"))

## Assess the nutrient concentration and the fishery biomass in each community
source(here::here("nutrients", "analyses", "nutrients_analyses.R"))

## Assess the  in each community
source(here::here("", "analyses", "productivity_analyses.R"))

## Assess the  in each community
source(here::here("", "analyses", "productivity_analyses.R"))

## Assess the  in each community
source(here::here("", "analyses", "productivity_analyses.R"))

