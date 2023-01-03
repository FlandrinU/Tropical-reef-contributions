#-----------------Loading packages-------------------
pkgs <- c("here", "dplyr", "stringr", "rfishprod", "stats", "googledrive", "tidyverse",
          "devtools")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
devtools::install_github("renatoamorais/rfishprod")
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

#-----------------renv---------------------
renv::install("r-lib/devtools") #to install packages from github
renv::install() #install all packages noted in file DESCRIPTION
renv::snapshot()
#-----------------Loading all data---------------------

path = (here::here("data"))
setwd(path)
files <- list.files(path, pattern = ".RData|Rdata")
data_list <- lapply(files, load, .GlobalEnv)


#-----------------Loading all functions---------------------
path = (here::here("analyses"))
setwd(path)
sapply(list.files(path), source)

#----------------Run project------------------------
setwd(here::here())


#----------------Loading results------------------------

path = (here::here("outputs"))
setwd(path)
files <- list.files(path, pattern = ".Rdata|RData")
data_list <- lapply(files, load, .GlobalEnv)

#----------------Analyse results------------------------
path = (here::here("analyses"))
setwd(path)
files.source <- list.files(path)
sapply(files.source, source)




