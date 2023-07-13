################################################################################
#' 
#' Extract the RLS data of the period 2006-2019 gathered for the REEF 
#'  FUTURE project (this needs permission access from S. Villéger)
#' Functions used in the scripts of the script `analysis/preping_data_rls.R`
#' 
#'
#' @author Sebastien Villéger
#'
#' 
################################################################################
## cleaning memory
rm(list=ls())

## preparing environment ######
# token to authorize googledrive to reach online files
drive_auth(email="sebastien.villeger@univ-montp2.fr")


# directories where to store data (=inputs) and results (=outputs) of R scripts
dir_data_source<-"data_raw/source"

## downloading data from googledrive repository on 2021-02-02 ######

# downloading Fish data
link_fish<-"https://drive.google.com/file/d/1WAeDFL22p1H-oOUntDhmemn2Z5yyJqZ2/"
drive_download( file=link_fish, path = file.path(dir_raw_data, "RLS_fish.rds"),
                type = NULL, overwrite = TRUE, verbose = TRUE)

# downloading species data
link_species<-"https://drive.google.com/file/d/1YpTGA4tqq4YUjqpmhYbMPOtr_3niYNSY/"
drive_download( file=link_species, path = file.path(dir_raw_data, "Species_List_RLS_Jan2021.csv"),
                type = NULL, overwrite = TRUE, verbose = TRUE)


# downloading sites data
link_sites<-"https://drive.google.com/file/d/1u7LJaZ_CC3iC_Sce2QWDUz3fNMp39t55/"
drive_download( file=link_sites, path = file.path(dir_raw_data, "RLS_sites.rds"),
                type = NULL, overwrite = TRUE, verbose = TRUE)

# downloading environmental data
link_env<-"https://drive.google.com/file/d/1gT_hye6jQZPe0sT3Lzrb_bx35sTx4CDH/"
drive_download( file=link_env, path = file.path(dir_raw_data, "RLS_env_spatio_temporal.rds"),
                type = NULL, overwrite = TRUE, verbose = TRUE)


# file with trait data 
# downloaded from Nicolas L github on 2021-03-18
# "traits.Rdata"

# files with data to run fishflux 
# shared by Nina on 2021-09-01
# "species_parameters.csv" including diet from Parravicini et al 2020 paper
# "metpar_fam_smr.csv" , "kmax_combined.csv", "ae_dietcat.csv"

## END of source data importation #####







