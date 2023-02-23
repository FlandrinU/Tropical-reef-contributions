### filtering to keep only RLS data for tropical sites and detectable fishes ######

## preparing environment ######

# cleaning memory
rm(list=ls())

## loading data from RLS and from fishflux ####
load(here::here("data_raw", "rls_all_raw.Rdata"))
head(rls_all_raw)


## filtering to keep only transects from tropical sites #####
# i.e with min(SST)>=17°C

rls_trop<-rls_all_raw |> 
  dplyr::filter( SiteMinSST >= 17)

dplyr::n_distinct(rls_trop$SurveyID) # 3 794 surveys
dplyr::n_distinct(rls_trop$SiteCode) # 2 057 sites
dplyr::n_distinct(rls_trop$taxa) # 2 321 taxa
dplyr::n_distinct(rls_trop$family) # 123 families


##  keeping only actinopterygii ####

unique(rls_trop$class)
# => taxa without class information

# looking at names of NO family and NO class
rls_trop |>
  dplyr::filter(class=="" | family=="") |> 
  dplyr::distinct(taxa, taxo_level, family, class)
# 83 taxa => check, only Actinopterygii except "Dasyatid"  ###### @@@@@ some xxid spp from families not removed in next scripts (e.g. Acanthurid spp, Scarid Scarus spp)

# dplyr::filtering to remove all Elasmobranchii   ###### @@@@@ also removing some Actinopterygii Xxx spp. with missing values for class and family (see above)
rls_trop_fish<- rls_trop |>
  dplyr::filter(class=="Actinopterygii") |> 
  dplyr::filter(taxa!="Dasyatid spp.") |>    
  dplyr::select(-"class")

head(rls_trop_fish)

dplyr::n_distinct(rls_trop_fish$family) # 105 families


## removing species from anguilliforms families ####
# because hard to spot everywhere: 

  # checking Anguilliforms families present
  sort( unique(rls_trop_fish$family) ) # only 3
  fam_anguilliform<-c("Congridae", "Muraenidae", "Ophichthidae")
  
  # filtering
  rls_trop_fish<- rls_trop_fish |>
    dplyr::filter(! family %in% fam_anguilliform)
  
  # diversity remaining
  dplyr::n_distinct(rls_trop_fish$SurveyID) # 3794 surveys
  dplyr::n_distinct(rls_trop_fish$family) # 102 families
  dplyr::n_distinct(rls_trop_fish$taxa) # 2135 taxa
  
  # small cryptobenthic from 17 families
  # 17 families selected by Mattia Ghilardi
  fam_smallcrypto<-c("Aploactinidae", "Apogonidae", "Blenniidae", 
                     "Bythitidae", "Callionymidae", "Chaenopsidae", 
                     "Creediidae", "Dactyloscopidae", "Gobiesocidae", 
                     "Gobiidae", "Grammatidae", "Labrisomidae", 
                     "Opistognathidae", "Plesiopidae", "Pseudochromidae", 
                     "Syngnathidae", "Tripterygiidae")

  # filtering
  rls_trop_fish<- rls_trop_fish |>
    dplyr::filter(! family %in% fam_smallcrypto)
  
  # diversity remaining
  dplyr::n_distinct(rls_trop_fish$SurveyID) # 3794 surveys
  dplyr::n_distinct(rls_trop_fish$family) # 87 families
  dplyr::n_distinct(rls_trop_fish$taxa) # 1713 taxa
  
  
## saving as Rdata ####
save( rls_trop_fish, file=here::here("data_raw", "rls_trop_fish.Rdata")  )
  
  