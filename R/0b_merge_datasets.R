#### exploring RLS data, recoding some variables, merging datasets ######

## preparing environment ######

# cleaning memory
rm(list=ls())

## metadata of surveys ####

  # loading raw dataset ----
  rls_surveys<-readRDS(here::here("data_raw", "source", "RLS_sites.rds"))
  head(rls_surveys)
    
  # recoding ID as character string ----
  class(rls_surveys$SurveyID) # integer but better as a character
  rls_surveys$SurveyID<-as.character(rls_surveys$SurveyID)
  n_distinct(rls_surveys$SurveyID) # 7100 surveys
  
  # site code ----
  class(rls_surveys$SiteCode) # factor
  n_distinct(rls_surveys$SiteCode) # 3544 sites


## environment #####
  
  # loading raw dataset ----
  rls_env<-readRDS(here::here("data_raw", "source", "RLS_env_spatio_temporal.rds"))
  head(rls_env)

 #plot(rls_env$min_sst_1year, rls_env$min_sst_5year)
 rls_env %>% filter(min_sst_1year>=17) %>% nrow()
 rls_env %>% filter(min_sst_5year>=17) %>% nrow()
  # keeping 5 years average, to be more representative
    
  # keeping only SST ----
  rls_env_sst<- rls_env %>%
    select(SurveyID, 
           SiteMinSST=min_sst_5year, 
           SiteMeanSST=mean_sst_5year, 
           SiteMaxSST=max_sst_5year)
  rls_env_sst$SurveyID<-as.character(rls_env_sst$SurveyID) # recoding variable
  
  n_distinct(rls_env_sst$SurveyID) # 7100 surveys
  identical( sort( unique(rls_env_sst$SurveyID)), 
             sort( unique(rls_surveys$SurveyID) ) ) # OK
  
## UVC database #####

  # loading raw dataset ----
  rls_uvc<-readRDS(here::here("data_raw", "source", "RLS_fish.rds"))
  head(rls_uvc)
  
  # ID of surveys recoded and checked ---
  class(rls_uvc$SurveyID) # integer but better as a character
  rls_uvc$SurveyID<-as.character(rls_uvc$SurveyID) # recoding variables
  n_distinct(rls_uvc$SurveyID) # 7100 surveys as in rls_surveys
  
  identical( sort( unique(rls_uvc$SurveyID)), 
             sort( unique(rls_surveys$SurveyID) ) ) # OK
  
  # Check for duplicates ---
  nrow(rls_uvc) # 380132
  n_distinct(rls_uvc) # 380118
  # There are 14 repeated observations
  rls_uvc %>% 
    group_by(across()) %>%
    count() %>% 
    filter(n > 1)
  # Mainly observations from a species of wrasse (Ophthalmolepis lineolata)
  # Remove these
  rls_uvc <- unique(rls_uvc)
  nrow(rls_uvc) # 380 118 observations
  
  # taxonomy ----
  n_distinct(rls_uvc$SPECIES_NAME) # 2832 "species" in original dataset), useless
  class(rls_uvc$TAXONOMIC_NAME) # factor
  n_distinct(rls_uvc$TAXONOMIC_NAME) # 2801 valid taxonomic names
  nrow (filter(rls_uvc, is.null(TAXONOMIC_NAME)) ) # all obs have a valid name
  
  # number of individuals ----
  class(rls_uvc$Num) # number of individuals is a numeric variable (integer)
  summary(rls_uvc$Num) # min=0, no NAs
  rls_uvc %>% filter(Num==0) # only one obs with no individual, Error ?
  
  # size class ----
  class(rls_uvc$Sizeclass) # size classes are categories
  levels(rls_uvc$Sizeclass) # numeric values as character + "NULL" category
  # recoding as numeric (median of size class), NULL replaced by NA
  rls_uvc$Sizeclass <- as.numeric(as.character(rls_uvc$Sizeclass))
  summary(rls_uvc$Sizeclass) # 236 NAs
  
  
  # biomass ----
  class(rls_uvc$Biomass) # factor but should be numeric
  levels(rls_uvc$Biomass) # biomass as character
  # recoding as numeric variable, NULL replaced by NA
  rls_uvc$Biomass <- as.numeric(as.character(rls_uvc$Biomass))
  summary(rls_uvc$Biomass) # 0 and 3308 NAs
  rls_uvc %>% filter(Biomass==0) # only one obs with no individual, Error ?
  
  nrow(rls_uvc) # 380 118 entries
  
## Taxonomy of species #####
  
  # loading raw dataset ----
  rls_species<-read.csv(here::here("data_raw", "source", "Species_List_RLS_Jan2021.csv") )
  head(rls_species)
  
  unique(rls_species$Level) # "species" "genus"   "higher"  "family" 
  
  # unique valid latin name from TAXONOMIC_NAME variable
  nrow(rls_species) # 2847
  n_distinct(rls_species$SPECIES_NAME) # 2 836 species
  n_distinct(rls_species$TAXONOMIC_NAME) # 2 805 species
  rls_species<-distinct(rls_species, TAXONOMIC_NAME,.keep_all = TRUE)
  nrow(rls_species) # 2 805 species
  
  # valid species names recoded ----
  rls_species$TAXONOMIC_NAME<-gsub(" ", replacement ="_", x=rls_species$TAXONOMIC_NAME) # replacing space by "_"
  rls_uvc$TAXONOMIC_NAME<-gsub(" ", replacement ="_", x=rls_uvc$TAXONOMIC_NAME)
  
  # checking match with uvc data ----
  which( (rls_uvc$TAXONOMIC_NAME %in% rls_species$TAXONOMIC_NAME)==FALSE ) 
  # => OK all species in surveys are in species database

## Traits of species #####

  # loading raw dataset ----
  load(here::here("data_raw", "source", "traits.RData"))
  rls_traits<-traits
  head(rls_traits)

  
  # recoding species name and keeping unique valid names ----
  rls_traits$CURRENT_TAXONOMIC_NAME<-gsub(" ", replacement ="_", x=rls_traits$CURRENT_TAXONOMIC_NAME)
  nrow(rls_traits) # 3 189 species
  rls_traits<-distinct(rls_traits, CURRENT_TAXONOMIC_NAME,.keep_all = TRUE)
  nrow(rls_traits) # 3 073 unique valid species
  
  
  # Maxlength ----
  summary(rls_traits$MaxLength) # 268 NA
  
  # recoding variables for 3 traits as factor ----
  # Trophic.group2, Water.column, Diel.Activity
  is.factor(rls_traits$Trophic.group2)
  rls_traits$Trophic.group2<-as.factor(rls_traits$Trophic.group2)
  summary(rls_traits$Trophic.group2)
  
  is.factor(rls_traits$Water.column)
  rls_traits$Water.column<-as.factor(rls_traits$Water.column)
  summary(rls_traits$Water.column)
  
  is.factor(rls_traits$Diel.Activity)
  rls_traits$Diel.Activity<-as.factor(rls_traits$Diel.Activity)
  summary(rls_traits$Diel.Activity)


## Merging data about species #####

  # taxonomy + trait  ----
  rls_species_data<- left_join(
                              select(rls_species, 
                                                  taxa = TAXONOMIC_NAME, 
                                                  taxo_level = Level,
                                                  genus = Genus_FishBase, 
                                                  family = Family_FishBase,
                                                  class = Class_FishBase)
                              ,
                              select(rls_traits,
                                                CURRENT_TAXONOMIC_NAME,
                                                MaxLength,
                                                Diet=Trophic.group2, 
                                                Position=Water.column,
                                                Activity=Diel.Activity
                                                ), 
                                
                                  by=c("taxa"="CURRENT_TAXONOMIC_NAME") 
                                 )

  head(rls_species_data )
  nrow(rls_species_data)  # 2805 taxa
  n_distinct(rls_species$TAXONOMIC_NAME)  # 2805 taxa
  n_distinct(rls_species_data$taxa) # 2805 taxa
  
  
  
  # checking species with missing MaxLength ----
  nrow( rls_species_data %>% 
    filter(is.na(MaxLength)) %>%
    select(taxa, family) )
  # 307 taxa with unknown max length




## Building a single database with all relevant variables from RLS datasets ####

  # keeping only relevant variable from uvc data ----
  rls_all_raw<- rls_uvc %>%
                  select(SurveyID,
                          taxa = TAXONOMIC_NAME, 
                          number = Num,
                          size_class = Sizeclass,
                          biomass= Biomass
                         )
  
  # left merging with species data ----
  rls_all_raw<- rls_all_raw %>% 
                  left_join( rls_species_data,
                             by="taxa")

  # left merging with survey geography data ----
  rls_all_raw <- rls_all_raw %>%
                left_join( select(rls_surveys, 
                                  SurveyID, SurveyDate, SurveyDepth=Depth,
                                  SiteCode, SiteLatitude, SiteLongitude, 
                                  SiteCountry=Country, SiteEcoregion=Ecoregion),
                           by="SurveyID")


  # left merging with SST data ----
  rls_all_raw <- rls_all_raw %>%
                left_join( select(rls_env_sst, 
                                  SurveyID,
                                  SiteMinSST, SiteMeanSST, SiteMaxSST),
                           by="SurveyID")
             

head(rls_all_raw)
summary(rls_all_raw)
nrow(rls_uvc) == nrow(rls_all_raw) # OK
nrow(rls_all_raw) # 380 118 observations
n_distinct(rls_all_raw$SurveyID) # 7100 surveys
n_distinct(rls_all_raw$taxa) # 2 801 taxa
unique(rls_all_raw$taxo_level) #  "genus"   "species" "higher"  "family" 
length(which(rls_all_raw$taxo_level=="species"))/nrow(rls_uvc) # 98.8% of observations at species level

## saving as Rdata ####
save(rls_all_raw, file=here::here("data_raw", "rls_all_raw.Rdata")  )

############## END #############
