### filtering to keep only RLS data for detectable fishes ######

## preparing environment ######

# cleaning memory
rm(list=ls())

## loading data from RLS and from fishflux ####
load(here::here("data_raw", "rls_trop_fish.Rdata"))
head(rls_trop_fish)


## checking maximum fish size ####        

# looking for large fish to check for errors ----
# as they can have large effect on biomass
fish_mlsupp80cm<-rls_trop_fish |> 
  dplyr::filter(MaxLength>80) |>
  dplyr::distinct(taxa, MaxLength) |>
  dplyr::arrange(taxa)
nrow(fish_mlsupp80cm)
# 146 cases
# have a look to fish_mlsupp80cm table

# groupers species for which max size is too low ----
# Epinephelus lanceolatus is not 85.10 but 270
# Epinephelus spp is not 142 but 270
rls_trop_fish[which(rls_trop_fish$taxa=="Epinephelus_lanceolatus"),"MaxLength"]<-270
rls_trop_fish[which(rls_trop_fish$taxa=="Epinephelus_spp."),"MaxLength"]<-270


# filling missing values of Max length ----

# number of cases
sp_MaxLengthNA<-rls_trop_fish |>
  dplyr::filter (is.na(MaxLength) ) |>
  dplyr::distinct(taxa, genus, family, MaxLength)
# open sp_MaxLengthNA 
nrow(sp_MaxLengthNA) # => 110 taxa, many Genus spp.


# filling missing values of Maximum length
# based on maxlength of species from same genus or family
for (k in 1:nrow(sp_MaxLengthNA) ) { 
  
  # median of species from same genus
  samegenus_k<-rls_trop_fish |>
    dplyr::filter( genus==sp_MaxLengthNA[k,"genus"] &
              (!is.na(MaxLength))  ) |>
    dplyr::distinct(taxa, MaxLength)
  
  if (length(samegenus_k$MaxLength)>=1) {
    sp_MaxLengthNA[k,"MaxLength"]<-median(samegenus_k$MaxLength)
  }else {
    # median of species from same family
    samefamily_k<-rls_trop_fish |>
      dplyr::filter( family==sp_MaxLengthNA[k,"family"] &
                (!is.na(MaxLength))  ) |>
      dplyr::distinct(taxa, MaxLength)
    
    if (length(samefamily_k$MaxLength)>=1) {
      sp_MaxLengthNA[k,"MaxLength"]<-median(samefamily_k$MaxLength)
    } else {
      # NA
      sp_MaxLengthNA[k,"MaxLength"]<-NA}
  }
  
}# end of k    

# checking remaining NA  
sp_MaxLengthNA |> 
  dplyr::filter(is.na(MaxLength)) # 3 taxa

# Latropiscis_purpurissatus filling based on data from Fishbase  
sp_MaxLengthNA[which(sp_MaxLengthNA$taxa=="Latropiscis_purpurissatus"),"MaxLength"]<-60

# for Aseraggodes_spp. based on observed size *2
rls_trop_fish |>
  dplyr::filter( family=="Soleidae")
sp_MaxLengthNA[which(sp_MaxLengthNA$taxa=="Aseraggodes_spp."),"MaxLength"]<-20

# for Pseudomugil_spp. based on observed size *2
rls_trop_fish |>
  dplyr::filter( family=="Pseudomugilidae")
sp_MaxLengthNA[which(sp_MaxLengthNA$taxa=="Pseudomugil_spp."),"MaxLength"]<-5


# filling missing values
for (k in 1:nrow(sp_MaxLengthNA) ) { 
  rls_trop_fish[which(rls_trop_fish$taxa==sp_MaxLengthNA[k,"taxa"]), 
                "MaxLength"]<-sp_MaxLengthNA[k,"MaxLength"]
}# end of k


# checking
summary(rls_trop_fish$MaxLength) # OK no NA remaining



## removing small fish ####

# species with size == 2.5cm ----
taxa_small<-rls_trop_fish |>
  dplyr::filter(rls_trop_fish$size_class==2.5)
dplyr::n_distinct(taxa_small$taxa) # 642 taxa
dplyr::n_distinct(taxa_small$family) # 39 families
# => small cryptobenthic or juveniles of larger fish

# removing fish<2.5cm, and juveniles (5cm) of fish larger than 25cm ----
# keeping missing size (except if Maxlength<5cm)

# keeping fish larger than 5cm and fish of 5cm if Maxsize <25cm
rls_trop_fish_nosmall<-rls_trop_fish |>
  dplyr::filter ( (is.na(size_class) & MaxLength>=5) | 
             (size_class>5 & MaxLength >=5 ) | 
             (size_class==5 & MaxLength<25) )

# diversity remaining
dplyr::n_distinct(rls_trop_fish_nosmall$family) # 84 families
dplyr::n_distinct(rls_trop_fish_nosmall$taxa) # 1 697 taxa
summary(rls_trop_fish_nosmall$size_class) # from 5 to 300cm with 146 NA


## surveys with missing size for not small fishes ####   
sizeNA<-rls_trop_fish_nosmall |>
  dplyr::filter ( is.na(size_class) )
summary(sizeNA$number) # from 1 to 905 individuals, median =3
summary(sizeNA$MaxLength) # Q1=12cm ; Q3>30cm

surveys_sizeNA <- sizeNA |>
  dplyr::distinct(SurveyID)
nrow(surveys_sizeNA) # 43 surveys

# removing transects with missing size data ----
rls_trop_fish_nosmall <- rls_trop_fish_nosmall |>
  dplyr::filter ( !(SurveyID %in% surveys_sizeNA$SurveyID ) )

# diversity remaining
dplyr::n_distinct(rls_trop_fish_nosmall$SurveyID) # 3 751 surveys
dplyr::n_distinct(rls_trop_fish_nosmall$family) # 83 families
dplyr::n_distinct(rls_trop_fish_nosmall$taxa) # 1 694 taxa



## looking at fish size too big to not be error ####

# names of all taxa
taxa_trop<- unique(rls_trop_fish_nosmall$taxa)
length(taxa_trop) # 1 694 taxa

#  outliers for size within each species ----

# summary of maxLength
summary(rls_trop_fish_nosmall$MaxLength) 
# from 2.5 (while size filtered to be >5 ) to 333

rls_trop_fish_nosmall |> dplyr::filter(MaxLength<5) |> dplyr::distinct(family, taxa)
# small fish from only 8 taxa

# threshold = 2*(maxLength) or (maxlength+50)
list_outlier_size<-list()

# loop on taxa
for (k in taxa_trop) {
  # observations
  ind_k<- rls_trop_fish_nosmall |>
    dplyr::filter(taxa==k)
  
  ind_k<- ind_k |>
    dplyr::select(taxa, family, MaxLength,  SurveyID, size_class, number)
  
  list_outlier_size[[k]]<-dplyr::filter(ind_k, 
                                 (size_class > (MaxLength)*2 )  | 
                                 (size_class > (MaxLength+50) )
  )
  
} # end of k

# list to dataframe and exporting as csv ----
outlier_size<-do.call(rbind.data.frame, list_outlier_size)
row.names(outlier_size)<-NULL
write.csv( outlier_size, here::here("data_raw", "outlier_size.csv") )
# open outlier_size.csv and reading all rows
# => many small fishes (likely SL in mm not well converted to cm)


# summary of potential errors
nrow(outlier_size) # 132 cases
dplyr::n_distinct(outlier_size$SurveyID) # 112 transects
unique(outlier_size$family) # 22 families
outlier_size |> dplyr::filter(MaxLength>30) # 14 cases of error on large fish
summary(outlier_size$number) # >50% cases with 1 indiv
outlier_size |> dplyr::filter(number>10) # 9 cases with >10 indiv for small species


# removing observations with erroneous size data ----
rls_trop_fish_sizeok <- rls_trop_fish_nosmall |>
  dplyr::filter( ! (SurveyID %in% unique(outlier_size$SurveyID) ) )

# diversity remaining
dplyr::n_distinct(rls_trop_fish_sizeok$SurveyID) # 3 639 surveys
dplyr::n_distinct(rls_trop_fish_sizeok$family) # 81 families
dplyr::n_distinct(rls_trop_fish_sizeok$taxa) # 1 684 taxa


## saving as Rdata ####
save( rls_trop_fish_nosmall, file=here::here("data_raw", "rls_trop_fish_nosmall.Rdata")  )
save( rls_trop_fish_sizeok, file=here::here("data_raw", "rls_trop_fish_sizeok.Rdata")  )

