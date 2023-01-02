#' Aggregating productivity at transect level, and logging it 
#' 
#' @param data_with_prod Output from the calc_prod function
#' @param transect_info transect area
#' 
#' @return dataframe with productivity calculated at site level or community level
#' 

calc_prod_transect <- function(data_with_prod,transect_info){
  
  #transect_info <- read.table("data/RLS_transect_info.txt")
  
  transect_info = transect_info %>% 
    dplyr::select(SurveyID, SiteLatitude, SiteLongitude, SiteCode, Depth = SurveyDepth, Country = SiteCountry)
  
  transect_site = transect_info %>% dplyr::select(SurveyID, SiteCode) %>%
    dplyr::filter(SurveyID %in% data_with_prod$SurveyID) 
  
  data_with_prod = data_with_prod %>% left_join(transect_site, by = "SurveyID")
  
  
  # At the scale of the community (transect)
  data_prod_brut = data_with_prod %>%
    #Sum for each transect
    group_by(SurveyID) %>%
    mutate(Biom = sum(Biom)/500,
           Prod = sum(Prod)/500,
           Productivity = (Prod/Biom)*100) %>%
    ungroup() %>%
    # #Mean for each site
    # group_by(SiteCode) %>%
    # mutate(Biom = mean(Biom),
    #        Prod = mean(Prod),
    #        Productivity = mean(Productivity)) %>%
    # ungroup() %>%
    #joinin with transect data
    dplyr::select(SurveyID, Biom, Prod, Productivity) %>%
    distinct(SurveyID, .keep_all = T) %>%
    left_join(transect_info, by ="SurveyID") %>%
    #Transforming data
    mutate(log10ProdB = Productivity,
           log10Biom = log10(Biom+1),
           log10Prod = log10(Prod+1),
           SiteLatitude = as.numeric(as.character(SiteLatitude)),
           SiteLongitude = as.numeric(as.character(SiteLongitude))) %>%
    dplyr::rename(site_code = SiteCode) %>%
    dplyr::select(SurveyID, site_code, Biom, Prod, Productivity, SiteLatitude, SiteLongitude, Depth, Country, log10ProdB, log10Biom, log10Prod)

  
  # data_prod_brut <-  data_with_prod[,c("SurveyID","Biom","Prod")]
  # data_prod_brut_sites <-  data_with_prod[,c("SurveyID","SiteCode")]
  # data_prod_brut <- aggregate(. ~ SurveyID, data = data_prod_brut, sum, na.rm=T)
  # data_prod_brut$Productivity  <- data_prod_brut$Prod/data_prod_brut$Biom
  # data_prod_brut <- left_join(data_prod_brut,data_prod_brut_sites,by="SurveyID")
  # data_prod_brut <- aggregate(. ~ SiteCode, data = data_prod_brut[-1], mean, na.rm=T)
  # 
  # # # # #At the scale of the community (transect)
  # # data_prod_brut <-  data_with_prod[,c("SiteCode","Biom","Prod")]
  # # data_prod_brut_sites <-  data_with_prod[,c("SurveyID","SiteCode")]
  # # data_prod_brut <- aggregate(. ~ SiteCode, data = data_prod_brut, sum, na.rm=T)
  # # data_prod_brut$Productivity  <- data_prod_brut$Prod/data_prod_brut$Biom
  # # # data_prod_brut <- left_join(data_prod_brut,data_prod_brut_sites,by="SiteCode")
  # # # data_prod_brut <- aggregate(. ~ SiteCode, data = data_prod_brut[-5], mean, na.rm=T)
  # # 
  # data_prod_brut = left_join(data_prod_brut, transect_info, by ="SiteCode")
  # 
  # # data_prod_brut$log10ProdB <- log10(data_prod_brut$Productivity+1)
  # data_prod_brut$log10ProdB <-data_prod_brut$Productivity
  # data_prod_brut$log10Biom <-log10(data_prod_brut$Biom+1)
  # data_prod_brut$log10Prod<-log10(data_prod_brut$Prod+1)
  # data_prod_brut$SiteLatitude <- as.numeric(as.character(data_prod_brut$SiteLatitude))
  # data_prod_brut$SiteLongitude <- as.numeric(as.character(data_prod_brut$SiteLongitude)) 
  # 
  # data_prod_brut = data_prod_brut %>% dplyr::rename(site_code = "SiteCode")

  return(data_prod_brut)

}



calc_prod_site <- function(data_with_prod,transect_info){
  
  #transect_info <- read.table("data/RLS_transect_info.txt")
  
  transect_info = transect_info %>% 
    dplyr::select(SurveyID, SiteLatitude, SiteLongitude, SiteCode, Depth = SurveyDepth, Country = SiteCountry)
  
  transect_site = transect_info %>% dplyr::select(SurveyID, SiteCode) %>%
    dplyr::filter(SurveyID %in% data_with_prod$SurveyID) 
  
  data_with_prod = data_with_prod %>% left_join(transect_site, by = "SurveyID")
  
  
  # At the reef scale (site)
  data_prod_site_brut = data_with_prod %>%
    #Sum for each transect
    group_by(SurveyID) %>%
    mutate(Biom = sum(Biom)/500,
           Prod = sum(Prod)/500,
           Productivity = (Prod/Biom)*100) %>%
    ungroup() %>%
    #Mean for each site
    group_by(SiteCode) %>%
    mutate(Biom = mean(Biom),
           Prod = mean(Prod),
           Productivity = mean(Productivity)) %>%
    ungroup() %>%
    #joinin with transect data
    dplyr::select(SiteCode, Biom, Prod, Productivity) %>%
    distinct(SiteCode, .keep_all = T) %>%
    left_join(transect_info, by ="SiteCode") %>%
    #Mean depth
    group_by(SiteCode) %>%
    mutate(Depth = mean(Depth)) %>%
    ungroup() %>%
    #Transforming data
    mutate(log10ProdB = Productivity,
           log10Biom = log10(Biom+1),
           log10Prod = log10(Prod+1),
           SiteLatitude = as.numeric(as.character(SiteLatitude)),
           SiteLongitude = as.numeric(as.character(SiteLongitude))) %>%
    dplyr::rename(site_code = SiteCode) %>%
    dplyr::select(site_code, Biom, Prod, Productivity, SurveyID, SiteLatitude, SiteLongitude, Depth, Country, log10ProdB, log10Biom, log10Prod) %>%
    distinct(site_code, .keep_all = T)

  
  return(data_prod_site_brut)
  }