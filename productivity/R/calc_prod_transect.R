#' Aggregating productivity at transect level, and logging it 
#' 
#' @param data_with_prod Output from the calc_prod function
#' @param transect_info transect area
#' 
#' @return dataframe with productivity calculated at survey level 
#' 


calc_prod_transect <- function(data_with_prod,transect_info){
  
  #transect_info <- read.table("data/RLS_transect_info.txt")
  
  transect_info = transect_info |> 
    dplyr::select(SurveyID, SiteLatitude, SiteLongitude, SiteCode, Depth = SurveyDepth, Country = SiteCountry)
  
  transect_site = transect_info |> dplyr::select(SurveyID, SiteCode) |>
    dplyr::filter(SurveyID %in% data_with_prod$SurveyID) 
  
  data_with_prod = data_with_prod |> dplyr::left_join(transect_site, by = "SurveyID")
  
  
  # At the scale of the community (transect)
  data_prod_brut = data_with_prod |>
    #Sum for each transect
    dplyr::group_by(SurveyID) |>
    dplyr::mutate(Biom = sum(Biom)/500,
                  Prod = sum(Prod)/500,
                  Productivity = (Prod/Biom)*100) |>
    dplyr::ungroup() |>
    #joinin with transect data
    dplyr::select(SurveyID, Biom, Prod, Productivity) |>
    dplyr::distinct(SurveyID, .keep_all = T) |>
    dplyr::left_join(transect_info, by ="SurveyID") |>
    #Transforming data
    dplyr::mutate(log10ProdB = Productivity,
                  log10Biom = log10(Biom+1),
                  log10Prod = log10(Prod+1),
                  SiteLatitude = as.numeric(as.character(SiteLatitude)),
                  SiteLongitude = as.numeric(as.character(SiteLongitude))) |>
    dplyr::rename(site_code = SiteCode) |>
    dplyr::select(SurveyID, site_code, Biom, Prod, Productivity, SiteLatitude, SiteLongitude, Depth, Country, log10ProdB, log10Biom, log10Prod)
  
  
  return(data_prod_brut)
  
}




#' Aggregating productivity at site level, and logging it 
#' 
#' @param data_with_prod Output from the calc_prod function
#' @param transect_info transect area
#' 
#' @return dataframe with productivity calculated at community (site) level
#' 



calc_prod_site <- function(data_with_prod,transect_info){
  
  #transect_info <- read.table("data/RLS_transect_info.txt")
  
  transect_info = transect_info |> 
    dplyr::select(SurveyID, SiteLatitude, SiteLongitude, SiteCode, Depth = SurveyDepth, Country = SiteCountry)
  
  transect_site = transect_info |> dplyr::select(SurveyID, SiteCode) |>
    dplyr::filter(SurveyID %in% data_with_prod$SurveyID) 
  
  data_with_prod = data_with_prod |> dplyr::left_join(transect_site, by = "SurveyID")
  
  
  # At the reef scale (site)
  data_prod_site_brut = data_with_prod |>
    #Sum for each transect
    dplyr::group_by(SurveyID) |>
    dplyr::mutate(Biom = sum(Biom)/500,
                  Prod = sum(Prod)/500,
                  Productivity = (Prod/Biom)*100) |>
    dplyr::ungroup() |>
    #Mean for each site
    dplyr::group_by(SiteCode) |>
    dplyr::mutate(Biom = mean(Biom),
                  Prod = mean(Prod),
                  Productivity = mean(Productivity)) |>
    dplyr::ungroup() |>
    #joinin with transect data
    dplyr::select(SiteCode, Biom, Prod, Productivity) |>
    dplyr::distinct(SiteCode, .keep_all = T) |>
    dplyr::left_join(transect_info, by ="SiteCode") |>
    #Mean depth
    dplyr::group_by(SiteCode) |>
    dplyr::mutate(Depth = mean(Depth)) |>
    dplyr::ungroup() |>
    #Transforming data
    dplyr::mutate(log10ProdB = Productivity,
                  log10Biom = log10(Biom+1),
                  log10Prod = log10(Prod+1),
                  SiteLatitude = as.numeric(as.character(SiteLatitude)),
                  SiteLongitude = as.numeric(as.character(SiteLongitude))) |>
    dplyr::rename(site_code = SiteCode) |>
    dplyr::select(site_code, Biom, Prod, Productivity, SurveyID, SiteLatitude, SiteLongitude, Depth, Country, log10ProdB, log10Biom, log10Prod) |>
    dplyr::distinct(site_code, .keep_all = T)
  
  
  return(data_prod_site_brut)
}