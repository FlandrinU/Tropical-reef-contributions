#' Estimating biomass production of fishes according to their size, the SST, and
#'  the Kmax of their species
#'
#' @param data_final Dataframe of RLS observation. Each line is an observation 
#'   of fish with its sizeclass and max length, species, the SST,...
#'                    
#'
#' @return return the equivalent dataframe with additional column for biomass 
#'   production and productivity
#' 
#' 


calc_prod_rfishprod <- function(data_final){
  
  # Formula from Morais and Bellwood (2018) #
  fmod <- stats::formula(~ sstmean + MaxSizeTL )
  
  data_final_prod = data_final |> 
    #Selecting columns of interest
    dplyr::select(SurveyID,Num,Family, Species, Sizeclass,Temperature,lwa,lwb,MaxLength,Area) |>
    #Renaming for rfishprod compatibliy
    dplyr::rename(MaxSizeTL = "MaxLength",
                  sstmean = "Temperature",
                  Size = "Sizeclass",
                  a = "lwa",
                  b = "lwb") |>
    dplyr::mutate(Size = ifelse(Size >= MaxSizeTL,MaxSizeTL,Size))
  
  
  # Check dataset repdata #
  # (repdata <- rfishprod:::repdata)
  # 
  # # Getting levels ready #
  # repdata <- rfishprod::tidytrait (repdata, db)
  # data_final_prod = repdata
  
  datagr <- rfishprod::predKmax(data_final_prod,
                                # dataset = db, 
                                fmod = fmod,
                                niter = 100,     #number of xgboost models run for the bootstrap procedure -> 1000 iterations recommended
                                return = 'pred')
  
  datagr <- datagr$pred
  
  # Predicting M/Z: the instantaneous mortality rate (Recommendation: see help file for) #
  datagr$Md <- with (datagr, rfishprod::predM(Lmeas = Size,
                                              Lmax = MaxSizeTL,
                                              Kmax = Kmax,
                                              t = 1,
                                              method = 'Gislason'))
  
  # Positioning your fish in their growth trajectory #
  # aka. what's the size they're supposed to have on the next day? #
  datagr$Size_nextday= with(datagr, rfishprod::applyVBGF (Lmeas = Size,
                                                          Lmax = MaxSizeTL,
                                                          Kmax = Kmax,
                                                          t = 1))
  
  # Estimating gross somatic growth (g) #
  datagr$somatic_growth = with(datagr, rfishprod::somaGain (a = a,
                                                            b = b,
                                                            Lmeas = Size,
                                                            Lmax = MaxSizeTL,
                                                            Kmax = Kmax,
                                                            t = 1))
  
  range(datagr$somatic_growth)
  
  # Applying stochastic mortality #
  datagr$mortality = rfishprod::applyMstoch(datagr$Md,t=1)
  
  
  #Alternatively, estimating per capita mass loss due to mortality #
  datagr$soma_loss = with(datagr, rfishprod::somaLoss (M = Md,
                                                       Lmeas = Size,
                                                       a = a,
                                                       b = b,
                                                       t = 1))
  
  
  datagr_prod = datagr |>
    #Calculating biomass turnover
    dplyr::mutate(W = a*(Size^b),
                  Biom = (W*Num),
                  Prod = ifelse(mortality == T, (somatic_growth * Num),0))
  
  return(datagr_prod)
  
}
