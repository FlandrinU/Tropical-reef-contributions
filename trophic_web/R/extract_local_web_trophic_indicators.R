################################################################################
##
## Extract local trophic web of each surveys from the global metaweb and assess 
##  trophic indicators
##
## extract_local_web_trophic_indicators.R
##
## 15/02/2023
##
## Ulysse Flandrin
##
################################################################################

# #-----------------Loading packages-------------------
# pkgs <- c("here", "parallel", "igraph", "NetIndices", "ggplot2", "tibble", "dplyr",
#           "gtools")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading and preping data-------------
#Metaweb
load(file = here::here("trophic_web","outputs", "final_metaweb.Rdata"))
# Presence/absence matrix
load(file=here::here("biodiversity", "outputs", "occurrence_matrix_sp_survey_01.Rdata"))
# relative biomass matrix
load(file=here::here("biodiversity", "outputs", "occurrence_matrix_sp_survey_relative_biomass.Rdata"))

## Add producers in PA matrix: they are omnipresent
PA_matrix_survey <- cbind(as.data.frame(surveys_sp_occ),
                        primary_producers = rep(1,nrow(surveys_sp_occ)),
                        secondary_producers = rep(1,nrow(surveys_sp_occ)))

  
##------------- Functions to extract trophic indicators (from Albouy et al. 2019) -------------

## Principal function P_A_data: presence absence matrix
#                     Lniche: Metaweb
#                     mc.cores: number of core for parrallel computing

get_indic_cells  <-  function(P_A_data = PA_matrix_survey, 
                              biomass_data = surveys_sp_pbiom,
                              Lniche = final_metaweb, 
                              mc.cores=1, bin_threshold= 0.6){ 
  
  pkg <- c("NetIndices","igraph","gtools")
  sapply(pkg,require,character.only = TRUE)  
  
  reseau_cell <-  parallel::mclapply(1:nrow(P_A_data),mc.cores=mc.cores,function(i){
    if((i%%100)==0) cat("i=",i,"\n")
    Names <- names(P_A_data[i,which(P_A_data[i,]>0)])
    Calc_indic_proba(x=Lniche[Names,Names], biomass_data, bin_threshold = bin_threshold,
                     local_web = i)
  })
  reseau_cell
} # get_indic_cells 



Calc_indic_proba  <- function(x=Lniche[Names,Names], biomass_data, bin_threshold=0.8, local_web = 1) {
  
  #cat("### Indicator proba calculation ###", "\n")
  Species <- nrow(x)
  Connectance_p <- sum(x)/Species^2
  
  #cat("### Binary calculation indicators ###","\n")
  bin_net <- round(x,4); bin_net[bin_net<bin_threshold] <- 0 
  bin_net[bin_net>=bin_threshold] <- 1
  binary_indic <- get_binary_indic(web=bin_net, S=Species, biomass_data, local_web=local_web)
  
  web=data.matrix(bin_net)
  #cat("### Path calculation indicators ###","\n")
  Path_stat <- get_path_stats(web=web)
  
  #cat("### Igraph indicator calculation ###", "\n")
  igraph_res <- Calc_indic_igraph(web=web)
  
  c(Species=Species,Connectance_p=Connectance_p,binary_indic,igraph_res,Path_stat)
  
}  # end of function Calc_ind



get_binary_indic <- function(web=mat, S=Species, biomass_data = surveys_sp_pbiom, local_web = 1){
  TL <- NetIndices::TrophInd(Flow =web,Tij = t(web))
  length_chain <- ceiling(round(max(TL[,1]),1))
  TL_moy  <- mean(TL[,1])
  
  sp <- rownames(TL)[-which(rownames(TL) %in% c("primary_producers", "secondary_producers"))]
  biom_weighted_mTL <- sum(TL[sp, "TL"] * biomass_data[local_web, sp])
  
  Omn_moy <- mean(TL[,2])
  if ( nrow(TL)>20){
    Plankt <- length(which(TL[,1] > 2.4 & TL[,1] < 3.7 )) / nrow(TL) #proportion of plantivores
    HTI <- length(which(TL[,1] > 4 )) / nrow(TL) #High Trophic level Indicator (Bourdeau et al. 2016)
    MTI <- mean(TL[which(TL[,1] > 3.25 ),1])
  }else{ Plankt <- HTI <- MTI <- NA}
  
  Link <- sum(web); Link_max <- nrow(web)^2
  Connectance <- Link/Link_max
  b <- log2(Link) / (log2(S)-1) #parameter of the power law between L and S #Carpentier et al. 2021
  
  web <- web[-grep("producers",rownames(web), fixed=T),-grep("producers", colnames(web), fixed=T)]
  if( length(web)>1){
    Nprey <- colSums(web)  ### Nombre de proie qu'a chaque espèce
    Npred <- rowSums(web)  ### Nombre de prédateur qu'a chaque espece
    S <- nrow(web)
  }else{ Npred <- Nprey <- NA ; S <- 1}
  Ntop  <- ((sum(Npred == 0 & Nprey>=1)) / S)  #proportion of top predator
  Nbas  <- ((sum(Npred >= 0 & Nprey<=1)) / S) #proportion of herbivorous+planktivorous
  Nint  <- 1 - (Ntop + Nbas)
  
  Nprey <- Nprey[Nprey!=0] #espèces avec des proies
  Npred <- Npred[Npred!=0] #espèces avec des prédateurs
  
  Vul <-  mean(Npred) ; Vulsd<- sd(Npred)
  Gen <- mean(Nprey); Gensd <- sd(Nprey)
  
  c(Link=Link,Link_max=Link_max,Connectance=Connectance, b_power_law=b, Ntop=Ntop,
    Nbas=Nbas,Nint=Nint,Vul=Vul,Vulsd=Vulsd,Gen=Gen,Gensd=Gensd,TL_moy=TL_moy, weighted_mTL = biom_weighted_mTL,
    length_chain=length_chain,Omn_moy=Omn_moy, Planktivores=Plankt, HTI=HTI, MTI=MTI)
  
} # end of get_binary_indic


get_path_stats <- function(web=Path_net){
  web <-  igraph::graph.adjacency(web,mode="directed", weighted=NULL)
  Averpth_length <- igraph::average.path.length(web) #mean shortest path to each vertices
  
  b <- table( igraph::shortest.paths(web,mode="out"))
  b <- b[-which(names(b)==0 | names(b)==Inf)]
  c <- rep(NA,5) ; names(c) <- seq(1,5,1)
  d <- b/ sum(b)
  for (i in 1:length(d)){ c[i] <- d[i]}
  names(c) <- paste("Shortest_path_",names(c),sep="")
  
  c(Mean_path_length=Averpth_length,c)
  
} # end of get_path_stats


Calc_indic_igraph<- function(web=mat){
  web <-  igraph::graph.adjacency(web,weighted=NULL)
  
  Mod <-  igraph::modularity( igraph::walktrap.community(web))
  Diam <-  igraph::diameter(web)
  
  Din <-  igraph::degree(web, mode=c("in"))
  Din_stat <- c(Din_mean=mean(Din),Din_min = min(Din),Din_max = max(Din))
  
  Dout <-  igraph::degree(web, mode=c("out"))
  Dout_stat <- c(Dout_mean=mean(Dout),Dout_min = min(Dout),Dout_max = max(Dout))
  
  #Bet <- betweenness(web)
  #Bet_stat <- c(Betweenness_mean=mean(Bet),Betweenness_min = min(Bet),Betweenness_max = max(Bet))
  
  Clos <-  igraph::closeness(web,normalized = F)
  Clos_stat <- c(Closeness_mean=mean(Clos, na.rm=T),Closeness_min = min(Clos, na.rm=T),Closeness_max = max(Clos, na.rm=T))
  
  Nb_artpt <- length(articulation.points(web))  
  Trans <-   igraph::transitivity(web) # clustering coef
  
  Cor <-  igraph::graph.coreness(web)
  Cor_stat <- c(Coreness_mean=mean(Cor, na.rm=T),Coreness_min=min(Cor, na.rm=T),Coreness_max=max(Cor, na.rm=T))
  
  #Centrality <-  evcent(web)$vector
  #c(Modularity=Mod, Diameter= Diam, Din_stat,Dout_stat,Bet_stat,Clos_stat,Nb_articulate_point = Nb_artpt,Transitivity=Trans,Cor_stat) 
  c(Modularity=Mod,Diameter= Diam, Din_stat,Dout_stat,Clos_stat,Nb_articulate_point = Nb_artpt,Transitivity=Trans,Cor_stat) 
  
} # end of function Calc_indic_igraph



##------------- Calculation of network indicators for each survey -------------
trophic_indicators <- get_indic_cells(P_A_data= PA_matrix_survey ,
                                      biomass_data = surveys_sp_pbiom,
                                      Lniche= final_metaweb,
                                      mc.cores=parallel::detectCores()-5, bin_threshold=0.6) 

trophic_indicators_survey <- do.call(rbind, trophic_indicators)
trophic_indicators_survey <- as.data.frame(cbind(SurveyID = rownames(PA_matrix_survey), 
                                                 trophic_indicators_survey)) |>
  dplyr::mutate(across(.cols = everything(), .fns = as.numeric, .names = "{.col}" ))

save(trophic_indicators_survey, 
     file= here::here("trophic_web", "outputs", "trophic_indicators_survey.Rdata"))
