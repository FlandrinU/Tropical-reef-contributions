## function to compute all biodiversity metrics ##
# Author: Sébastien Villéger (sebastien.villeger@cnrs.fr)

## Inputs:
# - asb_sp_weight: a matrix with weight (biomass, coverage) of species (columns) in assemblages (rows)
# - sp_functcoord: a matrix with species(rows) coordinates along axes (columns) of a multidimensional functional space
# - sp_phylotree: an object "phylo" with phylogenetic relationships between species
# NB: species in the 3 objects should be the same

# Output:
# a matrix with for each assemblage (row), indices of taxonomic, functional and phylogenetic diversity (columns) indices rounded to 4 decimals

computing_biodiv<-function( asb_sp_weight, 
                            sp_functcoord,
                            sp_phylotree = NULL
                            ) {
  

## parameters and variables for computing biodiv indices ####
asb_nm<-row.names(asb_sp_weight)

# relative weights of species in assemblages
asb_sp_relw<-asb_sp_weight/apply(asb_sp_weight,1,sum)

# occurrences of species in assemblages
asb_sp_occ<-asb_sp_relw
asb_sp_occ[which(asb_sp_occ>0)]<-1


## Taxonomic diversity ####

asb_biodiv<-matrix(NA, nrow=length(asb_nm), ncol=2, 
                    dimnames=list(asb_nm, c("Taxo_Ric", "Taxo_Ent") ) 
                   )

# richness
asb_biodiv[,"Taxo_Ric"]<-apply(asb_sp_occ, 1, sum)

# entropy based on the Hill numbers framework (q=1 ; exp(Shannon) )
  for(i in asb_nm){
    asb_biodiv[i,"Taxo_Ent"]<-entropart::Diversity(Ps=asb_sp_relw[i,], 
                                                   q=1,
                                                   Correction="None",
                                                   CheckArguments=FALSE
                                                   ) 
  } # end of i



## Functional diversity ####
if (! is.null(sp_functcoord) ) {
  
  # Euclidean distance between species in the multidimensional space
  sp_dist_funct<-dist(sp_functcoord)
  
  # matrix to store indices values
  fd<-matrix(NA, nrow=length(asb_nm), ncol=2, 
             dimnames=list(asb_nm, c("Func_Ric", "Func_Ent") ) 
             )
  
  # functional richness with Hill numbers ( q=0 on species occurrences)
  fd[,"Func_Ric"] <- mFD::alpha.fd.hill (asb_sp_w = asb_sp_occ,
                                          sp_dist = sp_dist_funct,
                                          q = 0,
                                          tau= "mean",
                                          details_returned = FALSE
  )
  # functional entropy with Hill numbers (q=1 on species relative weights)
  fd[,"Func_Ent"] <- mFD::alpha.fd.hill (asb_sp_w = asb_sp_relw,
                                          sp_dist = sp_dist_funct,
                                          q = 1,
                                          tau= "mean",
                                          details_returned = FALSE
                                             )
  
  
  # merging with other metrics
  asb_biodiv<-cbind(asb_biodiv, fd)
  
  
}# end of funct div



## Phylogenetic diversity ####

if (! is.null(sp_phylotree) ) {
  
  # richness and entropy based on the Hill numbers framework
  pd<-matrix(NA, nrow=length(asb_nm), ncol=2, 
             dimnames=list(asb_nm, c("Phyl_Ric", "Phyl_Ent") ) 
             )
  
  # loop on assemblages          
  for(i in asb_nm)
  {
    # Phylo richness based on species occurrences and q=0
    pd[i,"Phyl_Ric"]<-entropart::ChaoPD(Ps=asb_sp_occ[i,]/sum(asb_sp_occ[i,]), 
                                        q=0, 
                                        PhyloTree=sp_phylotree, 
                                        Normalize=TRUE, 
                                        CheckArguments=FALSE)
    
    # Phylo entropy based on species weights and q=1
    pd[i,"Phyl_Ent"]<-entropart::ChaoPD(Ps=asb_sp_relw[i,], 
                                        q=1, 
                                        PhyloTree=sp_phylotree, 
                                        Normalize=TRUE, 
                                        CheckArguments=FALSE)
  } # end of i
  
  # merging with other metrics
  asb_biodiv<-cbind(asb_biodiv, pd)
  
}# end of phylo div

# rounding and returning
return(round(asb_biodiv,4) )

}# end of function
# 
# ###########################
# run_example<-FALSE
# if (run_example) {
#   # fake dataset
#   
#   ## creating fake datasets
#   
#   # parameters
#   nbasb<-6
#   nbsp<-50
#   nbdim<-4
#   
#   # assemblages with contrasted species richness and abundance evenness
#   asb_sp_weight<-matrix( 0, nbasb, nbsp, 
#                         dimnames=list( paste0("asb", 1:nbasb), paste0("sp",1:nbsp) )
#   )
#   asb_sp_weight[1,1:50]<-10 # 50 species, even abundance
#   asb_sp_weight[2,1:50]<-c(51,rep(1,49)) # 50 species, 1 over dominating species
#   asb_sp_weight[3,1:50]<-runif(50,1,20) # 50 species, random abundance
#   asb_sp_weight[4,1:10]<-10 # 10 species, even abundance
#   asb_sp_weight[5,1:10]<-c(51,rep(1,9)) # 10 species, 1 over dominating species
#   asb_sp_weight[6,1:10]<-runif(10,1,20) #  10 species, random abundance
#   
#   asb_sp_weight
#   
#   # species coordinates in a 4-dim space from a normal law (mean = 0, variance = 2)
#   sp_functcoord<-matrix( rnorm(nbsp*nbdim,0,2), nbsp, nbdim, 
#                     dimnames=list( paste0("sp", 1:nbsp), paste0("pc",1:nbdim) )
#   )
#   range(sp_functcoord)
#   
#   
#   # random phylogeny
#   sp_phylotree<-ape::rtree(nbsp, tip.label=paste0("sp",1:nbsp))
#   # plot.phylo(sp_phylotree)
#   
#   computing_biodiv(asb_sp_weight=asb_sp_weight, 
#                    sp_functcoord=sp_functcoord, 
#                    sp_phylotree=sp_phylotree)
#   
# }# end of example


