################################################################################
##
## Construction of the metaweb, with size based probability of interaction,
##  corrected i=with trophic guild, and water position 
##
## metaweb_construction.R
##
## 14/02/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse", "GenSA", "gtools", "parallel", "car")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
#interaction data from barnes et al. 2008
barnes_data <- read.csv(here::here("trophic_web", "data","size_barnes2008_FF.csv"), stringsAsFactors = T)

#species traits
load(file = here::here("trophic_web", "outputs", "data_species_trophic_web.Rdata"))

##-------------function loading -------------
source(here::here('trophic_web', 'R', 'Model_GenSA.R'))

##-------------calibrate the niche model-------------

# We used a global interaction data set from Barnes, et al. 2008 to calibrate a 
# model of trophic interactions using log of the observed body size for predator 
# (Mpred) and prey (Mprey). This data set was composed of 34931 marine predator and prey 
# interactions from 27 locations covering a wide range of environmental conditions, 
# from the tropics to the poles, for 93 predator species with sizes ranging from
# 0.3 cm to 309.69 cm and 174 prey species with sizes from 4.16 µm to 122.66 cm. 
# Interactions were compiled from published literature and if predator- or prey- 
# length was not measured in the original study, the length was calculated using 
# length-mass relationships.

### preping interaction data ###
MPred <- log10(barnes_data$standardised_predator_length)
MPrey <- log10(barnes_data$si_prey_length)

# data description #
levels(barnes_data$predator)
summary(barnes_data$predator_length)
levels(barnes_data$prey)
summary(barnes_data$prey_length)
length( table(paste(barnes_data$predator, barnes_data$prey)) )
#

# filter data: remove duplicates
interaction <- paste(barnes_data$predator, barnes_data$prey,
                     barnes_data$standardised_predator_length, 
                     barnes_data$si_prey_length)
unique_interaction <- unique(interaction)
ind <- c()
for (i in 1:nrow(barnes_data)){
  observation <- paste(barnes_data$predator[i], barnes_data$prey[i], 
                       barnes_data$standardised_predator_length[i],
                       barnes_data$si_prey_length[i])
  if(is.element(observation, unique_interaction)){
    ind <- c(ind,i)
    unique_interaction <- unique_interaction[unique_interaction != observation]
  }
}
unique_barnes_data <- barnes_data[ind,]


# filter data: reduce bias of over represented couples
final_barnes_data <- unique_barnes_data

for (pred in levels(final_barnes_data$predator)) {
  df_pred <- final_barnes_data[final_barnes_data$predator == pred, ]
  preys <- unique(df_pred$prey)
  
  for (prey in preys) {
    df_prey_pred <- df_pred[df_pred$prey == prey,]
    if (nrow(df_prey_pred) > 50) {
      remove <- sample(df_prey_pred$X,nrow(df_pred[df_pred$prey == prey,])-50 )
      final_barnes_data <- final_barnes_data[
        -which(is.element(final_barnes_data$X, remove)),]
}}}

MPred <- log10(final_barnes_data$standardised_predator_length)
MPrey <- log10(final_barnes_data$si_prey_length)


# description of new data #
table(paste(final_barnes_data$predator, final_barnes_data$prey))
plot((MPrey~MPred))
lm_M <- lm(MPrey~MPred)
plot(lm_M,2)
plot(lm_M,3)
plot(lm_M,1)
summary(lm_M)
#


### set model parameters ###
lm_M <- lm(MPrey~MPred)
pars <- c(a0 = lm_M$coefficients[1],a1 = lm_M$coefficients[2],b0 = sd(lm_M$residuals),b1 = 0) 
#Starting paramaters for the calibration
# a0.(Intercept)       a1.MPred             b0             b1 
# 0.1541158            0.5217897      0.2685798      0.0000000

# Setting the boundaries for the algorithm of parameters estimation.
par_lo <- c(a0 = -10, a1 = 0, b0 = -10, b1 = -10)
par_hi <- c(a0 = 10, a1 = 10, b0 = 10, b1 = 10)

# calibration data
data <- data.frame(MPrey = MPrey, MPred = MPred)

### Maximum likelihood estimation
estim.pars = GenSA(par = pars, fn = model, lower = par_lo, upper= par_hi, 
                   control = list(verbose =TRUE, max.time = 1000, smooth=FALSE),
                   data = data) #Search for parameters maximizing the posteriori probability of these observed interactions

### Saving the model parameters
write.table(estim.pars$par, file = here::here("trophic_web", "outputs", "niche_model_parameters.txt") )



##------------- Metaweb making -------------
Size <- log10(data_species_trophic_web$common_length)
names(Size) <- data_species_trophic_web$species

# Expand all pairs of species 
MPred <- log10(final_barnes_data$standardised_predator_length)
MPrey <- log10(final_barnes_data$si_prey_length)

M <- expand.grid(Size,Size) 
MPrey <- M[,1]; MPred <- M[,2]

Names <- expand.grid(names(Size),names(Size))
PreyNames <- Names[,1];PredNames <- Names[,2]

rm(M);rm(Names)

# Compute interaction probability
pars <- read.table(here::here("trophic_web", "outputs", "niche_model_parameters.txt"))
pLM <- pLMFitted(MPrey,MPred,pars) # function extracted from Model_GenSA.R code

# Transform into an adgency matrix
initial_metaweb <- matrix(pLM, nr = length(Size), nc = length(Size), byrow = FALSE)
colnames(initial_metaweb)<-rownames(initial_metaweb)<- names(Size)

# Save the initial metaweb
save(initial_metaweb, file = here::here("trophic_web", "outputs", "initial_metaweb.Rdata"))



##------------- Metaweb correction -------------

## Herbivory correction
herb_corrected_MW <- initial_metaweb

sp_herb <- data_species_trophic_web$species[
  which(data_species_trophic_web$Diet == "herbivores_microvores_detritivores")]

for(herb in sp_herb){
  herb_corrected_MW[which(herb_corrected_MW[, herb] > 0), herb] <- 0
  }


## self interaction: no cannibalism in general case
diag(herb_corrected_MW) <- 0 


## Small fish correction
small_corrected_MW <- herb_corrected_MW
small_species <- data_species_trophic_web$species[which(data_species_trophic_web$common_length<10)]
small_corrected_MW[, small_species] <- 0


## Link the primary and secondary producers
# Primary producer box
PP_prey <- rep(0,nrow(small_corrected_MW))
names(PP_prey) <- rownames(small_corrected_MW)
PP_prey[sp_herb] <- 1
PP_predator <- rep(0,nrow(small_corrected_MW)+2)

# Secondary producer box
low_TL_species <- data_species_trophic_web$species[
  which(data_species_trophic_web$Diet %in% c("sessile_invertivores",
                                             "corallivores",
                                             "microinvertivores",
                                             "macroinvertivores",
                                             "crustacivores",
                                             "planktivores"))]
second_prod <- unique(c(low_TL_species, small_species)) #small fishes + low trophic level species
PS_prey <- rep(0,nrow(small_corrected_MW))
names(PS_prey) <- rownames(small_corrected_MW)
PS_prey[second_prod] <- 1
PS_predator <- c(rep(0,nrow(small_corrected_MW)),1,1)

producers_MW <- rbind(small_corrected_MW, primary_producers=PP_prey, secondary_producers=PS_prey)
producers_MW <- cbind(producers_MW ,primary_producers =PP_predator, secondary_producers=PS_predator)


##------------- Save corrected metaweb -------------
final_metaweb <- producers_MW
save(final_metaweb, file = here::here("trophic_web","outputs", "final_metaweb.Rdata"))
