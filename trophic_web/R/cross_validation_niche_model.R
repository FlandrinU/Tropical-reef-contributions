#################################################################################
#'
#'This script makes cross validations of the niche model. 
#' 
#'
#'@author Ulysse Flandrin, \email{ulysse.flandrin@@gmail.com}
#'
#' @date 2023/03/03
#'
################################################################################

# ###--------------------- Library ---------------------
# pkgs <- c("GenSA", "gtools", "parallel", "PresenceAbsence" )
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))
# 

###--------------------- functions --------------------- 
source(here::here('trophic_web', 'R', 'Model_GenSA.R'))

###---------------------  Boyce index ---------------------
# (Hirzel et al. 2006) by Blaise Petitpierre (09.08.2012)

#### internal function calculating predicted-to-expected ratio for each class-interval
boycei<-function(interval,obs,fit){
  
  fit.bin<-fit
  obs.bin<-obs
  fit.bin[fit[]>=interval[1]&fit[]<=interval[2]]<-"i";fit.bin[fit.bin!="i"]<-0
  obs.bin[obs[]>=interval[1]&obs[]<=interval[2]]<-"i";obs.bin[obs.bin!="i"]<-0
  
  pi<-length(which(obs.bin=="i"))/length(obs)
  ei<-length(which(fit.bin=="i"))/length(fit.bin)
  fi<-pi/ei
  
  return(fi)
}

####' Calculating Boyce index as in Hirzel et al. 2006
####' fit: vector containing the suitability values of the study area (or extent used for the validation)
####' obs: vector containing the suitability values of the validation points
####' nclass : number of classes or vector with classes threshold. 
####'          If nclass=0, Boyce index is calculated with a moving window (see next parameters)
####' windows.w : width of the moving window (by default 1/10 of the suitability range)
####' res : resolution of the moving window (by default 100 focals)
####' PEplot : if True, plot the predicted to expected ratio along the suitability class

boyce<-function(fit,obs,nclass=0,window.w="default",res=100,PEplot=T){
  if(window.w=="default"){window.w<-(max(fit)-min(fit))/10}
  interval<-c(min(fit),max(fit))
  mini<-interval[1];maxi<-interval[2]
  
  if(nclass==0){
    vec.mov<- seq(from = mini, to = maxi - window.w, by = (maxi-mini- window.w)/res)
    interval<-cbind(vec.mov,vec.mov+window.w)
  }else if (length(nclass)>1){
    vec.mov<-c(mini,nclass)
    interval<-cbind(vec.mov,c(vec.mov[-1],maxi))
  }else if(nclass>0 & length(nclass)<2){
    vec.mov<-seq(from = mini, to = maxi, by = (maxi-mini)/nclass)
  }
  
  f<-apply(interval,1,boycei,obs,fit)
  
  b<-cor(f[which(f!="NaN")],vec.mov[which(f!="NaN")],method="spearman")
  
  if (PEplot==T)plot((apply(interval[which(f!="NaN"),],1,sum)/2),f[which(f!="NaN")],
                     xlab="Habitat suitability",ylab="Predicted/Expected ratio")
  
  results<-round(b,3);names(results)<-c("B")
  
  return(results)
}


###---------------------  Data treatment --------------------- 
#interaction data from barnes et al. 2008
barnes_data <- read.csv(here::here("trophic_web", "data","size_barnes2008_FF.csv"), stringsAsFactors = T)

### preping interaction data ###
MPred <- log10(barnes_data$standardised_predator_length)
MPrey <- log10(barnes_data$si_prey_length)

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


###--------------------- Simulation loop ---------------------
accuracy_model <- mclapply (1:999,mc.cores=12, function(i) { 
  
  # Sample 80% data #
  n80 <- sample(1:length(MPred), round(length(MPred)*0.8))
  MPred80 <- MPred[n80] ; MPrey80 <- MPrey[n80]
  MPred20 <- MPred[-n80] ; MPrey20 <- MPrey[-n80]
  data80 <- data.frame(MPrey = MPrey80, MPred = MPred80)
  data20 <- data.frame(MPrey20 = MPrey20, MPred20 = MPred20)
  # save(data80, file="Modèle niche/Cross validation/Sample0.8_MPred_MPrey.rdata")
  # save(data20, file="Modèle niche/Cross validation/Sample0.2_MPred_MPrey.rdata")
  
  
  # Model calibration
  lm_M80 <- lm(MPrey80~MPred80)
  pars80 <- c(a0 = lm_M80$coefficients[1],a1 = lm_M80$coefficients[2],b0 = sd(lm_M80$residuals),b1 = 0) 

  # Setting the boundaries for the algorithm of parameters estimation.
  par_lo <- c(a0 = -10, a1 = 0, b0 = -10, b1 = -10)
  par_hi <- c(a0 = 10, a1 = 10, b0 = 10, b1 = 10)
  
  # Computation of bayesian model
  estim.pars = GenSA(par = pars80, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, max.time = 50, smooth=FALSE), data = data80) #Search for parameters maximizing the posteriori probability of these observed interactions
  #write.table(estim.pars$par,file="Modèle niche/Cross validation/pars_0.8.txt")
  # a0.(Intercept)     a1.MPred80             b0             b1 
  # -0.5661537      0.9656987      0.1361257      0.1271494 
  
  ###### Cross validation #######
  Pars <- data.frame(par=estim.pars$par)
  pLM <- pLMFitted(MPrey20,MPred20,Pars) # function extracted from Model_GenSA.R code
  df_cross <- cbind(data20, pLM)
  

  # Reconstruction of modelled and observed data
  obsvdata  <- data.frame(MPred20,MPrey20,inter=rep(1,length(MPrey20)))
  simuldata <- data.frame(MPred20,MPrey20,inter=pLM)
  data_valid <- cbind(ID=c(1:length(MPrey20)),obsvdata[,3],simuldata[,3])
  
  # Accuracy calculation
  CMX1 <- cmx(data_valid, na.rm=TRUE, threshold=0.5)
  boyce_index <- boyce(simuldata,obsvdata,PEplot=F)
  data.frame(Sens=sensitivity(CMX1, st.dev=FALSE) *100,Boyce=boyce_index)
  
}) # end of simulation loop

### Save the data

accuracy_model <- do.call(rbind,accuracy_model)
save(file= here::here("trophic_web", "outputs", "accuracy_model_boyce_index.Rdata"),accuracy_model)

summary(accuracy_model$Boyce)


###---------------------  Model evaluation ---------------------
#### Test on some species ###
pars = read.table( here::here("trophic_web", "outputs", "niche_model_parameters.txt") )
data = final_barnes_data
pred <- c("Gadus morhua", "Zeus faber")
name_file = "Eval_model_"

for (i in pred){
  jpeg(here::here("trophic_web", "outputs", "niche_model_evaluation", paste0(name_file,i,".jpg")), width = 20, height = 10, units = "cm", pointsize = 12, res = 400)
  get_fig_eval(data=data,target_sp=i, pars=pars)
  dev.off()
}

