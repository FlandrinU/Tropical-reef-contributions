############################################
#
# Bayesian model of trophic interactions
# based on body size, using observation of interactions 
# only (presence)
# 
# Dominique Gravel
# June 25th, 2014
#
############################################


######################
# The model
model = function(pars,data) {
  MPred = data$MPred #list of log10 length of observed predator
	MPrey = data$MPrey
	a0 = pars[1]
	a1 = pars[2]
	b0 = pars[3]
	b1 = pars[4]
	meanMprey = mean(unique(MPrey))
	sdMprey = sd(unique(MPrey))
	

	# Optimum and range
	o = a0 + a1*MPred 
	r = b0 + b1*MPred
		
	# Compute the conditional
	pLM = exp(-(o-MPrey)^2/2/r^2)

	# Compute the marginal
	pM = dnorm(x=MPrey,mean=meanMprey,sd=sdMprey)
	
	#  Integrate the denominator
	pL = r/(r^2+sdMprey^2)^0.5*exp(-(o-meanMprey)^2/2/(r^2+sdMprey^2))	
	
	# Compute the posterior probability
	pML = pLM*pM/pL
	
	pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues

	return(-sum(log(pML)))		
}

######################

pLMFitted = function(MPrey,MPred,Pars) {
	with(Pars, {
        a0 = Pars[1,]; a1= Pars[2,]
        b0 = Pars[3,]; b1=Pars[4,]
        o = a0 + a1*MPred 
        #else if(isTemp == 1) o = a0 +(a1 + a2*Temp)*MPred + a3*Temp 
        r = b0 + b1*MPred 
        
        # Compute the conditional pLM for each predator
        exp(-(o-MPrey)^2/2/r^2)
	})
}




##############
### Model evaluation
require(ggplot2)

get_fig_eval <- function( data=data_end,target_sp = "Zeus faber",pars=pars) {

  par(mar = c(5,6,2,1),mfcol = c(1,2)) # plot parameters
  
  MPrey <- log10(data[,"si_prey_length"])
  MPred <- log10(data[,"standardised_predator_length"])
  
  target_preys = log10(data$si_prey_length[data$predator==target_sp])
  hist_target = hist(target_preys,plot = FALSE, breaks = 10)
  hist_target$counts = hist_target$counts/length(target_preys)
  hist_preys = hist(MPrey,plot = FALSE,breaks=20)
  hist_preys$counts = hist_preys$counts/length(MPrey)
  
  col_target <- alpha("mediumpurple1", alpha=0.3)
  col_prey <- alpha("darkgoldenrod3", alpha=0.3)
  
  plot(hist_target,col=col_target, xlim = range(MPrey),ylim = c(0,0.4),ylab = "Density",cex.axis = 1.25, cex.lab = 1.5,xlab = "log10 Body Size (cm)",main=paste("Predator:",target_sp)) #size of targeted preys
  plot(hist_preys,col=col_prey, add = TRUE) #size of all preys
  target_size = mean(log10(data$standardised_predator_length[data$predator == target_sp])) #taille moyenne du pred
  legend("topleft", legend=c("Observed preys", "All preys", "Predicted preys"), fill=c(col_target, col_prey, "black"), cex=1) 
  
  o = pars[1,] + pars[2,]*target_size
  r = pars[3,] + pars[4,]*target_size
  seqM = seq(min(MPrey)-0.5,max(MPrey),0.01)
  pLM = exp(-(o-seqM)^2/2/r^2)
  lines(seqM,0.2*pLM,lwd = 2) #interaction probability for this mean size of predator
  
  # Figure 1C
  seqX = seq(min(MPred),max(MPred),0.01)
  seqY = seq(min(MPrey),max(MPrey),0.01)
  XY = expand.grid(seqX,seqY)
  
  # Optimum and range
  o = pars[1,] + pars[2,]*XY[,1]
  r = pars[3,] + pars[4,]*XY[,1]
  
  # Compute the conditional
  pLM = exp(-(o-XY[,2])^2/2/r^2)
  Z = matrix(pLM,nr = length(seqX), nc = length(seqY))
  
  image(seqX,seqY,Z,xlab = "log10 Predator Body Size (cm)",ylab = "log10 Prey Body Size (cm)",col=heat.colors(10000),cex.axis = 1.25, cex.lab = 1.5)
  points(MPred,MPrey,pch = 19, cex = 0.1)

} # end of function get_fig_eval



### plot heatmap
require(ggplot2)

get_fig_eval_heatmap <- function( data=data_end,pars=pars) {
  
  MPrey <- log10(data[,"si_prey_length"])
  MPred <- log10(data[,"standardised_predator_length"])
  
  # Figure 1C
  seqX = seq(min(MPred),max(MPred),0.01)
  seqY = seq(min(MPrey),max(MPrey),0.01)
  XY = expand.grid(seqX,seqY)
  
  # Optimum and range
  o = pars[1,] + pars[2,]*XY[,1]
  r = pars[3,] + pars[4,]*XY[,1]
  
  # Compute the conditional
  pLM = exp(-(o-XY[,2])^2/2/r^2)
  Z = matrix(pLM,nr = length(seqX), nc = length(seqY))
  #Z[which(Z<0.8)] <- 0
  
  image(seqX,seqY,Z,xlab = "log10( Predator Body Size ) (cm)",ylab = "log10( Prey Body Size ) (cm)",col=heat.colors(10000),cex.axis = 1.25, cex.lab = 1.5,
        main="Probability of interaction in allometric model \n (threshold = 0.8)")
  
  col_out_threshold <- alpha("grey", alpha=0.5)
  Z[which(Z>0.8)]<- NA
  image(seqX,seqY,Z,col=col_out_threshold,cex.axis = 1.25, cex.lab = 1.5, add=T)
  
  points(MPred,MPrey,pch = 19, cex = 0.2)

} # end of function get_fig_eval

