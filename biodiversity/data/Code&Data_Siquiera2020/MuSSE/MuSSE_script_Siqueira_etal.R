#------------- Siqueira et al -------------#
#------- Reef fish trophic evolution ------#
#----------------- MuSSE ------------------#

#----- Loading packages -----#
library(ape)  # Version 5.3
library(phytools)  # Version 0.6.99
library(diversitree)  # Version 0.9.11
library(scales)  # Version 1.0
library(HDInterval)  # Version 0.2.0
library(coda)  # Version 0.19.3


#----- Setting WD -----#
wd <- "~/Desktop"   # Set the working directory
setwd(wd)


#----- Reading tree -----#
tr <- read.tree("TACT/Reef_fish_all_combined.trees")


#----- Reading data -----#
data <- read.csv("Data/Data_final_Siqueira_etal.csv", header=T, sep = ",", stringsAsFactors = T)

data.troph <- subset(data, data$Trophic_ID != "NA", select = c(Species, Trophic_ID))
rownames(data.troph) <- data.troph$Species

tr.all <- list()
for(i in 1:length(tr)){
  tr.all[[i]] <- keep.tip(tr[[i]], rownames(data.troph))
}

class(tr.all) <- "multiPhylo"


#----- Preparing data for analysis -----#
fish.troph <- as.numeric(droplevels(data.troph$Trophic_ID)); names(fish.troph) <- data.troph$Species


#----- Loop to run MuSSE in each tree and get the results -----#
# This will take very long if you want to run all 100 trees at once
# I suggest running the loop in parallel (e.g. for (i in 1:5)) and combining the results from the samp object below (this is the approach I used within a HPC environment)
# Alternatively, if you just want to check the results from a couple of trees, run the the loop with for (i in 1:2) (still takes quite a while) and plot;
# Results do not seem to change much between trees, so this approach should be enough to replicate the general pattern


samp <- data.frame(lambda1 = NA, lambda2 = NA, lambda3 = NA, lambda4 = NA, lambda5 = NA, lambda6 = NA, 
                   mu1 = NA, mu2 = NA, mu3 = NA, mu4 = NA, mu5 = NA, mu6 = NA, p = NA)

for (i in 1:length(tr.all)){
  troph.musse=make.musse(tr.all[[i]],fish.troph,length(unique(fish.troph)), sampling.f=0.92)
  troph.start=starting.point.musse(tr.all[[i]],length(unique(fish.troph)))
  
  #----- Fitting the model -----#
  troph.fituncons=find.mle(troph.musse,troph.start,control=list(maxit=100000))
  
  #----- MCMC -----#
  
  prior <- make.prior.exponential(1/2)
  prelim <- diversitree::mcmc(troph.musse, coef(troph.fituncons, full=TRUE), nsteps=100, prior=prior, w=1, print.every=10)
  
  w <- diff(sapply(prelim[2:(ncol(prelim)-1)], quantile, c(0.05, 0.95)))
  samples <- diversitree::mcmc(troph.musse, coef(troph.fituncons, full=TRUE), nsteps=2000, prior = prior, w=w, print.every=50)
  
  # Removing burnin
  samples <- samples[201:2000,]
  
  sam <- samples[,c(2:13,44)]
  
  samp <- rbind(samp, sam)
  
}

samp_final <- samp[-1,]


#----- Extracting rates -----#

#Diversification rates
dr_GC <- samp_final$lambda1-samp_final$mu1
dr_HD <- samp_final$lambda2-samp_final$mu2
dr_MI <- samp_final$lambda3-samp_final$mu3
dr_OM <- samp_final$lambda4-samp_final$mu4
dr_PK <- samp_final$lambda5-samp_final$mu5
dr_SI <- samp_final$lambda6-samp_final$mu6

dr_all <- cbind(dr_GC,dr_MI,dr_OM,dr_PK,dr_SI,dr_HD)

#Speciation rates
sr_GC <- samp_final$lambda1
sr_HD <- samp_final$lambda2
sr_MI <- samp_final$lambda3
sr_OM <- samp_final$lambda4
sr_PK <- samp_final$lambda5
sr_SI <- samp_final$lambda6

sr_all <- cbind(sr_GC,sr_MI,sr_OM,sr_PK,sr_SI,sr_HD)

#Extinction rates
xr_GC <- samp_final$mu1
xr_HD <- samp_final$mu2
xr_MI <- samp_final$mu3
xr_OM <- samp_final$mu4
xr_PK <- samp_final$mu5
xr_SI <- samp_final$mu6

xr_all <- cbind(xr_GC,xr_MI,xr_OM,xr_PK,xr_SI,xr_HD)


#----- Plotting results -----#
mode <- function(s) {
  d <- density(s)
  d$x[which.max(d$y)]
}


layout(matrix(c(1,1,1,2,3,4,5,6,7), 3, 3, byrow = TRUE), heights = c(1.5,1,1,1,1,1))

colours <- c("#BA3F1D","#F6AE2D","#99A1A6","#173F5F","#88498F","#388659")

dt = density(dr_all[,1])

par(mar = c(1.5, 0, 0, 0))
plot(dt, bty = "n", xlim = c(0.03,0.15), ylim = c(-15,190), col = "Black",
     xlab="", ylab="", yaxt="n", xaxt="n", main="", lwd=1, zero.line = F)
axis(side = 1, at = seq(0.03,0.15,0.04), lwd = 1.5, cex.axis = 1.2, tcl = -0.3, mgp = c(3, 0.5, 0))
mtext(side = 1, text = "Net Diversification", line = 1.8, cex = 1)


x_coord <- c(hdi(dr_all[,1])["lower"], dt$x[dt$x>=hdi(dr_all[,1])["lower"] & dt$x<=hdi(dr_all[,1])["upper"]], hdi(dr_all[,1])["upper"])
y_coord <- c(0,  dt$y[dt$x>=hdi(dr_all[,1])["lower"] & dt$x<=hdi(dr_all[,1])["upper"]], 0)

polygon(x_coord, y_coord, col=alpha("#BA3F1D", 0.7), lty = 0)
z <- seq(0, 1, length.out=ncol(dr_all) + 2)[-1] * par("usr")[3]
arrows(hdi(dr_all[,1])["lower"], z[1], hdi(dr_all[,1])["upper"], z[1], code=0, angle=90, col="#BA3F1D",lwd = 2)
points(mode(dr_all[,1]), z[1], col = "#BA3F1D", pch = 19, cex = 1)

for (i in 2:ncol(dr_all)){
  dt <- density(dr_all[,i])
  lines(dt, lwd=1, col="Black")
  x_coord <- c(hdi(dr_all[,i])["lower"], dt$x[dt$x>=hdi(dr_all[,i])["lower"] & dt$x<=hdi(dr_all[,i])["upper"]], hdi(dr_all[,i])["upper"])
  y_coord <- c(0,  dt$y[dt$x>=hdi(dr_all[,i])["lower"] & dt$x<=hdi(dr_all[,i])["upper"]], 0)
  polygon(x_coord, y_coord, col=alpha(colours[i], 0.7), lty = 0)
  arrows(hdi(dr_all[,i])["lower"], z[i], hdi(dr_all[,i])["upper"], z[i], code=0, angle=90, col=colours[i],lwd = 2)
  points(mode(dr_all[,i]), z[i], col = colours[i], pch = 19, cex = 1)
}

colours <- c("#BA3F1D","#F6AE2D","#99A1A6","#173F5F","#88498F","#388659")

for(i in 1:ncol(sr_all)){
  dt.s <- density(sr_all[,i])
  dt.x <- density(xr_all[,i])
  par(mar = c(3, 1.5, 0, 1.5))
  plot(dt.x, bty = "n", xlim = c(0,0.18), ylim = c(-4,190), col = colours[i],
       xlab="", ylab="", main="", lwd=1.5, yaxt="n", xaxt="n", zero.line = F)
  axis(side = 1, at = seq(0,0.18,0.06), lwd = 1.5, cex.axis = 1.2, tcl = -0.3, mgp = c(3, 0.5, 0))  
  
  polygon(dt.s, col=colours[i], lty = 0)
  
  arrows(hdi(sr_all[,i])["lower"], z[1], hdi(sr_all[,i])["upper"], z[1], code=0, angle=90, col=colours[i],lwd = 1.8)
  points(mode(sr_all[,i]), z[1], col = colours[i], pch = 19, cex = 0.8)
  
  arrows(hdi(xr_all[,i])["lower"], z[2], hdi(xr_all[,i])["upper"], z[2], code=0, angle=90, col=colours[i],lwd = 1.8)
  points(mode(xr_all[,i]), z[2], col = colours[i], pch = 19, cex = 0.8)
  if(i == 5){
    mtext(side = 1, text = "Extinction & Speciation", line = 1.8, cex = 1) 
  }
}


##########################################
