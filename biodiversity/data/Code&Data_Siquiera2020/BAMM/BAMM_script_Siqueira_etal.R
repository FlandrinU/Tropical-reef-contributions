#------------- Siqueira et al -------------#
#------- Reef fish trophic evolution ------#
#--------------- BAMM-rates ---------------#

#----- Loading packages -----#
library(ape)  # Version 5.3
library(BAMMtools)  # Version 2.1.6
library(coda)  # Version 0.19.3


#----- Setting WD -----#
wd <- "~/Desktop"   # Set the working directory
setwd(wd)

#----- Reading data -----#
data <- read.csv("Data/Data_final_Siqueira_etal.csv", header=T, sep = ",", stringsAsFactors = T)
rownames(data) <- data$Species


#----- Loop to get the results from 100 runs -----#

df.lambda <- as.data.frame(matrix(NA, ncol = 100, nrow = nrow(data)))
rownames(df.lambda) <- rownames(data)

df.mu <- as.data.frame(matrix(NA, ncol = 100, nrow = nrow(data)))
rownames(df.mu) <- rownames(data)

for(i in 1:100){
  tree <- read.tree(paste("BAMM/bamm",i,"/Reef_fish_all.tacted.newick.tre", sep = ""))
  
  events <- read.csv(paste("BAMM/bamm",i,"/event_data.txt", sep = ""))
  
  ed <- getEventData(tree, events, burnin = 0.1)
  
  tip.rates <- getTipRates(ed)
  
  lambda <- tip.rates$lambda.avg
  lambda <- lambda[rownames(data)]
  
  df.lambda[,i] <- lambda
  
  mu <- tip.rates$mu.avg
  mu <- mu[rownames(data)]
  
  df.mu[,i] <- mu
  
  cat (paste ("Tree", i), "\n")
}

df.lambda$Median <- apply(df.lambda, 1, median, na.rm = T)

df.mu$Median <- apply(df.mu, 1, median, na.rm = T)

median_df <- data.frame(lambda_median = df.lambda$Median, mu_median = df.mu$Median, net_div_median = df.lambda$Median - df.mu$Median,
                       row.names = rownames(data))

#write.csv(median_df, "Rates.csv")

#----- Checking for convergence -----#

mcmc_df <- data.frame(logLik = NA, N_shifts = NA)

ess <- vector()

for(i in 1:100){
  mcmc <- read.csv(paste("BAMM/bamm",i,"/mcmc_out.txt", sep = ""), header=T) 
  
  burnstart <- floor(0.1 * nrow(mcmc))
  mcmcpost <- mcmc[burnstart:nrow(mcmc), ]
  
  ess[i] <- effectiveSize(mcmcpost$logLik)
  
  mcmc_df <- rbind(mcmc_df, mcmcpost[,c("logLik", "N_shifts")])
}

mcmc_df <- mcmc_df[-1,]

mean(ess) 
