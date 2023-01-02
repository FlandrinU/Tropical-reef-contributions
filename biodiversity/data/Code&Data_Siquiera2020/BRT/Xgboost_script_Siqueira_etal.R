#------------- Siqueira et al -------------#
#------- Reef fish trophic evolution ------#
#---------------- xgboost -----------------#

# This Script was adapted from Morais and Bellwood (2018) Global drivers of reef fish growth. Fish and Fisheries, 19:874–889.

library (plyr)  # Version 1.8.4
library (tidyverse)  # Version 1.2.1
library (pdp)  # Version 0.7.0
library (xgboost)  # Version 0.82.1
library (ape)  # Version 5.3
library (phytools)  # Version 0.6.99
library (scales)  # Version 1.0


#-----------------#
### Loading data ##
#-----------------#

wd <- "~/Desktop"   # Set the working directory
setwd(wd)

db <- read.csv ("Data/Data_final_Siqueira_etal.csv", h = T, stringsAsF = F)
db <- db[complete.cases(db),]

db <- db %>% mutate (Trophic = factor (Trophic_ID, c ("GC", "HD", "MI", "OM", "PK", "SI")),
                     Activity = factor (Activity, c ("Both", "Day", "Night")),
                     Position = factor (Position, c ("Benthic", "Benthopelagic", "Pelagic")),
                     Ocean = factor (Ocean, c("IP", "Atl", "Both")))

rownames(db) <- db$Species

db$net_div_median <- ifelse(db$net_div_median < 0, 0, db$net_div_median)

#-----------------------#
### Setting the model  ##
#-----------------------#

fmod <- formula (~ Trophic + Size + SST + DistIAA + Activity + 
                   Position + Range + Pprod + Ocean)


#---------------------------------------#
#### Tuning XGBOOST: Systematic way #####
#---------------------------------------#

modmat <- model.matrix (fmod, db) [, -1]
dtrain <- xgb.DMatrix (data = modmat, label = db$net_div_median)

All_rmse.sist <- c ()
Param_group.sist <- list ()
Best_iter.sist <- c ()

set.seed (31)

pardb <- expand.grid (eta = seq (0.1, 0.4, by = 0.05),
			 gamma = seq (0, 0.2, by = 0.05),
			 max_depth = 4:7,
			 subsample = seq (0.5, 0.9, by = 0.1))


for (j in 1:nrow (pardb)) {

	params <- list (booster = 'gbtree',
					objective = 'reg:gamma',
					eta = pardb [j, 'eta'],
					gamma = pardb [j, 'gamma'],
					max_depth = pardb [j, 'max_depth'],
					subsample = pardb [j, 'subsample'])

	xgb.tune <- xgb.cv (params = params,
							  data = dtrain,
							  nrounds = 150,
							  nfold = 5,
							  metrics = list ('rmse'),
							  showsd = T,
							  stratified = T,
							  verbose = F,
							  early_stopping_rounds = 50,
							  maximize = F)
			min_rmse <- min (xgb.tune$evaluation_log$test_rmse_mean)
			All_rmse.sist <- append (All_rmse.sist, min_rmse)
			Param_group.sist [[j]] <- params
			Best_iter.sist <- append (Best_iter.sist, xgb.tune$best_iteration)
	
	cat (paste ("Tuning step", j), "\n")
}

(params.sist = Param_group.sist [[which.min (All_rmse.sist)]])
(best.iter <- Best_iter.sist [[which.min (All_rmse.sist)]])
min (All_rmse.sist)


#-----------------------------------#
#### Tuning XGBOOST: Random way #####
#-----------------------------------#


All_rmse.rand <- c ()
Param_group.rand <- list ()
Best_iter.rand <- c ()


for (i in 1:1000) {
	set.seed (i)
	params <- list (booster = 'gbtree',
						  objective = 'reg:gamma',
						  eta = runif (1,  0.28, 0.32),
						  gamma = runif (1, 0.03, 0.08),
						  max_depth = sample (6:8,1),
						  subsample = runif (1, 0.88, 0.92),
						  colsample_bytree = 1,
						  min_child_weight = 1)
	xgb.tune <- xgb.cv (params = params,
							  data = dtrain,
							  nrounds = 200,
							  nfold = 5,
							  metrics = list ('rmse'),
							  showsd = T,
							  stratified = T,
							  verbose = F,
							  early_stopping_rounds = 50,
							  maximize = F)
			min_rmse <- min (xgb.tune$evaluation_log$test_rmse_mean)
			All_rmse.rand <- append (All_rmse.rand, min_rmse)
			Param_group.rand [[i]] <- params
			Best_iter.rand <- append (Best_iter.rand, xgb.tune$best_iteration)
	
	cat (paste ("Tuning step", i), "\n")
}

(params.rand = Param_group.rand [[which.min (All_rmse.rand)]])
(best.iter.rand <- Best_iter.rand [[which.min (All_rmse.rand)]])
min (All_rmse.rand)

if (min (All_rmse.sist) < min (All_rmse.rand)) params = params.sist else params = params.rand
#saveRDS (params, 'xgboostparams.rds')
#params = readRDS ('xgboostparams.rds')

#------------------------------------------------------#
#### Now cross-validating predictions with XGBOOST #####
#------------------------------------------------------#

niter <- 1000
set.seed (31)

meanbias <- vector ()
medianbias <- vector ()
lowci <- vector ()
upci <- vector ()
R2 <- vector ()

for (i in 1:niter) {

samp <- sample (1:nrow (db), round (0.2*nrow(db)), replace = FALSE)

db.train <- db [-samp, ]
db.test <- db [samp, ]

modmat.train <- model.matrix (fmod, db.train) [, -1]
modmat.test <- model.matrix (fmod, db.test) [, -1]

xgbmod = xgboost (modmat.train, label = db.train$net_div_median,
						nrounds = best.iter, params = params, verbose = 0, print_every = 1000)
						
		pred <- predict (xgbmod, newdata = modmat.test, ntreelimit = best.iter)
		meas <- db.test$net_div_median
		
		bias <- meas - pred
		
		meanbias [i] <- mean (bias)
		medianbias [i] <- median (bias)
			q <- qt (.975, df = length (bias) -1)
			se <- sd (meas - pred) / sqrt (length (pred))
		lowci [i] <- mean (bias) - q*se
		upci [i] <- mean (bias) + q*se
		R2 [i] <- summary (lm (pred ~ meas))$r.squared

		cat (paste ("Model and prediction", i), "\n")

}


XGBoutput <- data.frame (medianbias = medianbias, meanbias = meanbias, lowcibias = lowci, upcibias = upci, R2 = R2)
#write.csv (XGBoutput, 'XGBoutput.csv', row.names = F)
#XGBoutput <- read.csv ('XGBoutput.csv', h = T)


#-------------------------------------#
#### Checking relevant indicators #####
#-------------------------------------#
ma <- mean (XGBoutput$meanbias)
ma
r <- mean (XGBoutput$R2)
round(r,2)


##############################################
########### PREDICTING VALUES ################
##############################################

#-----------------------------#
### Categorizing variables  ###
#-----------------------------#

grid <- expand.grid (Trophic = levels (db$Trophic), Size = 	seq (min (db$Size), 150, length = 20), SST = mean(db$SST), 
                     DistIAA = mean(db$DistIAA), Activity = levels (db$Activity), Position = levels (db$Position), 
                     Range = mean (db$Range), Pprod = mean (db$Pprod), Ocean = levels (db$Ocean)) 		

		
#----------------------------------#
### Creating grid of predictions ###
#----------------------------------#

modmatgrid <- model.matrix (fmod, data = grid) [, -1]


##----------------------------------##
#### Bootstrapping and predicting ####
##----------------------------------##

modmat <- model.matrix (fmod, db) [,-1]

niter = 1000
set.seed (31)

xgbmod <- list ()
rel.imp <- list ()
pred.grid <- as.data.frame (matrix (ncol = niter, nrow = nrow (grid), 
									dimnames = list (NULL, paste0 ('Boot',1:niter))))


for (i in 1:niter) {
	
	set.seed (i)
	xgbmod [[i]] <- xgboost (modmat, label = db$net_div_median,
									nrounds = best.iter, params = params, verbose = 0, print_every = 1000)
	rel.imp [[i]] <- xgb.importance (colnames (modmat), model = xgbmod[[i]])
	
	pred.grid [, i] <- predict (xgbmod [[i]], newdata = modmatgrid)
	
	cat (paste ("Bootstrapping the model, round", i), "\n")
	
}


#--------------------------------------------------------#
### Exploring the relative importance of the variables ###
#--------------------------------------------------------#

relimpdf <- do.call (rbind, rel.imp) %>% 
              mutate (FeatureProc = substr (Feature, 1,5)) %>%
              group_by (Feature, FeatureProc) %>%
                summarise (meanGain = mean (Gain),
                           uppGain = quantile(Gain)['75%'],
                           lowGain = quantile(Gain)['25%']) %>% 
                    group_by (FeatureProc) %>%
                      summarise (meanGain = sum (meanGain),
                                 uppGain = sum (uppGain),
                                 lowGain = sum (lowGain)) %>% ungroup %>%
                      arrange (desc(meanGain)) %>% as.data.frame


#----------------------------------#
### Getting the predicted values ###
#----------------------------------#

pred.grid.trim <- pred.grid[grid$Activity == "Day" & grid$Position == "Benthic" & grid$Ocean == "IP",]

pred.grid.trim$Trophic <- grid[grid$Activity == "Day" & grid$Position == "Benthic" & grid$Ocean == "IP","Trophic"]
pred.grid.trim$Size <- grid[grid$Activity == "Day" & grid$Position == "Benthic" & grid$Ocean == "IP","Size"]

Size_df <- pred.grid.trim %>% gather (key = 'Iteration', value = 'Diversification', Boot1:Boot1000)
troph_df <- Size_df[Size_df$Size > 24 & Size_df$Size < 25,] 

Size_df_plot <- Size_df %>% group_by(Trophic, Size) %>% 
  summarize(median  = median(Diversification), low = quantile(Diversification)["25%"], upp = quantile(Diversification)["75%"]) %>% as.data.frame()

#---------------------------#
### Plotting predictions  ###
#---------------------------#


layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = T), heights = c(1.4,1))
par (oma = c(1,1.5,1,1))

#Plot1
names.arg = as.character(rev(relimpdf$FeatureProc))
col.bars = c(rep("gainsboro", 7), "deepskyblue3", "deepskyblue3")

x <- barplot (sort(relimpdf$meanGain), names.arg = names.arg, xlab = "", xaxt="n", lwd = 1,
         col = col.bars, border = NA,las = 1, horiz = T, xlim = c(0,0.50), cex.names = 1.5)
axis(side = 1, at = seq(0,0.50,0.1), lwd = 1.7, cex.axis = 1.2, tcl = -0.3, mgp = c(3, 0.5, 0), labels = c("0","10","20","30","40","50"))
mtext(side = 1, text = "Relative importance (%)", line = 1.8, cex = 1)
lines(c(0.11, 0.11), c(-1,11), col = "deepskyblue4", lwd = 2, lty = 2)

for(i in 1:nrow(relimpdf)){
  lines(c(relimpdf$uppGain[i], relimpdf$lowGain[i]),
        c(rev(x)[i],rev(x)[i]), 
        col = "grey8", lwd = 3)
  points(x = relimpdf$meanGain[i],rev(x)[i] , cex = 1.5, bg = rev(col.bars)[i], pch = 21)
}


#Plot2
colours <- c("#BA3F1D","#388659","#F6AE2D","#99A1A6","#173F5F","#88498F")

plot(jitter(as.numeric(as.factor(troph_df$Trophic)), 0.8), troph_df$Diversification, ylim = c(0.05,0.15), 
     pch = 16,col=alpha(colours,0.08), xaxt = 'n', yaxt = 'n', ylab = "", xlab = "", tcl = -0.3, mgp = c(3, 0.5, 0), bty = "n")
axis (1, at = 1:6, labels = levels(troph_df$Trophic), tcl = -0.3, mgp = c(3, 0.5, 0), cex.axis = 1.2, lwd = 1.7)
axis (2, at = seq(0.05,0.16,0.03), labels = as.character(seq(0.05,0.16,0.03)), tcl = -0.3, mgp = c(3, 0.5, 0), cex.axis = 1.2, las = 1, lwd = 1.7)
mtext(side = 1, text = "Trophic identity", line = 1.8, cex = 1)
mtext(side = 2, text = "Predicted tip diversification rate", line = 2.7, cex = 1)

for(i in 1:length(unique(troph_df$Trophic))){
  lines(x = c(i,i), c(quantile(troph_df$Diversification[troph_df$Trophic == unique(troph_df$Trophic)[i]])["75%"],
                      quantile(troph_df$Diversification[troph_df$Trophic == unique(troph_df$Trophic)[i]])["25%"]), lwd = 7, col = "Black")
  points(x = i, summary(troph_df$Diversification[troph_df$Trophic == unique(troph_df$Trophic)[i]])["Median"], cex = 3, bg = colours[i], pch = 21)
}


#Plot 3
plot(1, type="n", xlab="", ylab="", xlim=c(0, 100), ylim=c(0.05, 0.18), tcl = -0.3, mgp = c(3, 0.5, 0), 
     las = 1, yaxt = 'n', xaxt = 'n', bty = "n", lwd = 1.7)
axis (1, at = seq(0,100,25), labels = as.character(seq(0,100,25)), tcl = -0.3, mgp = c(3, 0.5, 0), cex.axis = 1.2, las = 1, lwd = 1.7)
axis (2, at = seq(0.05,0.18,0.04), labels = as.character(seq(0.05,0.18,0.04)), tcl = -0.3, mgp = c(3, 0.5, 0), cex.axis = 1.2, las = 1, lwd = 1.7)
mtext(side = 1, text = "Max Body length", line = 1.8, cex = 1)


lines(c(10,10), c(0.05, 0.17), col = "grey8", lwd = 1.7, lty = 2)
lines(c(30,30), c(0.05, 0.17), col = "grey8", lwd = 1.7, lty = 2)

text(3, 0.16, "d", cex = 1.8)
text(20, 0.16, "e", cex = 1.8)
text(60, 0.16, "f", cex = 1.8)


for (i in 1:length(unique(Size_df_plot$Trophic))) {
  Size_df_group = Size_df_plot[Size_df_plot$Trophic == unique(Size_df_plot$Trophic)[i],]
  lines(Size_df_group$Size, Size_df_group$median, col = colours[i], lwd = 3)
  
  polygon(c(Size_df_group$Size, rev(Size_df_group$Size)), c(Size_df_group$upp, rev(Size_df_group$low)), 
          col = alpha(colours[i],0.2), border = NA)

}


#Plot 4
plot(1, type="n", xlab="", ylab="", xlim=c(1,6), ylim=c(-0.02, 0.03), tcl = -0.3, mgp = c(3, 0.5, 0), 
     cex.axis = 1.2, las = 1, yaxt = 'n', bty = "n", xaxt = "n")
axis (1, at = 1:6, labels = levels(troph_df$Trophic), tcl = -0.3, mgp = c(3, 0.5, 0), cex.axis = 1.2, lwd = 1.7)
axis (2, at = seq(-0.02, 0.03,0.01), labels = as.character(seq(-0.02, 0.03,0.01)), tcl = -0.3, mgp = c(3, 0.5, 0), cex.axis = 1.2, las = 1, lwd = 1.7)
mtext(side = 2, text = "Effect size", line = 2.7, cex = 1)

sub = subset(Size_df, Size_df$Size <= 10)
med = median(sub$Diversification)

lines(x = c(0,6.5), c(0,0), lwd = 2, col = "Grey",lty = 2)

for(i in 1:length(unique(troph_df$Trophic))){
  sub_troph = subset(sub, sub$Trophic == unique(troph_df$Trophic)[i])
  
  effect_size_troph = median(sub_troph$Diversification) - med
  upp_quant = quantile(sub_troph$Diversification)['75%'] - med
  low_quant = quantile(sub_troph$Diversification)['25%'] - med
  
  lines(x = c(i,i), c(upp_quant,low_quant), lwd = 4, col = "Black")
  points(x = i, effect_size_troph, cex = 2.5, bg = colours[i], pch = 21)
}


#Plot 5
plot(1, type="n", xlab="", ylab="", xlim=c(1,6), ylim=c(-0.02, 0.03), tcl = -0.3, mgp = c(3, 0.5, 0), 
     cex.axis = 1.2, las = 1, yaxt = 'n', bty = "n", xaxt = "n")
axis (1, at = 1:6, labels = levels(troph_df$Trophic), tcl = -0.3, mgp = c(3, 0.5, 0), cex.axis = 1.2, lwd = 1.7)
axis (2, at = seq(-0.02, 0.03,0.01), labels = as.character(seq(-0.02, 0.03,0.01)), tcl = -0.3, mgp = c(3, 0.5, 0), cex.axis = 1.2, las = 1, lwd = 1.7)
mtext(side = 1, text = "Trophic identity", line = 1.8, cex = 1)

sub = subset(Size_df, Size_df$Size > 10 & Size_df$Size <= 30)
med = median(sub$Diversification)

lines(x = c(0,6.5), c(0,0), lwd = 2, col = "Grey",lty = 2)

for(i in 1:length(unique(troph_df$Trophic))){
  sub_troph = subset(sub, sub$Trophic == unique(troph_df$Trophic)[i])
  
  effect_size_troph = median(sub_troph$Diversification) - med
  upp_quant = quantile(sub_troph$Diversification)['75%'] - med
  low_quant = quantile(sub_troph$Diversification)['25%'] - med
  
  lines(x = c(i,i), c(upp_quant,low_quant), lwd = 4, col = "Black")
  points(x = i, effect_size_troph, cex = 2.5, bg = colours[i], pch = 21)
}


#Plot 6
plot(1, type="n", xlab="", ylab="", xlim=c(1,6), ylim=c(-0.02, 0.1), tcl = -0.3, mgp = c(3, 0.5, 0), 
     cex.axis = 1.2, las = 1, yaxt = 'n', bty = "n", xaxt = "n")
axis (1, at = 1:6, labels = levels(troph_df$Trophic), tcl = -0.3, mgp = c(3, 0.5, 0), cex.axis = 1.2, lwd = 1.7)
axis (2, at = seq(-0.02, 0.1,0.04), labels = as.character(seq(-0.02, 0.1,0.04)), tcl = -0.3, mgp = c(3, 0.5, 0), cex.axis = 1.2, las = 1, lwd = 1.7)

sub = subset(Size_df, Size_df$Size > 30)
med = median(sub$Diversification)

lines(x = c(0,6.5), c(0,0), lwd = 2, col = "Grey",lty = 2)

for(i in 1:length(unique(troph_df$Trophic))){
  sub_troph = subset(sub, sub$Trophic == unique(troph_df$Trophic)[i])
  
  effect_size_troph = median(sub_troph$Diversification) - med
  upp_quant = quantile(sub_troph$Diversification)['75%'] - med
  low_quant = quantile(sub_troph$Diversification)['25%'] - med
  
  lines(x = c(i,i), c(upp_quant,low_quant), lwd = 4, col = "Black")
  points(x = i, effect_size_troph, cex = 2.5, bg = colours[i], pch = 21)
}


#########################################################################
