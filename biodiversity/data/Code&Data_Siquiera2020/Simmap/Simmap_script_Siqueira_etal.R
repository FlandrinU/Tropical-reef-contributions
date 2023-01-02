#------------- Siqueira et al -------------#
#------- Reef fish trophic evolution ------#
#---------------- simmaps -----------------#

#----- Loading packages -----#
library(ape)  # Version 5.3
library(phytools)  # Version 0.6.99
library(circlize)  # Version 0.4.6
library(scales)  # Version 1.0
library(diversitree)  # Version 0.9.11


#----- Setting WD -----#
wd <- "~/Desktop"   # Set the working directory
setwd(wd)


#----- Reading tree -----#
trees <- read.tree("TACT/Reef_fish_all_combined.trees")


#----- Reading data -----#
data <- read.csv("Data/Data_final_Siqueira_etal.csv", header=T, sep = ",", stringsAsFactors = T)

data.troph <- subset(data, data$Trophic_ID != "NA", select = c(Species, Trophic_ID))
rownames(data.troph) <- data.troph$Species

tr.all <- list()
for(i in 1:length(trees)){
  tr.all[[i]] = keep.tip(trees[[i]], rownames(data.troph))
}

class(tr.all) <- "multiPhylo"


#----- Preparing data for analysis -----#
fish.troph <- as.numeric(droplevels(data.troph$Trophic_ID)); names(fish.troph) <- data.troph$Species


#----- Sourcing simmap.musse function -----#
source(paste0(wd,"/simmap/make.simmap.musse.R"))


#----- Making simmap -----#
# This will take very long if you want to run all 100 trees at once
# I suggest running a number of trees in parallel (e.g. tr.all[1:5]) and combining the results from the Tot_matrix below (this is the approach I used within a HPC environment)
# Alternatively, if you just want to check the results from a couple of trees, run the analysis with tr.all[1:2] (still takes quite a while) and plot;
# Results do not seem to change much between trees, so this approach should be enough to replicate the general pattern

mtree <- make.simmap.musse(tr.all,fish.troph,nsim = 10, sampling.f = 0.92, Q = "musse") 


#----- Extracting mean transition values -----#
cols <- setNames(c("#BA3F1D","#388659","#F6AE2D","#99A1A6","#173F5F","#88498F"),sort(unique(fish.troph)))
plot(mtree[[1]],cols,ftype="i",lwd=1,mar=c(4.1,1.1,1.1,1.1), fsize = 0.001, type = "fan", part = 0.98)
trans_total <- sapply(mtree,markChanges, plot = F)
dev.off()

lis <- list()
for (x in 1:length(trans_total)){
  sub <- trans_total[[x]]
  Total_matrix <- matrix(0,nrow = length(unique(fish.troph)), ncol = length(unique(fish.troph)), dimnames = list(as.character(sort(unique(fish.troph))),as.character(sort(unique(fish.troph)))))
  for (i in 1:nrow(sub)){
    from = sub('-.*', '', rownames(sub)[i])
    to = sub('.*>', '', rownames(sub)[i])
    Total_matrix[from,to] = Total_matrix[from,to] + 1
  }
  lis[[x]] <- Total_matrix 
}

Tot_matrix <- Reduce('+', lis)/length(trans_total); rownames(Tot_matrix) <- levels(data.troph$Trophic_ID); colnames(Tot_matrix) <- levels(data.troph$Trophic_ID)


#----- Chord diagrams -----#
cols<-setNames(c("#BA3F1D","#388659","#F6AE2D","#99A1A6","#173F5F","#88498F"),levels(data.troph$Trophic_ID))

par(mfrow = c(2,3))
arr.col = data.frame(c("GC", "GC", "GC"), c("PK", "MI", "OM"),
                     c("black", "black", "black"))
chordDiagram(Tot_matrix, row.col = c(alpha("#BA3F1D", 0.8),alpha("#388659", 0.05),alpha("#F6AE2D", 0.05),alpha("#99A1A6", 0.05),
                                       alpha("#173F5F", 0.05),alpha("#88498F", 0.05)), annotationTrackHeight = c(0.03, 0.015),
             column.col = cols, grid.col = cols, annotationTrack = c("grid","axis"), scale = F,
             directional = 1, direction.type = "arrows",link.arr.col = arr.col,link.arr.length = 0.1, link.arr.type = "big.arrow")

arr.col = data.frame(c("MI", "MI", "MI", "MI"), c("GC","PK", "SI", "OM"),
                     c("black", "black", "black", "black"))
chordDiagram(Tot_matrix, row.col = c(alpha("#BA3F1D", 0.05),alpha("#388659", 0.05),alpha("#F6AE2D", 0.8),alpha("#99A1A6", 0.05),
                                       alpha("#173F5F", 0.05),alpha("#88498F", 0.05)), annotationTrackHeight = c(0.03, 0.015),
             column.col = cols, grid.col = cols, annotationTrack = c("grid","axis"), scale = F,
             directional = 1, direction.type = "arrows",link.arr.col = arr.col,link.arr.length = 0.1, link.arr.type = "big.arrow")

arr.col = data.frame(c("OM", "OM", "OM", "OM", "OM"), c("GC","PK", "SI", "MI", "HD"),
                     c("black", "black", "black", "black", "black"))
chordDiagram(Tot_matrix, row.col = c(alpha("#BA3F1D", 0.05),alpha("#388659", 0.05),alpha("#F6AE2D", 0.05),alpha("#99A1A6", 0.8),
                                       alpha("#173F5F", 0.05),alpha("#88498F", 0.05)), annotationTrackHeight = c(0.03, 0.015),
             column.col = cols, grid.col = cols, annotationTrack = c("grid","axis"), scale = F,
             directional = 1, direction.type = "arrows",link.arr.col = arr.col,link.arr.length = 0.1, link.arr.type = "big.arrow")

arr.col = data.frame(c("PK", "PK", "PK", "PK"), c("GC", "OM", "MI", "HD"),
                     c("black", "black", "black", "black"))
chordDiagram(Tot_matrix, row.col = c(alpha("#BA3F1D", 0.05),alpha("#388659", 0.05),alpha("#F6AE2D", 0.05),alpha("#99A1A6", 0.05),
                                       alpha("#173F5F", 0.8),alpha("#88498F", 0.05)), annotationTrackHeight = c(0.03, 0.015),
             column.col = cols, grid.col = cols, annotationTrack = c("grid","axis"), scale = F,
             directional = 1, direction.type = "arrows",link.arr.col = arr.col,link.arr.length = 0.1, link.arr.type = "big.arrow")

arr.col = data.frame(c("SI", "SI", "SI", "SI"), c("PK", "OM", "MI", "HD"),
                     c("black", "black", "black", "black"))
chordDiagram(Tot_matrix, row.col = c(alpha("#BA3F1D", 0.05),alpha("#388659", 0.05),alpha("#F6AE2D", 0.05),alpha("#99A1A6", 0.05),
                                       alpha("#173F5F", 0.05),alpha("#88498F", 0.8)), annotationTrackHeight = c(0.03, 0.015),
             column.col = cols, grid.col = cols, annotationTrack = c("grid","axis"), scale = F,
             directional = 1, direction.type = "arrows",link.arr.col = arr.col,link.arr.length = 0.1, link.arr.type = "big.arrow")

arr.col = data.frame(c("HD", "HD", "HD", "HD"), c("PK", "OM", "MI", "SI"),
                     c("black", "black", "black", "black"))
chordDiagram(Tot_matrix, row.col = c(alpha("#BA3F1D", 0.05),alpha("#388659", 0.8),alpha("#F6AE2D", 0.05),alpha("#99A1A6", 0.05),
                                       alpha("#173F5F", 0.05),alpha("#88498F", 0.05)), annotationTrackHeight = c(0.03, 0.015),
             column.col = cols, grid.col = cols, annotationTrack = c("grid","axis"), scale = F,
             directional = 1, direction.type = "arrows",link.arr.col = arr.col, link.arr.length = 0.1, link.arr.type = "big.arrow")


############################################################