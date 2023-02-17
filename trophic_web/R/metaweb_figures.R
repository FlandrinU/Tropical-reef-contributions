################################################################################
###                            Analysis Metaweb                              ###
###                            Ulysse Flandrin                               ###
###                               14/02/23                                   ###
################################################################################

### Loading of library
pkgs <- c("here", "parallel", "igraph", "NetIndices", "ggplot2", "ggsignif",
          "igraph", "graphlayouts", "ggraph")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

#data
load(here::here("data","metadata_surveys.Rdata"))
load(file = here::here("trophic_web", "outputs", "data_species_trophic_web.Rdata"))
load(file = here::here("trophic_web", "outputs", "final_metaweb.Rdata"))
load(file=here::here("biodiversity", "outputs", "occurrence_matrix_sp_survey_01.Rdata"))

Traits <- data_species_trophic_web
MW <- data.frame(final_metaweb)
PA_matrix_site <- as.data.frame(surveys_sp_occ) |>
  tibble::rownames_to_column(var= "SurveyID") |>
  dplyr::left_join( dplyr::select(metadata_surveys, SiteCode, SurveyID)) |>
  dplyr::select(-SurveyID) |>
  dplyr::group_by( SiteCode) |>
  dplyr::summarise(across(.cols = everything(), .fns = max, .names = "{.col}")) |>
  tibble::column_to_rownames(var= "SiteCode") |>
  dplyr::bind_cols(primary_producers = rep(1, length(unique(metadata_surveys$SiteCode))),
                   secondary_producers = rep(1, length(unique(metadata_surveys$SiteCode))))


## Threshold influence
thres <- seq(0.1,0.975,0.025)

### obs metaweb
TLm <- c()
for(t in thres){
  bin_net <- round(MW, 4); bin_net[bin_net<t] <- 0 
  bin_net[bin_net>=t] <- 1
  
  Troph <- NetIndices::TrophInd(Flow =bin_net,Tij = t(bin_net))
  sp <- which(is.element(rownames(Troph), Traits$species))
  TL <- abs(Troph$TL[sp] - Traits$trophic_level[sp])
  TLm <- c(TLm, sum(TL, na.rm=T))
}


jpeg(here::here("trophic_web", "outputs", "Deviation from observed trophic levels.jpg"),
     quality = 100, width = 20, height = 15, units = "cm",   pointsize = 12, res = 400)
plot(TLm~thres, pch=1, xlab="Binary threshold", ylab= "Sum of (observed TL - inferred TL)" )
points(thres[which(TLm == min(TLm))], min(TLm), col= "red", pch=20, add=T)
abline(h=min(TLm), col="red", lty=2, lwd=0.5)
dev.off()


### Trophic level distribution in metaweb
bin_threshold=0.6
bin_net <- round(MW,4); bin_net[bin_net<bin_threshold] <- 0 
bin_net[bin_net>=bin_threshold] <- 1
Troph <- NetIndices::TrophInd(Flow =bin_net,Tij = t(bin_net))
jpeg(here::here("trophic_web", "outputs", "Hist_metaweb_trophic_level.jpg"), quality = 100, width = 10, height = 10, units = "cm",   pointsize = 12, res = 200)
hist(Troph$TL, main=paste0("Inferred trophic level, threshold = ", bin_threshold), freq=T, breaks=30)
dev.off()


## analyse of predicted trophic level
t=bin_threshold
bin_net <- round(MW, 4); bin_net[bin_net<t] <- 0 
bin_net[bin_net>=t] <- 1
Troph <- NetIndices::TrophInd(Flow =bin_net,Tij = t(bin_net))
sp_nb <- which(is.element(rownames(Troph), Traits$species))

# --> for all species
sp <- rownames(Troph)[sp_nb]
TL <- Troph[sp, "TL"] - Traits[which(Traits$species %in% sp), 'trophic_level']
jpeg(here::here("trophic_web", "outputs", "Histogramm_of_deviation_in_trophic_level.jpg"), quality = 100, width = 13,      height = 10, units = "cm",   pointsize = 12, res = 200)
hist(TL$trophic_level, xlab="TL inferred - TL obs ", main= paste("Ecart entre TL observé et TL calculé, seuil =", t ), breaks=30)
dev.off()

# --> per size class
df1 <- as.data.frame(cbind(Traits$common_length, rep("Observed TL", length(sp)), Traits$trophic_level , rep(NA, length(sp))))
colnames(df1) <- c( "Length", "dataset", "TL", "Category"); rownames(df1) <- sp
df2 <- as.data.frame(cbind(Traits$common_length, rep("Inferred TL", length(sp)), Troph[sp, 'TL'] , rep(NA, length(sp))))
colnames(df2) <- c( "Length", "dataset", "TL", "Category"); rownames(df2) <- sp

df <- rbind(df1,df2)
df$TL <- as.numeric(df$TL)
df$Length <- as.numeric(df$Length)

for (i in 1:nrow(df)){
  if (df$Length[i]<10){ df$Category[i] <- ("< 10cm")
  }else if (df$Length[i]<40){ df$Category[i] <- ("10-40")
  }else if (df$Length[i]<70){ df$Category[i] <- ("40-70")
  }else if (df$Length[i]<100){ df$Category[i] <- ("70-100")
  }else if (df$Length[i]>=100){ df$Category[i] <- ("> 100cm")
  }}

df$Category <- factor(df$Category, levels = c("< 10cm", "10-40", "40-70", "70-100", "> 100cm"))
library(ggplot2)
P <- ggplot(df) +
  aes(x = Category, y = TL, fill = dataset) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  labs(x = "Size class (cm)", y = "Trophic level", title = "Trophic level observed and inferred in each size class") +
  #geom_text(data = data.frame(Category = c("< 10cm", "10-40", "40-70", "70-100", "> 100cm"), TL = c(4.6,4, 4.7,4.8,5.5)), label = "***")+
  theme_light() +
  theme(
    plot.title = element_text(size = 14L,
                              face = "bold",
                              hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold") )
P
ggsave(P, width = 20, height = 15,
       filename=here::here("trophic_web", "outputs", "distribution of observed and inferred TL per size class.png"))

df_gap <- cbind(df1,TL); colnames(df_gap) <- c( "Length", "dataset", "TL observed", "Category", "TL difference"); rownames(df_gap) <- sp
df_gap$`TL difference` <- as.numeric(df_gap$`TL difference`)
df_gap$Length <- as.numeric(df_gap$Length)

for (i in 1:nrow(df_gap)){
  if (df_gap$Length[i]<10){ df_gap$Category[i] <- ("< 10cm")
  }else if (df_gap$Length[i]<40){ df_gap$Category[i] <- ("10-40")
  }else if (df_gap$Length[i]<70){ df_gap$Category[i] <- ("40-70")
  }else if (df_gap$Length[i]<100){ df_gap$Category[i] <- ("70-100")
  }else if (df_gap$Length[i]>=100){ df_gap$Category[i] <- ("> 100cm")
  }}

df_gap$Category <- factor(df_gap$Category, levels = c("< 10cm", "10-40", "40-70", "70-100", "> 100cm"))
P_gap <- ggplot(df_gap) +
  aes(x = Category, y = `TL difference`) +
  geom_boxplot(fill = "#919396") +
  labs(x = "Size class (cm)", y = "TL inferred - TL observed", title = "Difference between trophic level inferred and observed \n in each size class") +
  theme_light() +
  theme(
    plot.title = element_text(size = 14L,
                              face = "bold",
                              hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold")
  )
P_gap
ggsave(P_gap, width = 20, height = 15,
       filename=here::here("trophic_web", "outputs", "difference between observed and inferred TL per size class.png"))



###----------------metaweb observation-------------------
###sub metaweb
i=2
Names <- sample(rownames(MW)[1:(length(rownames(MW))-2)], 98)
Names <- c(Names, "primary_producers", "secondary_producers")
sub_MW <- MW[Names,Names]

bin_net <- round(sub_MW,4); bin_net[bin_net<bin_threshold] <- 0
bin_net[bin_net>=bin_threshold] <- 1
Troph <- NetIndices::TrophInd(Flow =bin_net,Tij = t(bin_net))

graph <- igraph::graph.adjacency(data.matrix(bin_net),weighted=TRUE)

layout.matrix<-matrix(nrow=length(V(graph)),ncol=2)  # Rows equal to the number of vertices
layout.matrix[,1]<-stats::runif(length(V(graph))) # randomly assign along x-axis
layout.matrix[,2] <- Troph$TL # y-axis value based on trophic level


Degree_in <- colSums(bin_net)
Degree <- colSums(bin_net) + rowSums(bin_net)
Trophic_Level <- as.character(round(Troph$TL,0))

P <- ggraph::ggraph(graph, layout = "stress")+
  ggraph::geom_edge_link(aes(edge_alpha = 0.1), edge_colour = "grey80", arrow=arrow(ends="last", angle=20, length=unit(0.15, "inches"), type="closed"), show.legend=F)+
  ggraph::geom_node_point(aes(size = Degree), col = "white")+ #enlève les pointe de flèche visible par transparence
  ggraph::geom_node_point(aes( fill = Trophic_Level, size = Degree, alpha=0), shape = 21, show.legend = c(fill=T, Degree_in =T, alpha=F)) +
  ggraph::geom_node_text(aes(label = name, size=20, fontface="italic"), family = "Helvetica-Narrow", repel=T, show.legend = F)+
  labs(size="Degree", fill= "Trophic level")+
  guides(fill = guide_legend(override.aes = list(size = 5)))+
  theme_graph()+
  theme(text=element_text(size=15, face="bold"),legend.position=c(0.96,0.05), legend.box="horizontal", 
        plot.background =element_rect(fill="white", colour="black", linewidth = 2),
        legend.box.background =element_rect(fill="white", colour="black", size=0.5))

ggsave(P, width = 20, height = 15,
       filename=here::here("trophic_web", "outputs", paste0("Trophic_web_", i, ".png")))


P_tree <- ggraph(graph, layout = layout.matrix)+
  geom_edge_link(aes(edge_alpha = 0.1), edge_colour = "grey66", arrow=arrow(ends="last", angle=20, length=unit(0.15, "inches"), type="closed"), show.legend=F)+
  geom_node_point(aes(fill = Trophic_Level, size = Degree_in), shape = 21) +
  geom_node_text(aes(label = name), family = "serif")+
  theme_graph()+
  theme()
ggsave(P_tree, width = 20, height = 15,
       filename=here::here("trophic_web", "outputs", paste0("Trophic_web_tree_", i, ".png")))

###local web
plot_local_web <- function(mat_PA = PA_matrix_site, MW = final_metaweb, site_code, bin_threshold ){
  Names <- names(mat_PA[site_code,which(mat_PA[site_code,]>0)])
  local_web = MW[Names,Names]
  
  bin_net <- round(local_web,4); bin_net[bin_net<bin_threshold] <- 0
  bin_net[bin_net>=bin_threshold] <- 1
  Troph <- NetIndices::TrophInd(Flow =bin_net,Tij = t(bin_net))
  
  graph <- igraph::graph.adjacency(data.matrix(bin_net),weighted=TRUE)
  
  set.seed(6)
  layout.matrix<-matrix(nrow=length(V(graph)),ncol=2)  # Rows equal to the number of vertices
  layout.matrix[,1]<-stats::runif(length(V(graph))) # randomly assign along x-axis
  layout.matrix[,2] <- Troph$TL # y-axis value based on trophic level
  
  Degree_in <- colSums(bin_net)
  Degree <- colSums(bin_net) + rowSums(bin_net)
  Trophic_Level <- as.character(round(Troph$TL,0))
  
  P_tree <- ggraph(graph, layout = layout.matrix)+
    geom_edge_link(aes(edge_alpha = 0.1), edge_colour = "grey66", arrow=arrow(ends="last", angle=20, length=unit(0.15, "inches"), type="closed"), show.legend=F)+
    geom_node_point(aes(fill = Trophic_Level, size = Degree_in), shape = 21) +
    geom_node_text(aes(label = name), family = "serif")+
    theme_graph()+
    theme()
  
  P_tree
}

site_code = c("RAJA26", "FP42", "QLD135", "ETP155", "CAR3", "LHI33", "CAN88", "GBR18", "USEC24")
for(i in site_code){
  P_tree <- plot_local_web(mat_PA = PA_matrix_site, MW = final_metaweb, site_code=i, bin_threshold = 0.6 )
  ggsave(P_tree, width = 15, height = 10,
         filename=here::here("trophic_web", "outputs", paste0("Local_web_", i, ".png")))
}

