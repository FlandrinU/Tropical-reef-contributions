################################################################################
###                            Analysis Metaweb                              ###
###                            Ulysse Flandrin                               ###
###                               14/02/23                                   ###
################################################################################
# 
# ### Loading of library
# pkgs <- c("here", "parallel", "igraph", "NetIndices", "ggplot2", "ggsignif",
#           "igraph", "graphlayouts", "ggraph")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
# ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

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
binary_web <- MW
binary_web["secondary_producers", "secondary_producers"] <- 0 #NetIndices can't process self-interaction
TLm <- c()

TL_diff <- parallel::mclapply( thres, mc.cores=15, function(t){
  binary_web[binary_web<t] <- 0 
  binary_web[binary_web>=t] <- 1
  
  Troph <- NetIndices::TrophInd(Flow =binary_web,Tij = t(binary_web))
  sp <- rownames(Troph)[is.element(rownames(Troph), Traits$species)]
  TL <- abs(Troph[sp, "TL"] - Traits[which(Traits$species == sp),"trophic_level"])
  
  #sum(TL, na.rm=T)
  data.frame(threshold = rep(t, nrow(TL)), TL_deviation = TL$trophic_level)
})

TL_difference <- do.call(rbind, TL_diff)
TL_diff_mean <- TL_difference |> 
  dplyr::group_by(threshold) |>
  dplyr::summarise(sum_deviation = sum(TL_deviation, na.rm=T),
                   mean_deviation = mean(TL_deviation, na.rm=T),
                   median_deviation = median(TL_deviation, na.rm=T))


jpeg(here::here("trophic_web", "outputs", "Deviation from observed trophic levels.png"),
     quality = 100, width = 30, height = 10, units = "cm",   pointsize = 12, res = 400)
  par(mfrow=c(1,3))
  plot(TL_diff_mean$sum_deviation ~ TL_diff_mean$threshold, pch=1, xlab="Binary threshold",
       ylab= "Sum of (fishbases's TL - inferred TL)" )
  points(thres[which(TL_diff_mean$sum_deviation == min(TL_diff_mean$sum_deviation))], 
         min(TL_diff_mean$sum_deviation), col= "red", pch=20, add=T)
  
  plot(TL_diff_mean$mean_deviation ~ TL_diff_mean$threshold, pch=1, xlab="Binary threshold",
       ylab= "Mean of (fishbases's TL - inferred TL)" )
  points(thres[which(TL_diff_mean$mean_deviation == min(TL_diff_mean$mean_deviation))], 
         min(TL_diff_mean$mean_deviation), col= "red", pch=20, add=T)
  
  plot(TL_diff_mean$median_deviation ~ TL_diff_mean$threshold, pch=1, xlab="Binary threshold",
       ylab= "Median of (fishbases's TL - inferred TL)" )
  points(thres[which(TL_diff_mean$median_deviation == min(TL_diff_mean$median_deviation))], 
         min(TL_diff_mean$median_deviation), col= "red", pch=20, add=T)
dev.off()


library(ggplot2)
ggplot(data = TL_difference, mapping = aes(x = as.factor(threshold), y = TL_deviation)) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  labs(x = "Threshold", y = "distribution of (fishbases's TL - inferred TL)", title = "") +
  theme_light()
ggsave(width = 20, height = 15, filename=here::here("trophic_web", "outputs", "distribution of (fishbase's - inferred TL) per threshold.png"))



## Threshold and trophic level
thres_figures = c(0.3,0.6,0.8)
rm(bin_threshold, binary_web)

parallel::mclapply( thres_figures, mc.cores=5, function(bin_threshold){
  rm(binary_web)
  
  binary_web <- MW
  binary_web["secondary_producers", "secondary_producers"] <- 0 #NetIndices can't process self-interaction

  # Trophic level distribution in metaweb
  binary_web[binary_web<bin_threshold] <- 0 
  binary_web[binary_web>=bin_threshold] <- 1
  Troph <- NetIndices::TrophInd(Flow =binary_web,Tij = t(binary_web))
  jpeg(here::here("trophic_web", "outputs", paste0("Hist_metaweb_trophic_level_threshold=", bin_threshold,".jpg")),
                  quality = 100, width = 10, height = 10, units = "cm",   pointsize = 12, res = 200)
  hist(Troph$TL, main=paste0("Inferred trophic level, threshold = ", bin_threshold), freq=T, breaks=30)
  dev.off()
  
  
  # --> for all species
  sp <- rownames(Troph)[is.element(rownames(Troph), Traits$species)]
  TL <- Troph[sp, "TL"] - Traits[which(Traits$species == sp),"trophic_level"]
  jpeg(here::here("trophic_web", "outputs", paste0("Hist_of_deviation_in_TL_threshold=", bin_threshold,".jpg")),
       quality = 100, width = 13,height = 10, units = "cm",   pointsize = 12, res = 200)
  hist(TL$trophic_level, xlab="TL inferred - TL from fishbase ",
       main= paste("Difference between inferred TL and TL from fishbases, treshold =", bin_threshold ), breaks=30)
  dev.off()
  
  # --> per size class
  df1 <- as.data.frame(cbind(Traits$common_length, rep("fishbase's TL", length(sp)), Traits$trophic_level , rep(NA, length(sp))))
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
    labs(x = "Size class (cm)", y = "Trophic level", title = "Trophic level from fishbase and inferred in each size class") +
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
         filename=here::here("trophic_web", "outputs", paste0("distribution of fishbase's and inferred TL per size class _ threshold=", bin_threshold, ".png")))
  
  df_gap <- cbind(df1,TL); colnames(df_gap) <- c( "Length", "dataset", "fishbase's TL", "Category", "TL difference"); rownames(df_gap) <- sp
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
    labs(x = "Size class (cm)", y = "TL inferred - fishbase's TL", 
         title = "Difference between inferred TL  and fishbase's TL \n in each size class") +
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
         filename=here::here("trophic_web", "outputs", paste0("Difference between fishbase's and inferred TL per size class _ threshold=", bin_threshold, ".png")))

}) #end of mclapply studying threshold influence on TL


###----------------metaweb observation-------------------
bin_threshold = 0.6 # threshold minimizing mean differences of TL and avoiding large outliers

###sub metaweb
for (i in c(1:3)){
  rm(Names, binary_web)
  Names <- sample(rownames(MW)[1:(length(rownames(MW))-2)], 98)
  Names <- c(Names, "primary_producers", "secondary_producers")
  binary_web <- MW[Names,Names]
  binary_web["secondary_producers", "secondary_producers"] <- 0 #NetIndices can't process self-interaction
  binary_web[binary_web<bin_threshold] <- 0 
  binary_web[binary_web>=bin_threshold] <- 1
  
  Troph <- NetIndices::TrophInd(Flow =binary_web,Tij = t(binary_web))
  
  graph <- igraph::graph.adjacency(data.matrix(binary_web),weighted=TRUE)
  
  layout.matrix <- matrix(nrow=length(igraph::V(graph)),ncol=2)  # Rows equal to the number of vertices
  layout.matrix[,1] <- stats::runif(length(igraph::V(graph))) # randomly assign along x-axis
  layout.matrix[,2] <- Troph$TL # y-axis value based on trophic level
  
  Degree_in <- colSums(binary_web)
  Degree <- colSums(binary_web) + rowSums(binary_web)
  Trophic_Level <- as.character(round(Troph$TL,0))
  
  library(ggplot2)
  P <- ggraph::ggraph(graph, layout = "stress")+
    ggraph::geom_edge_link(aes(edge_alpha = 0.1), edge_colour = "grey80", arrow=arrow(ends="last", angle=20, length=unit(0.15, "inches"), type="closed"), show.legend=F)+
    ggraph::geom_node_point(aes(size = Degree), col = "white")+ #enlève les pointe de flèche visible par transparence
    ggraph::geom_node_point(aes( fill = Trophic_Level, size = Degree, alpha=0), shape = 21, show.legend = c(fill=T, Degree_in =T, alpha=F)) +
    ggraph::geom_node_text(aes(label = name, size=20, fontface="italic"), family = "Helvetica-Narrow", repel=T, show.legend = F)+
    labs(size="Degree", fill= "Trophic level")+
    guides(fill = guide_legend(override.aes = list(size = 5)))+
    ggraph::theme_graph()+
    theme(text=element_text(size=15, face="bold"),legend.position=c(0.96,0.05), legend.box="horizontal", 
          plot.background =element_rect(fill="white", colour="black", linewidth = 2),
          legend.box.background =element_rect(fill="white", colour="black", size=0.5))
  
  ggsave(P, width = 20, height = 15,
         filename=here::here("trophic_web", "outputs", "metaweb",  paste0("Trophic_web_", i, ".png")))
  
  
  P_tree <- ggraph::ggraph(graph, layout = layout.matrix)+
    ggraph::geom_edge_link(aes(edge_alpha = 0.1), edge_colour = "grey66", arrow=arrow(ends="last", angle=20, length=unit(0.15, "inches"), type="closed"), show.legend=F)+
    ggraph::geom_node_point(aes(fill = Trophic_Level, size = Degree_in), shape = 21) +
    ggraph::geom_node_text(aes(label = name), family = "serif")+
    ggraph::theme_graph()+
    theme()
  ggsave(P_tree, width = 20, height = 15,
         filename=here::here("trophic_web", "outputs", "metaweb", paste0("Trophic_web_tree_", i, ".png")))

} # end for(i in 1:3) about metaweb 



###----------------local webs observation-------------------
plot_local_web <- function(mat_PA = PA_matrix_site, MW = final_metaweb, site_code, bin_threshold ){
  Names <- names(mat_PA[site_code,which(mat_PA[site_code,]>0)])
  binary_web <- MW[Names,Names]
  binary_web["secondary_producers", "secondary_producers"] <- 0 #NetIndices can't process self-interaction
  binary_web[binary_web<bin_threshold] <- 0 
  binary_web[binary_web>=bin_threshold] <- 1
  
  Troph <- NetIndices::TrophInd(Flow =binary_web,Tij = t(binary_web))
  graph <- igraph::graph.adjacency(data.matrix(binary_web),weighted=TRUE)
  
  layout.matrix<-matrix(nrow=length(igraph::V(graph)),ncol=2)  # Rows equal to the number of vertices
  layout.matrix[,1]<-stats::runif(length(igraph::V(graph))) # randomly assign along x-axis
  layout.matrix[,2] <- Troph$TL # y-axis value based on trophic level
  
  Degree_in <- colSums(binary_web)
  Degree <- colSums(binary_web) + rowSums(binary_web)
  Trophic_Level <- as.character(round(Troph$TL,0))
  
  library(ggplot2)
  P_tree <- ggraph::ggraph(graph, layout = layout.matrix)+
    ggraph::geom_edge_link(aes(edge_alpha = 0.1), edge_colour = "grey66", arrow=arrow(ends="last", angle=20, length=unit(0.15, "inches"), type="closed"), show.legend=F)+
    ggraph::geom_node_point(aes(fill = Trophic_Level, size = Degree_in), shape = 21) +
    ggraph::geom_node_text(aes(label = name), family = "serif")+
    ggraph::theme_graph()+
    theme()
  
  ggsave(P_tree, width = 15, height = 10,
         filename=here::here("trophic_web", "outputs", "local_webs", paste0("Local_web_", site_code, ".png")))
}


plot_TL_size_class_local_web <- function(mat_PA = PA_matrix_site, MW = final_metaweb, site_code, bin_threshold ){
  Names <- names(mat_PA[site_code,which(mat_PA[site_code,]>0)])
  binary_web <- MW[Names,Names]
  binary_web["secondary_producers", "secondary_producers"] <- 0 #NetIndices can't process self-interaction
  binary_web[binary_web<bin_threshold] <- 0 
  binary_web[binary_web>=bin_threshold] <- 1
  
  Troph <- NetIndices::TrophInd(Flow =binary_web,Tij = t(binary_web))
  sp <- rownames(Troph)[is.element(rownames(Troph), Traits$species)]
  TL <- Troph[sp, "TL"] - Traits[which(Traits$species %in% sp),"trophic_level"]

  df1 <- as.data.frame(cbind(Traits[which(Traits$species %in% sp),"common_length"],
                             rep("fishbase's TL", length(sp)),
                             Traits[which(Traits$species %in% sp),"trophic_level"] ,
                             rep(NA, length(sp))))
  colnames(df1) <- c( "Length", "dataset", "TL", "Category"); rownames(df1) <- sp
  df2 <- as.data.frame(cbind(Traits[which(Traits$species %in% sp),"common_length"],
                             rep("Inferred TL", length(sp)), Troph[sp, 'TL'] ,
                             rep(NA, length(sp))))
  colnames(df2) <- c( "Length", "dataset", "TL", "Category"); rownames(df2) <- sp
  
  df <- rbind(df1,df2)
  df_gap <- cbind(df,TL)
  colnames(df_gap) <- c( "Length", "dataset", "TL", "Category", "TL difference")
  
  df_gap$`TL difference` <- as.numeric(df_gap$`TL difference`)
  df_gap$TL <- as.numeric(df_gap$TL) ; df_gap$Length <- as.numeric(df_gap$Length)
  
  for (i in 1:nrow(df_gap)){
    if (df_gap$Length[i]<10){ df_gap$Category[i] <- ("< 10cm")
    }else if (df_gap$Length[i]<40){ df_gap$Category[i] <- ("10-40")
    }else if (df_gap$Length[i]<70){ df_gap$Category[i] <- ("40-70")
    }else if (df_gap$Length[i]<100){ df_gap$Category[i] <- ("70-100")
    }else if (df_gap$Length[i]>=100){ df_gap$Category[i] <- ("> 100cm")
    }}
  df_gap$Category <- factor(df_gap$Category, levels = c("< 10cm", "10-40", "40-70", "70-100", "> 100cm"))
  
  library(ggplot2)
  P <- ggplot(df_gap) +
    aes(x = Category, y = TL, fill = dataset) +
    geom_boxplot() +
    scale_fill_hue(direction = 1) +
    labs(x = "Size class (cm)", y = "Trophic level", title = "Trophic level from fishbase and inferred in each size class",
         subtitle = site_code) +
    theme_light()
  ggsave(P, width = 20, height = 15,
         filename=here::here("trophic_web", "outputs", "local_webs", paste0("distribution of fishbase's and inferred TL per size class _ " , site_code, ".png")))
  
    P_gap <- ggplot(df_gap) +
    aes(x = Category, y = `TL difference`) +
    geom_boxplot(fill = "#919396") +
    labs(x = "Size class (cm)", y = "TL inferred - fishbase's TL", 
         title = "Difference between inferred TL  and fishbase's TL \n in each size class",
         subtitle = site_code) +
    theme_light()
  ggsave(P_gap, width = 20, height = 15,
         filename=here::here("trophic_web", "outputs", "local_webs", paste0("Difference between fishbase's and inferred TL per size class _ ", site_code, ".png")))
}

site_code = c("RAJA26", "FP42", "QLD135", "ETP155", "CAR3", "LHI33", "CAN88", "GBR18", "USEC24")
parallel::mclapply(site_code, mc.cores = 10, function(i){
  plot_local_web(mat_PA = PA_matrix_site, MW = final_metaweb, site_code=i, bin_threshold )
  
  if(i %in% c("RAJA26", "FP42","LHI33" )){
    plot_TL_size_class_local_web(mat_PA = PA_matrix_site, MW = final_metaweb, site_code=i, bin_threshold )
  }
})

