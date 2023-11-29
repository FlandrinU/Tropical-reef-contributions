################################################################################
##
## From all contributions gathered in '1a_merge_contributions_surveys.R', this 
##  script reproduces the figure 1 of the paper Flandrin et al. 2023.
##
## 2a_make_fig_1.R
##
## 19/04/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
# pkgs <- c("here", "tidyverse","FactoMineR", "corrplot",
#           "factoextra", "ggplot2", "harrypotter", "grid", "gridExtra")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
library(ggplot2)

rm(list=ls())

##-------------loading data-------------
load(here::here("outputs","all_Contrib_site_log_transformed.Rdata"))

source(here::here("R", "elbow.R"))

color_grad = c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788")

##-------------Contributions in categories-------------
## Classify variables in Nature for Nature (NN) and Nature for People (NP)
grp_NN_NP <- as.factor(c(N_Recycling = "NN",
                         P_Recycling = "NN",
                         Taxonomic_Richness = "NN",
                         Functional_Entropy = "NN", 
                         Phylogenetic_Entropy = "NN",
                         Trait_Distinctiveness = "NN",
                         Evolutionary_Distinctiveness = "NN",
                         Herbivores_Biomass = "NN",
                         Invertivores_Biomass = "NN",
                         Piscivores_Biomass = "NN",
                         Endemism = "NN",
                         Elasmobranch_Diversity = "NN",
                         Low_Mg_Calcite = "NN",
                         High_Mg_Calcite = "NN",
                         Aragonite = "NN",
                         Monohydrocalcite = "NN",
                         Amorphous_Carbonate = "NN",
                         Trophic_Web_Robustness = "NN",
                         Mean_Trophic_Level = "NN",
                         
                         Turnover_Available_Biomass = "NP",
                         Selenium = "NP",
                         Zinc = "NP",
                         Omega_3 = "NP",
                         Calcium = "NP",
                         Iron = "NP",
                         Vitamin_A = "NP",
                         Available_Biomass = "NP",
                         Aesthetic = "NP",
                         Public_Interest = "NP")) # /!\ the order matter
##-------------computing PCA-------------
Contrib_site_log_transformed <- dplyr::arrange(Contrib_site_log_transformed, Biomass) 

Contrib_site_selected <- Contrib_site_log_transformed |> 
  subset(select = -c(Biomass, SiteCode, SurveyDate,
                     SiteCountry, SiteEcoregion, SurveyDepth, 
                     SiteMeanSST, SiteLatitude, SiteLongitude,
                     HDI, MarineEcosystemDependency,
                     coral_imputation, gravtot2, mpa_name,
                     mpa_enforcement, protection_status, 
                     mpa_iucn_cat))


Contrib_site_for_pca <- scale(Contrib_site_selected)

colnames(Contrib_site_for_pca)
w <- c(1/7, 1/7, 1/5, 
       1/5, 1/5, 1/5,
       1/5, 1/5, 1/5,
       1/5, 1/5, 1/5,
       1/7, 1/7, 1/7,
       1/7, 1/7, 1/2,
       1/2, 1/2, 1/6,
       1/6, 1/6, 1/6,
       1/6, 1/6, 1/2,
       1/2, 1/2)
names(w) <- colnames(Contrib_site_for_pca)
w 

pca <- FactoMineR::PCA(Contrib_site_for_pca, col.w = w,
                       scale.unit = FALSE, graph=F, ncp=9)



#-----------------------------figure 1a: NN biplot ----------------------------
#eigenvalues
eig <- as.data.frame(factoextra::get_eig(pca))

#label position
data <- data.frame(obsnames=row.names(pca$ind$coord), pca$ind$coord)
datapc <- data.frame(varnames=rownames(pca$var$coord), pca$var$coord)
mult <- min(
  (max(data[,"Dim.2"]) - min(data[,"Dim.2"])/(max(datapc[,"Dim.2"])-min(datapc[,"Dim.2"]))),
  (max(data[,"Dim.1"]) - min(data[,"Dim.1"])/(max(datapc[,"Dim.1"])-min(datapc[,"Dim.1"])))
)
datapc <- transform(datapc,
                    v1 = .7 * mult * (datapc$Dim.1),
                    v2 = .7 * mult * (datapc$Dim.2) )


plot_1a <- factoextra::fviz_pca_biplot(
  pca, 
  col.var = "forestgreen", 
  select.var = list(cos2 = 0.25,
                    name = names(grp_NN_NP)[
                      which(grp_NN_NP== "NN")]),
  geom="point", pointshape=21, 
  pointsize = 3,
  geom.var = c("arrow"),
  fill.ind = Contrib_site_log_transformed$Biomass,
  col.ind =  Contrib_site_log_transformed$Biomass,
  alpha.ind = 0.6,
  gradient.cols =
    # RColorBrewer::brewer.pal(10,name="PuOr"),
    rev(color_grad),
  # harrypotter::hp(n=9,option="RonWeasley",begin = 0.3, direction = -1),
  repel = TRUE,
  ggtheme = theme_bw(base_line_size = 0)) +
  
  #Add labels
  ggrepel::geom_label_repel(
    data= datapc[which(grp_NN_NP== "NN" &
                         factoextra::facto_summarize(
                           pca, element= "var",axes = c(1,2))$cos2 > 0.25)
                 ,], #extract cos2 of variables
    aes(x= v1*1.02, y= v2*1.02, label= gsub("_"," ", varnames)),
    size=4, fill = "white", label.size = 0.1,
    fontface = "bold", 
    color = "forestgreen", alpha = 0.8,
    force_pull = 3,
    direction = "both",
    seed = 1968)+
  # #Add boxes
  # ggrepel::geom_label_repel(
  #   data= datapc[which(grp_NN_NP== "NN" &
  #                        factoextra::facto_summarize(
  #                          pca, element= "var",axes = c(1,2))$cos2 > 0.25), ], 
  #   aes(x= v1*1.02, y= v2*1.02, label= gsub("_"," ", varnames)),
  #   alpha = 1, size=4, label.size = NA,  force_pull = 3, color = "forestgreen",
  #   fontface = "bold", direction = "both", seed = 1968)+
  
  xlim(c(-4.5,4.5)) +
  ylim(c(-3.2,4.3)) +
  labs(title = "Nature for Nature",
       x = "",
       y = paste0("PC2 (", round(eig$variance.percent[2], 1), "%)")) + 
  theme( legend.position = "none",
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 14),
         plot.margin =unit(c(0,0,0,0), 'cm'),
         plot.title = element_text(colour = "forestgreen", 
                                   face = "bold", size = 17,
                                   margin=margin(t = 20, b = -20),
                                   hjust = 0.01))
print(plot_1a)


#-----------------------------figure 1b: NP biplot ----------------------------
plot_1b <- factoextra::fviz_pca_biplot( 
  pca, 
  col.var = "dodgerblue3", 
  select.var = list(cos2 = 0,
                    name = names(grp_NN_NP)[
                      which(grp_NN_NP== "NP")]),
  geom="point", pointshape=21,
  pointsize = 3,
  geom.var = c("arrow"),
  fill.ind = Contrib_site_log_transformed$Biomass,
  col.ind = Contrib_site_log_transformed$Biomass,
  gradient.cols = 
    rev(color_grad),
  # harrypotter::hp(n=9,option="RonWeasley",begin = 0.3, direction = -1),
  #RColorBrewer::brewer.pal(10,name="PuOr"), #rev(color_grad), 
  # scico::scico(9,palette = "lajolla", begin = 0.2, direction = -1),
  alpha.ind = 0.6,
  repel = TRUE,
  ggtheme = theme_bw(base_line_size = 0)) +
  
  #Add labels
  ggrepel::geom_label_repel(
    data= datapc[which(grp_NN_NP== "NP"),],
    aes(x= v1*1.02, y= v2*1.02, label= gsub("_"," ", varnames)),
    size=4, fill = "white", label.size = 0.1,
    fontface = "bold", 
    color = "dodgerblue3", alpha = 0.8,
    force_pull = 3,
    direction = "both",
    seed = 1968)+
  # #Add boxes
  # ggrepel::geom_label_repel(
  #   data= datapc[which(grp_NN_NP== "NP"), ], 
  #   aes(x= v1*1.02, y= v2*1.02, label= gsub("_"," ", varnames)),
  #   alpha = 0.2, size=4, label.size = NA, fill = "dodgerblue3", force_pull = 3,
  #   fontface = "bold", direction = "both", seed = 1968)+
  
  
  xlim(c(-4.5,4.5)) +
  ylim(c(-3.2,4.3)) +
  labs(title = "Nature for People",
       fill = "log(total biomass)",
       x = paste0("PC1 (", round(eig$variance.percent[1], 1), "%)"),
       y = paste0("PC2 (", round(eig$variance.percent[2], 1), "%)")) +
  guides(color = "none",
         fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  theme( legend.position =c(0.15,0.13),
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 14),
         plot.margin =unit(c(0,0,0,0), 'cm'),
         legend.direction = "horizontal",
         legend.title = element_text(hjust = 0.5, size = 14),
         legend.title.align = 0,
         legend.key.size = unit(0.8, 'cm'),
         legend.text = element_text(size = 10),
         plot.title = element_text(colour = "dodgerblue3",
                                   face = "bold", size = 17,
                                   margin=margin(t = 20, b = -20),
                                   hjust = 0.01))


plot_1b

##------------------------Insert figure 1a: variables eigenvalues ------------------------------
eig <- factoextra::get_eig(pca)
variance_explained <- data.frame(
  Axe = c(1:nrow(eig)), 
  cumulative_variance_explained = eig[,3],
  contribution_coefficient = eig[,2]/100)

elbow_values <- elbow(variance_explained[,c("Axe", "cumulative_variance_explained")])

elbow_point <- c( elbow_values$"Axe"[tail(which(elbow_values$"SelectedorNot" =="Selected"), n = 1)],
                  elbow_values$"cumulative_variance_explained"[tail(
                    which(elbow_values$"SelectedorNot" =="Selected"), n = 1)])


ndim=15
elbow_plot <- ggplot(variance_explained, 
                     aes(x = Axe, 
                         y = cumulative_variance_explained,
                         color = "black")) + 
  #barchart of eigenvalues
  geom_bar(data = as.data.frame(eig)[c(1:ndim),], 
           aes(x= c(1:ndim), y = variance.percent), 
           width = 0.8,
           stat= "identity",
           col = "grey50",
           fill = "grey50")+
  
  #cumulative curve
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5, alpha = 0.8) +
  stat_summary(fun = "mean", linewidth = 0.1) +
  scale_y_continuous(labels = c("", "25", "50", "75", "100%"))+
  
  # elbow point
  geom_segment(aes(x = elbow_point[1], xend = elbow_point[1],
                   y = 0 , yend = elbow_point[2]),
               color = "grey30", linetype = "dotted", linewidth = 1) +
  
  geom_segment(aes(y = elbow_point[2], yend = elbow_point[2],
                   x = 0, xend = elbow_point[1]),
               color = "grey30", linetype = "dotted", linewidth = 1) +
  
  geom_point(aes(y = elbow_point[2], x = elbow_point[1]),
             color = "grey30", size = 4, shape = 19)+

  labs(x = "", y = "") +
  theme_bw(base_line_size = 0) +
  theme(plot.margin = unit(c(1,0.3,-.5,-.5), 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color=NA),
        panel.background = element_rect(fill = alpha("white", 0.8), color=NA), 
        strip.text.x = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold", size = 16,
                                  margin=margin(t=-20, b = 5),
                                  hjust = -.1),
        legend.position = "none") +
  coord_cartesian(expand = FALSE, xlim = c(0, ndim), ylim = c(0, 100))+
  harrypotter::scale_colour_hp_d(option = "DracoMalfoy", alpha = 0.8)



elbow_plot # to see the complete figure of the elbow, refer to '1c_PCA_analyses_on_contributions.R'



#------------------------figure 1c: contributions ------------------------------
var <- factoextra::get_pca_var(pca)
contributions <- var$contrib
for( i in 1:ncol(contributions)){ 
  contributions[,i] <- contributions[,i] * 
    variance_explained$contribution_coefficient[i]
}

colnames(contributions) <- c(1:ncol(contributions))

#Order contributions according to the first dimension
contributions <- cbind(contributions, category = grp_NN_NP[rownames(contributions)])
contributions <- contributions[order( contributions[,"category"], -contributions[,1]),]
contributions <- contributions[, colnames(contributions) !="category"]

#set color background
color <- var$coord
color_bg <- color[rownames(contributions),]
color_bg <- ifelse(color_bg>0, "white", "grey90")


#divide into two plots
contributions_NN <- contributions
contributions_NN[names(grp_NN_NP)[ grp_NN_NP!="NN" ],] <- 0
row.names(contributions_NN) <- gsub("_", " ",  row.names(contributions_NN))


contributions_NP <- contributions
contributions_NP[names(grp_NN_NP)[ grp_NN_NP!="NP" ],] <- NA
row.names(contributions_NP) <- gsub("_", " ",  row.names(contributions_NP))



# plot the figure 1c
png(filename = here::here("outputs", "figures","FIG_1C_contribution_in_total_variance_NN_NP.png"), 
    width=11, height = 21.5, units = "cm", res = 1000)
par(mar = c(2,0,0,0), xaxs = "i", yaxs = "i")
corrplot::corrplot(contributions_NN, is.corr=FALSE, 
                   bg = color_bg,
                   col = c("forestgreen"), tl.col = "black", 
                   cl.pos = "n", 
                   tl.cex = 1,
                   tl.srt = 0, 
                   na.label.col = "transparent",
                   mar=c(2,0,2,0))

corrplot::corrplot(contributions_NP, is.corr=FALSE,
                   col = c("dodgerblue3"), tl.col = "transparent",
                   cl.pos = "n", tl.cex = 1.2, bg = "transparent",
                   tl.srt = 60, na.label.col = "transparent", 
                   add = TRUE)
lines(x=c(-9,9.5), y=rep(10.5, 2), lty = 1, lwd = 2)
text(5,30.7, label= "Dimensions (PC)", cex=1.2 , font=2 ) 
for(i in c(1:ncol(contributions))){
  text(x=-0.2+i, y= -0.5, label= paste(round(colSums(contributions)[i],0), "%"), 
       cex=0.9, srt=60, font=1)}
text(0,-0.6, label= "Variance explained:", cex=1 , font=3, pos=2 )

dev.off()



#----------------------------- Plot Figure 1  ----------------------------
#insert the elbow figure into fig 1a
plot_merged <- plot_1a +
  annotation_custom( ggplotGrob(elbow_plot),
                     xmin = -4.6,
                     xmax = -1,
                     ymin = -3.5,
                     ymax = 0)
plot_merged

#load figure 1c
corrplot_fig1c <- cowplot::ggdraw() + 
  cowplot::draw_image(here::here("outputs","figures",
                                 "FIG_1C_contribution_in_total_variance_NN_NP.png"))+
  theme(plot.margin =unit(c(0,0,0,-0.5), 'cm'))


#Assemble the figures
png(filename = here::here("outputs", "figures","FIGURE_1.png"), 
    width=40, height = 27, units = "cm", res = 1000)
gridExtra::grid.arrange(
  plot_merged, plot_1b, corrplot_fig1c, 
  ncol=2,
  layout_matrix = rbind(c(1,1,1,3,3),
                        c(1,1,1,3,3),
                        c(1,1,1,3,3),
                        c(2,2,2,3,3),
                        c(2,2,2,3,3),
                        c(2,2,2,3,3)))
# Add labels A, B and C
grid::grid.text("A", x=unit(0.01, "npc"), y=unit(0.98, "npc"), 
                just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
grid::grid.text("B", x=unit(0.01, "npc"), y=unit(0.48, "npc"), 
                just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
grid::grid.text("C", x=unit(0.65, "npc"), y=unit(0.98, "npc"), 
                just="left", gp=grid::gpar(fontsize=17, fontface="bold"))

dev.off()

