################################################################################
##
## Makes figure 1
##
## make_fig_1.R
##
## 19/04/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse","FactoMineR", "corrplot",
          "factoextra", "ggplot2", "harrypotter", "grid", "gridExtra")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load(here::here("outputs","all_NCP_site_log_transformed.Rdata"))

source(here::here("R", "elbow.R"))

##-------------NCP in categories-------------
## Classify variables in Nature for Nature (NN) and Nature for Society (NS)
grp_NN_NS <- as.factor(c(N_recycling = "NN",
                         P_recycling = "NN",
                         Taxonomic_Richness = "NN",
                         Functional_Entropy = "NN", 
                         Phylogenetic_Entropy = "NN",
                         Functional_Distinctiveness = "NN",
                         Evolutionary_distinctiveness = "NN",
                         Low_TL_Biomass = "NN",
                         Medium_TL_Biomass = "NN",
                         High_TL_Biomass = "NN",
                         # IUCN_Species = "NN",
                         Endemism = "NN",
                         Elasmobranch_Diversity = "NN",
                         Low_Mg_Calcite = "NN",
                         High_Mg_Calcite = "NN",
                         Aragonite = "NN",
                         Monohydrocalcite = "NN",
                         Amorphous_Carbonate = "NN",
                         Trophic_web_robustness = "NN",
                         mean_Trophic_Level = "NN",
                         
                         Productivity = "NS",
                         Selenium = "NS",
                         Zinc = "NS",
                         Omega_3 = "NS",
                         Calcium = "NS",
                         Iron = "NS",
                         Vitamin_A = "NS",
                         Fishery_Biomass = "NS",
                         Aesthetic = "NS",
                         Public_Interest = "NS")) # /!\ the order matter
##-------------computing PCA-------------
NCP_site_selected <- subset(NCP_site_log_transformed, 
                            select = -c(Biomass, SiteCode, SurveyDate,
                                        SiteCountry, SiteEcoregion, SurveyDepth, 
                                        SiteMeanSST, SiteLatitude, SiteLongitude,
                                        HDI, MarineEcosystemDependency,
                                        coral_imputation, gravtot2, mpa_name,
                                        mpa_enforcement, protection_status, 
                                        mpa_iucn_cat))


NCP_site_for_pca <- scale(NCP_site_selected)
pca <- FactoMineR::PCA(NCP_site_for_pca, scale.unit = FALSE, graph=F, ncp=9) 


#-----------------------------figure 1a: NN biplot ----------------------------
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

library(ggplot2)

plot_1a <- factoextra::fviz_pca_biplot(
  pca, 
  col.var = "forestgreen", 
  select.var = list(cos2 = 0.25,
                    name = names(grp_NN_NS)[
                      which(grp_NN_NS== "NN")]),
  geom="point", pointshape=21, 
  pointsize = 3,
  geom.var = c("arrow"),
  fill.ind = NCP_site_log_transformed$Biomass,
  col.ind =  NCP_site_log_transformed$Biomass,
  alpha.ind = 0.7,
  gradient.cols = rev(
    RColorBrewer::brewer.pal(10,name="RdYlBu")),
  repel = TRUE,
  ggtheme = theme_bw(base_line_size = 0)) +
  
  #Add labels
  ggrepel::geom_label_repel(
    data= datapc[which(grp_NN_NS== "NN" &
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
  #   data= datapc[which(grp_NN_NS== "NN" &
  #                        factoextra::facto_summarize(
  #                          pca, element= "var",axes = c(1,2))$cos2 > 0.25), ], 
  #   aes(x= v1*1.02, y= v2*1.02, label= gsub("_"," ", varnames)),
  #   alpha = 1, size=4, label.size = NA,  force_pull = 3, color = "forestgreen",
  #   fontface = "bold", direction = "both", seed = 1968)+

  xlim(c(-7,8)) +
  ylim(c(-6,6)) +
  labs(fill="log(Biomass)",
       title = "Nature to Nature") +
  xlab("")+
  theme( legend.position = "none",
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 14),
         plot.margin =unit(c(0,0,0,0), 'cm'),
         plot.title = element_text(colour = "forestgreen", 
                                   face = "bold", size = 17,
                                   margin=margin(t = 20, b = -20),
                                   hjust = 0.01))
print(plot_1a)


#-----------------------------figure 1b: NS biplot ----------------------------

plot_1b <- factoextra::fviz_pca_biplot( 
  pca, 
  col.var = "dodgerblue3", 
  select.var = list(cos2 = 0,
                    name = names(grp_NN_NS)[
                      which(grp_NN_NS== "NS")]),
  geom="point", pointshape=21,
  pointsize = 3,
  geom.var = c("arrow"),
  fill.ind = NCP_site_log_transformed$Biomass,
  col.ind = NCP_site_log_transformed$Biomass,
  alpha.ind = 0.7,
  gradient.cols = rev(RColorBrewer::brewer.pal(10,name="RdYlBu")),
  repel = TRUE,
  ggtheme = theme_bw(base_line_size = 0)) +
  
  #Add labels
  ggrepel::geom_label_repel(
    data= datapc[which(grp_NN_NS== "NS"),],
    aes(x= v1*1.02, y= v2*1.02, label= gsub("_"," ", varnames)),
    size=4, fill = "white", label.size = 0.1,
    fontface = "bold", 
    color = "dodgerblue3", alpha = 0.8,
    force_pull = 3,
    direction = "both",
    seed = 1968)+
  # #Add boxes
  # ggrepel::geom_label_repel(
  #   data= datapc[which(grp_NN_NS== "NS"), ], 
  #   aes(x= v1*1.02, y= v2*1.02, label= gsub("_"," ", varnames)),
  #   alpha = 0.2, size=4, label.size = NA, fill = "dodgerblue3", force_pull = 3,
  #   fontface = "bold", direction = "both", seed = 1968)+

  
  xlim(c(-7,8)) +
  ylim(c(-6,6)) +
  labs(title = "Nature to People",
       fill = "log(Biomass)") +
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
print(plot_1b)

##------------------------figure 1c: variables eigenvalues ------------------------------
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
           stat= "identity",
           col = "grey30",
           fill = "grey30")+
  
  #cumulative curve
  stat_summary(fun = "mean", geom = "line", size = 0.5, alpha = 1) +
  stat_summary(fun = "mean", size = 0.1) +
  scale_y_continuous(labels = c("", "25", "50", "75", "100%"))+
  
  # elbow point
  geom_segment(aes(x = elbow_point[1], xend = elbow_point[1],
                   y = 0 , yend = elbow_point[2]),
               color = "black", linetype = "dotted", linewidth = 1) +
  
  geom_segment(aes(y = elbow_point[2], yend = elbow_point[2],
                   x = 0, xend = elbow_point[1]),
               color = "black", linetype = "dotted", linewidth = 1) +
  
  geom_point(aes(y = elbow_point[2], x = elbow_point[1]),
             color = "black", size = 4, shape = 19)+
  # geom_label(aes(label = paste0("Elbow rule: ", elbow_point[1], " dimensions \n", 
  #                               "Variance explained = ", round(elbow_point[2],1) , " %"),
  #                y = elbow_point[2]-5, x = elbow_point[1]+1), size = 4,
  #            color = "black", hjust = 0)

  # labs(x = "Number of dimensions") + 
  # labs(y = "Variance explained (in %)") +
  labs(x = "", y = "") +
  theme_bw(base_line_size = 0) +
  theme(plot.margin = unit(c(1,0.3,-.5,-.5), 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color=NA),
        panel.background = element_rect(fill = alpha("white", 0.7), color=NA), 
        strip.text.x = element_blank(),
        axis.text = element_text(size = 12, face = "bold"),
        # axis.title.x = element_text(size = 15, face = "bold",
        #                             hjust = 1.12, vjust = 9),
        plot.title = element_text(face = "bold", size = 16,
                                  margin=margin(t=-20, b = 5),
                                  hjust = -.1),
        legend.position = "none") +
  coord_cartesian(expand = FALSE, xlim = c(0, ndim), ylim = c(0, 100))+
  harrypotter::scale_colour_hp_d(option = "DracoMalfoy")
    


elbow_plot

# ggsave(filename = here::here("outputs", "figures",
#                              "Variance explained by ACP_elbow rule.png"), elbow_plot, 
#        width = 15, height =10 )

##------------------------figure 1d: contributions ------------------------------
# #divide into two plot
# contributions_NN <- contributions[ names(grp_NN_NS)[ grp_NN_NS=="NN" ] ,] 
# contributions_NN <- contributions_NN[order(-contributions_NN[,"Dim.1"]), ]
# row.names(contributions_NN) <- gsub("_", " ",  row.names(contributions_NN))
# par(mar = c(0,0,0,0), new = T)
# corrplot::corrplot(contributions_NN, is.corr=FALSE, 
#                      col = c("forestgreen"), tl.col = "black", 
#                      cl.pos = "n", tl.cex = 0.7)
# NN_corr <- grDevices::recordPlot()
# 
# contributions_NS <- contributions[ names(grp_NN_NS)[ grp_NN_NS=="NS" ] ,] 
# contributions_NS <- contributions_NS[order(-contributions_NS[,"Dim.1"]), ]
# row.names(contributions_NS) <- gsub("_", " ",  row.names(contributions_NS))
# row.names(contributions_NS)[1] <- paste(
#   c(rep(" ", 
#       max(nchar(row.names(contributions_NN))) - nchar(row.names(contributions_NS)[1])+3),
#   row.names(contributions_NS)[1]),
#   collapse="") #change legend length to have same width in both plot
# par(mar = c(0,0,0,0), new = T)
# corrplot::corrplot(contributions_NS, is.corr=FALSE,
#                      col = c("dodgerblue3"), tl.col = "black", 
#                      cl.pos = "n", tl.pos = 'l', tl.cex = 0.7)
# NS_corr <- grDevices::recordPlot()
# 
# #merge the two plot
# png(filename = here::here("outputs", "figures","contribution_in_total_variance_NN_NS.png"), 
#     width=10, height = 24, units = "cm", res = 1000)
# print(cowplot::plot_grid(NN_corr, NS_corr, ncol=1, 
#                    rel_heights = c(nrow(contributions_NN)+5, nrow(contributions_NS)), #change 0.7 to ajust relative sizes
#                    align= "hv"))
# dev.off()


var <- factoextra::get_pca_var(pca)
contributions <- var$contrib
for( i in 1:ncol(contributions)){ 
  contributions[,i] <- contributions[,i] * variance_explained$contribution_coefficient[i]}

colnames(contributions) <- paste0(colnames(contributions),": ", round(colSums(contributions),0), "%")

#divide into two plot
contributions_NN <- contributions[ names(grp_NN_NS)[ grp_NN_NS=="NN" ] ,] 
contributions_NN <- contributions_NN[order(-contributions_NN[,1]), ]
row.names(contributions_NN) <- gsub("_", " ",  row.names(contributions_NN))

contributions_NS <- contributions[ names(grp_NN_NS)[ grp_NN_NS=="NS" ] ,] 
contributions_NS <- contributions_NS[order(-contributions_NS[,1]), ]
row.names(contributions_NS) <- gsub("_", " ",  row.names(contributions_NS))
row.names(contributions_NS)[1] <- paste(
  c(rep(" ", 
        max(nchar(row.names(contributions_NN))) - nchar(row.names(contributions_NS)[1])+4),
    row.names(contributions_NS)[1]),
  collapse="") #change legend length to have same width in both plot

# plot pannel
png(filename = here::here("outputs", "figures","contribution_in_total_variance_NN_NS.png"), 
    width=11, height = 20, units = "cm", res = 1000)
layout(matrix(c(rep(1, nrow(contributions_NN)+2), rep(2, nrow(contributions_NS))), #change proportion of both plot with +2
              ncol = 1))
par(mar = c(0,0,0,0))
corrplot::corrplot(contributions_NN, is.corr=FALSE, 
                   col = c("forestgreen"), tl.col = "black", 
                   cl.pos = "n", tl.cex = 1.2,
                   tl.srt = 60)
# text(-3,21,"D)",cex=1.5)
corrplot::corrplot(contributions_NS, is.corr=FALSE,
                   col = c("dodgerblue3"), tl.col = "black", 
                   cl.pos = "n", tl.pos = 'l', tl.cex = 1.2)  
dev.off()


#----------------------------- Plot panel Figure 1  ----------------------------
#merge 1a and 1c
plot_merged <- plot_1a +
  annotation_custom( ggplotGrob(elbow_plot),
                     xmin = -7.5,
                     xmax = -1.2,
                     ymin = -6.5,
                     ymax = -0.3)
plot_merged


corrplot_fig1c <- cowplot::ggdraw() + 
  cowplot::draw_image(here::here("outputs","figures",
                                 "contribution_in_total_variance_NN_NS.png"))+
  theme(plot.margin =unit(c(0,0,0,-0.5), 'cm'))
         

# fig1 <- gridExtra::grid.arrange(
#   plot_1a, plot_1b, elbow_plot, corrplot_fig1c, 
#   ncol=2,
#   layout_matrix = rbind(c(1,1,1,1,3,3,3),
#                         c(1,1,1,1,3,3,3),
#                         c(1,1,1,1,4,4,4),
#                         c(2,2,2,2,4,4,4),
#                         c(2,2,2,2,4,4,4),
#                         c(2,2,2,2,4,4,4)))                                  



png(filename = here::here("outputs", "figures","Panel_fig_1.png"), 
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
    # Add labels 
    grid::grid.text("A", x=unit(0.01, "npc"), y=unit(0.98, "npc"), 
                    just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
    grid::grid.text("B", x=unit(0.01, "npc"), y=unit(0.48, "npc"), 
                    just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
    grid::grid.text("C", x=unit(0.65, "npc"), y=unit(0.98, "npc"), 
                    just="left", gp=grid::gpar(fontsize=17, fontface="bold"))
    
dev.off()


