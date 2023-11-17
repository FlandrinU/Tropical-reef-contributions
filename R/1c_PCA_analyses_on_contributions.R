################################################################################
##
## This script study the dimensionality of contributions in the RLS reefs, using
##  Principal Component Analyses (PCA).
##
## 1c_PCA_analyses_on_contributions.R
##
## 03/01/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
# pkgs <- c("here", "tidyverse","FactoMineR", "tibble", "questionr", "corrplot",
#           "factoextra", "ggpubr", "scico", "RColorBrewer", "plotly", "fishualize", 
#           "ggplot2", "patchwork", "colormap", "grDevices", "ggnewscale", "sf",
#           "rdacca.hp", "harrypotter", "grid", "gridExtra", "ade4")
# nip <- pkgs[!(pkgs %in% installed.packages())]
# nip <- lapply(nip, install.packages, dependencies = TRUE)
library(ggplot2)
library(patchwork)

rm(list=ls())

##-------------loading data-------------
load(here::here("outputs","all_Contrib_site.Rdata"))
load(here::here("outputs","all_Contrib_site_log_transformed.Rdata"))
### To make test on other filtering condition, uncomment the wanted lines
# load(here::here("outputs","Contrib_site_log_coral_reef.Rdata"))
# load(here::here("outputs","Contrib_site_log_SST20.Rdata"))
# load(here::here("outputs","Contrib_site_log_coral_5_imputed.Rdata"))
# load(here::here("outputs","Contrib_site_log_wo_australia.Rdata"))
# load(here::here("outputs","Contrib_site_log_random.Rdata"))
# load(here::here("outputs","Contrib_site_log_only_australia.Rdata"))

# Contrib_site_log_transformed <- Contrib_site_condition
###


load(here::here("biodiversity", "outputs", "occurrence_matrix_family_survey_relative_biomass.Rdata"))
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

source(here::here("R", "elbow.R"))

##------------Classify contributions into 2 categories: 
##              Nature for Nature (NN) and Nature for People (NP) -------------
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
                         
                         Productivity = "NP",
                         Selenium = "NP",
                         Zinc = "NP",
                         Omega_3 = "NP",
                         Calcium = "NP",
                         Iron = "NP",
                         Vitamin_A = "NP",
                         Fishery_Biomass = "NP",
                         Aesthetic = "NP",
                         Public_Interest = "NP")) # /!\ the order matter



##-------------computing PCA-------------
Contrib_site_selected <- subset(Contrib_site_log_transformed, 
                                select = -c(SiteCode, SurveyDate,
                                            SiteCountry, SiteEcoregion, SurveyDepth, 
                                            SiteMeanSST, SiteLatitude, SiteLongitude,
                                            HDI, MarineEcosystemDependency,
                                            coral_imputation, gravtot2, mpa_name,
                                            mpa_enforcement, protection_status, 
                                            mpa_iucn_cat,
                                            Biomass))


Contrib_site_for_pca <- scale(Contrib_site_selected)


### PCA with weigths according categories and items (see table 1, Flandrin et al. 2023)
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


####----------plot PCA with all Contributions at the community scale -------------

##----------------------Dimensionality of the contributions -----------------
eig <- factoextra::get_eig(pca)
variance_explained <- data.frame(
  Axe = c(1:nrow(eig)), 
  cumulative_variance_explained = eig[,3],
  contribution_coefficient = eig[,2]/100)

elbow_values <- elbow(variance_explained[,c("Axe", "cumulative_variance_explained")])

elbow_point <- c( elbow_values$"Axe"[tail(which(elbow_values$"SelectedorNot" =="Selected"), n = 1)],
                  elbow_values$"cumulative_variance_explained"[tail(
                    which(elbow_values$"SelectedorNot" =="Selected"), n = 1)])


ndim=20
elbow_plot <- ggplot(variance_explained, aes(x = Axe, y = cumulative_variance_explained,
                                             color = "black")) + 
  #barchart of eigenvalues
  geom_bar(data = as.data.frame(eig)[c(1:ndim),], 
           aes(x= c(1:ndim), y = variance.percent), 
           stat= "identity",
           col = "grey30",
           fill = "grey30")+
  
  #cumulative curve
  stat_summary(fun = "mean", geom = "line", linewidth = 1, alpha = 0.4) +
  stat_summary(fun = "mean", linewidth = 0.5) +
  
  labs(x = "Number of dimensions") + 
  labs(y = "Variance explained (in %)") +
  theme_bw(base_line_size = 0) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        panel.background = element_blank(),
        legend.position = "none") +
  coord_cartesian(expand = FALSE, xlim = c(0, ndim), ylim = c(0, 100))+
  harrypotter::scale_colour_hp_d(option = "LunaLovegood") +
  
  # elbow point
  geom_segment(aes(x = elbow_point[1], xend = elbow_point[1],
                   y = 0 , yend = elbow_point[2]),
               color = "black", linetype = "dotted", linewidth = 1) +
  
  geom_segment(aes(y = elbow_point[2], yend = elbow_point[2],
                   x = 0, xend = elbow_point[1]),
               color = "black", linetype = "dotted", linewidth = 1) +
  
  geom_point(aes(y = elbow_point[2], x = elbow_point[1]),
             color = "black", size = 4, shape = 19)+
  geom_label(aes(label = paste0("Elbow rule: ", elbow_point[1], " dimensions \n", 
                                "Variance explained = ", round(elbow_point[2],1) , " %"),
                 y = elbow_point[2]-5, x = elbow_point[1]+1), size = 4,
             color = "black", hjust = 0)+
  geom_text(data = variance_explained[c(1:8),],
            aes(x=Axe + 0.5, y = cumulative_variance_explained,
                label = paste(round(cumulative_variance_explained, 1), "%")))


elbow_plot
ggsave(filename = here::here("outputs", "figures",
                             "Variance explained by ACP_elbow rule.png"), elbow_plot, 
       width = 15, height =10 )



##-----------------Importance of Contributions in dimensions ---------------
var <- factoextra::get_pca_var(pca)
contributions <- var$contrib
for( i in 1:ncol(contributions)){ 
  contributions[,i] <- contributions[,i] * variance_explained$contribution_coefficient[i]}

png(filename = here::here("outputs", "figures","Importance_of_Contributions_in_total_variance.png"), 
    width= 12, height = 15, units = "cm", res = 1000)
print( corrplot::corrplot(contributions, is.corr=FALSE, col=rev(RColorBrewer::brewer.pal(10,name="RdYlBu"))))
dev.off()


#--------------Study the covariations of the contributions in the PCA ---------------
#### PCA in the 2 first dimensions, with representation quality ($cos^{2}$) of each variables
png(filename = here::here("outputs", "figures","PCA_all_Contrib.png"), 
    width= 20, height = 17, units = "cm", res = 1000)
print( factoextra::fviz_pca_var(pca, col.var = "cos2",
                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                repel = TRUE 
))
dev.off()


## Classify variables in Nature for Nature (NN) and Nature for People (NP)
png(filename = here::here("outputs", "figures","PCA_NP_NN_axes1-2_blind.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
  factoextra::fviz_pca_biplot(pca, col.var = grp_NN_NP, 
                              title="",
                              palette = c("forestgreen", "dodgerblue3"),
                              labelsize = 4, 
                              repel = TRUE,
                              geom="point", pointshape=21,
                              stroke=0, pointsize=3,
                              alpha.ind = 0.7,
                              fill.ind = Contrib_site_log_transformed$Biomass,
                              gradient.cols = RColorBrewer::brewer.pal(9,name="YlOrRd"))+
  labs(col = "Reef Contributions", fill = "Total Biomass")+
  scale_color_discrete(type=c("forestgreen", "dodgerblue3"),labels = c("NN", "NP"))+
  theme(legend.position = "bottom")

dev.off()


#### PCA in dimensions 3 and 4
png(filename = here::here("outputs", "figures","PCA_NP_NN_axes3-4.png"), 
    width= 30, height = 25, units = "cm", res = 1000)
print( factoextra::fviz_pca_var(pca, col.var = grp_NN_NP, axe= c(3,4), repel = TRUE,
                                palette = c("forestgreen", "dodgerblue3"),
                                legend.title = "Nature Based Contributions",
                                select.var = list(cos2 = 0)))
dev.off() 



#------------------------ Study the sites distribution--------------------------
### Biomass: see "PCA_NP_NN_axes1-2_blind.png"
### -> which correlation between PC1_coordinates and Biomass ?
pc1_coord <- pca$ind$coord[,1]
biom <- Contrib_site_log_transformed$Biomass
plot(pc1_coord ~ biom)
cor.test(pc1_coord, biom)
# cor 
# 0.8559888 
l <- lm(pc1_coord ~ biom)
summary(l)


### Survey date
png(filename = here::here("outputs", "figures","PCA_date_pattern.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
print( factoextra::fviz_pca_biplot(pca,
                                   select.var = list(cos2 = 0.4),
                                   geom="point", pointshape=21,
                                   fill.ind = as.numeric(stringr::str_sub(Contrib_site_log_transformed$SurveyDate , start = -4)),
                                   col.ind = as.numeric(stringr::str_sub(Contrib_site_log_transformed$SurveyDate , start = -4)),
                                   gradient.cols = rev(RColorBrewer::brewer.pal(10,name="RdYlBu")),
                                   col.var = "black", repel = TRUE,
                                   legend.title = "Year"))
dev.off()



date_axis <- data.frame(y1=pca$ind$coord[,1],
                   y2=pca$ind$coord[,2],
                   x=as.numeric(stringr::str_sub(Contrib_site_log_transformed$SurveyDate , start = -4)))

plotPC1 <- ggplot(date_axis, aes(x =x, y=y1 , alpha = 0.5)) +
  geom_point()+
  geom_smooth(method='lm', formula= y~x) +
  labs(x = "year" , y = "PC" , title = "") +
  theme_minimal()+
  theme(legend.position = 'none')+
  ggpubr::stat_regline_equation()+
  ggpubr::stat_cor(label.y = 4)
  
plotPC2 <-  ggplot(date_axis, aes(x =x, y=y2 , alpha = 0.5)) +
  geom_point()+
  geom_smooth(method='lm', formula= y~x) +
  labs(x = "year" , y = "PC" , title = "") +
  theme_minimal()+
  theme(legend.position = 'none')+
  ggpubr::stat_regline_equation()+
  ggpubr::stat_cor(label.y = 4) 

date_plot <- plotPC1 + plotPC2 +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))
ggsave(filename = here::here("outputs", "figures", "PCA_dimensions_according_to_date.png"), 
       date_plot, width = 15, height = 12)

lm <- lm(date_axis$y2 ~ date_axis$x)
summary(lm)


### mpa categories
png(filename = here::here("outputs", "figures","PCA_mpa_protection_status.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
print( factoextra::fviz_pca_biplot(pca,
                                   select.var = list(cos2 = 0.4),
                                   geom="point", pointshape=21,
                                   addEllipses = F,
                                   fill.ind = Contrib_site_log_transformed$protection_status,
                                   col.ind = Contrib_site_log_transformed$protection_status,
                                   gradient.cols = rev(RColorBrewer::brewer.pal(10,name="RdYlBu")),
                                   col.var = "black", repel = TRUE,
                                   legend.title = "Protection status of MPAs"))
dev.off()


### By countries
ind.p <- factoextra::fviz_pca_ind(pca, geom = "point",
                                  pointshape=21,
                                  col.ind = "grey20",
                                  fill.ind = Contrib_site_log_transformed$SiteCountry)

png(filename = here::here("outputs", "figures","PCA_countries.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
print( ggpubr::ggpar(ind.p,
                     title = "Principal Component Analysis of Contributions",
                     subtitle = "RLS data set",
                     caption = "Source: factoextra",
                     xlab = "PC1", ylab = "PC2",
                     legend.title = "Country", legend.position = "top",
                     font.legend = c(8, "plain", "black"),
                     ggtheme = theme_minimal(),
                     palette = scico::scico(38, palette = "roma")))
dev.off() 



### Temperature pattern
png(filename = here::here("outputs", "figures","PCA_temperature_pattern.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
print( factoextra::fviz_pca_biplot(pca,
                                   select.var = list(cos2 = 0.4),
                                   geom="point", pointshape=21,
                                   fill.ind = Contrib_site_log_transformed$SiteMeanSST,
                                   col.ind = Contrib_site_log_transformed$SiteMeanSST,
                                   gradient.cols = rev(RColorBrewer::brewer.pal(10,name="RdYlBu")),
                                   col.var = "black", repel = TRUE,
                                   legend.title = "Mean temperature"))
dev.off()

df <- data.frame(x = Contrib_site_log_transformed$SiteMeanSST, y = pca$ind$coord[,"Dim.1"])
temp <- ggplot(df, aes(x , y , alpha = 0.5)) +
  geom_point()+
  geom_smooth() +
  labs(x = "Mean SST", y = "Dimension 1 in global PCA", title = "Importance of SST in PCA's axe 1") +
  theme_minimal()+
  theme(legend.position = 'none')
ggsave(filename = here::here("outputs", "figures", "PCA_temp_according_Dim1.png"), temp, width = 15, height = 12)



### Lattitude pattern
png(filename = here::here("outputs", "figures","PCA_lattitude_pattern.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
print( factoextra::fviz_pca_biplot(pca,
                                   geom="point", pointshape=21,
                                   fill.ind = abs(Contrib_site_log_transformed$SiteLatitude),
                                   col.ind = abs(Contrib_site_log_transformed$SiteLatitude),
                                   gradient.cols = RColorBrewer::brewer.pal(10,name="RdYlBu"),
                                   col.var = "black", repel = TRUE,
                                   legend.title = "abs(latitude)") )
dev.off()


##------- PCA co-inertia between NN and NP -------
pca_NP <- ade4::dudi.pca(Contrib_NP, center = T, scale = T, scannf=F, nf=10)
pca_NN <- ade4::dudi.pca(Contrib_NN, center = T, scale = T, scannf=F, nf=10)

coinertia <- ade4::coinertia(pca_NN, pca_NP, scannf = F, nf=3)
coinertia
summary(coinertia) #RV: 0.4567662 
plot(coinertia)
#strengh of the relationship
rv1 <- ade4::RV.rtest(pca_NP$tab, pca_NN$tab, 99)
plot(rv1)


#Mantel test between the two contributions tables
mantel_test <- vegan::mantel(dist(scale(Contrib_NP)), dist(scale(Contrib_NN)), permutations = 100)
mantel_test
#Mantel statistic r: 0.5443 
