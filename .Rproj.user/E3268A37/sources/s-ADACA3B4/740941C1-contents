################################################################################
##
## From a dataframe of all NCPs, analyse trade-offs and covariations of NCPs
##  in RLS sites
##
## NCP_PCA.R
##
## 31/10/2022
##
## Ulysse Flandrin
##
################################################################################
#----------------- cleaning memory -----------------
rm(list=ls())

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse","FactoMineR", "tibble", "questionr", "corrplot",
          "factoextra", "ggpubr", "scico")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

##-------------loading data-------------
load(here::here("outputs","all_NCP_site.Rdata"))
load(here::here("outputs","all_NCP_surveys.Rdata"))

##-------------Compute PCA at site scale-------------
rownames(NCP_site) <- NCP_site$SiteCode

pca <- FactoMineR::PCA(NCP_site[,-c(1:7)], scale.unit = T)
print(pca)
summary(pca)

png( here::here("figures", "Contribution of dimensions in variances.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 30))
dev.off()

png( here::here( "figures", "PCA_NCP_axes1-2.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte 
) #give the quality of representatin of each variables in the 2 first dimensions
dev.off()

png( here::here( "figures", "PCA_NCP_axes1-3.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_var(pca, col.var = "cos2",
             axe=c(1,3),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
dev.off()

png( here::here("figures", "PCA_NCP_axes1-4.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_var(pca, col.var = "cos2",
             axe=c(1,4),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
dev.off()

png( here::here("figures", "PCA_NCP_axes2-5.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_var(pca, col.var = "cos2",
             axe=c(2,5),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
dev.off()

png( here::here("figures", "PCA_NCP_axes3-4.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_var(pca, col.var = "cos2",
             axe=c(3,4),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
dev.off()

var <- get_pca_var(pca)
png( here::here("figures", "Contribution of NCP in dimensions.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
corrplot(var$contrib, is.corr=FALSE)    
dev.off()

#Countries distribution
png( here::here("figures", "Sites in ACP by contry.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
ind.p <- fviz_pca_ind(pca, geom = "point",
                      pointshape=21, 
                      fill.ind = NCP_site$SiteCountry)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis of NCPs",
              subtitle = "RLS data set",
              caption = "Source: factoextra",
              xlab = "PC1", ylab = "PC2",
              legend.title = "Country", legend.position = "top",
              ggtheme = theme_gray(), palette = scico(38, palette = "roma")
)
dev.off()

png( here::here("figures", "biplot sites_NCP.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_biplot (pca,
                 geom="point", pointshape=21,
                 col.ind = NCP_site$SiteCountry,
                 fill.ind = NCP_site$SiteCountry,
                 palette = scico(38, palette = "roma"),
                 addEllipses = F, label = "var",
                 col.var = "black", repel = TRUE,
                 legend.title = "Country")
dev.off()

png( here::here("figures", "biplot sites_NCP with Elipse by contry.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_biplot (pca,
                 geom="point", pointshape=21,
                 col.ind = NCP_site$SiteCountry,
                 fill.ind = NCP_site$SiteCountry,
                 palette = rainbow(38),
                 addEllipses = T, label = "var",
                 col.var = "black", repel = TRUE,
                 legend.title = "Country")
dev.off()

#Temperature pattern
png( here::here("figures", "biplot sites_NCP by temperature.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_biplot(pca, 
                geom="point", pointshape=21,
                fill.ind = NCP_site$SiteMeanSST,
                col.ind = NCP_site$SiteMeanSST,
                gradient.cols = c("blue", "yellow", "red"),
                col.var = "black", repel = TRUE,
                legend.title = "Mean temperature")
dev.off()

#Depth pattern
png( here::here("figures", "biplot sites_NCP by depth.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_biplot(pca, 
                geom="point", pointshape=21,
                fill.ind = NCP_site$SurveyDepth,
                col.ind = NCP_site$SurveyDepth,
                gradient.cols = c("blue", "yellow", "red"),
                col.var = "black", repel = TRUE,
                legend.title = "Mean depth")
dev.off()



#Ecoregions distribution
png( here::here("figures", "Sites in ACP by ecoregion.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
ind.p <- fviz_pca_ind(pca, geom = "point",
                      pointshape=21, 
                      fill.ind = NCP_site$SiteEcoregion)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis of NCPs",
              subtitle = "RLS data set",
              caption = "Source: factoextra",
              xlab = "PC1", ylab = "PC2",
              legend.title = "Ecoregion", legend.position = "top",
              ggtheme = theme_gray(), palette = scico(58, palette = "roma")
)
dev.off()


png( here::here("figures", "biplot sites_NCP with Elipse by ecoregion.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_biplot (pca,
                 geom="point", pointshape=21,
                 col.ind = NCP_site$SiteEcoregion,
                 fill.ind = NCP_site$SiteEcoregion,
                 palette = scico(58, palette = "roma"),
                 addEllipses = T, label = "var",
                 col.var = "black", repel = TRUE,
                 legend.title = "Ecoregion")
dev.off()





##-------------Compute PCA at surveys scale-------------
rownames(NCP_surveys) <- NCP_surveys$SurveyID

pca <- FactoMineR::PCA(NCP_surveys[,-c(1:7)], scale.unit = T)
print(pca)
summary(pca)

png( here::here("figures", "surveys_scale", "surveys_scale", "Contribution of dimensions in variances.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 30))
dev.off()

png( here::here( "figures", "surveys_scale", "surveys_scale", "PCA_NCP_axes1-2.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte 
) #give the quality of representatin of each variables in the 2 first dimensions
dev.off()

png( here::here( "figures", "surveys_scale", "surveys_scale", "PCA_NCP_axes1-3.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_var(pca, col.var = "cos2",
             axe=c(1,3),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
dev.off()

png( here::here("figures", "surveys_scale", "surveys_scale", "PCA_NCP_axes1-4.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_var(pca, col.var = "cos2",
             axe=c(1,4),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
dev.off()




var <- get_pca_var(pca)
png( here::here("figures", "surveys_scale", "surveys_scale", "Contribution of NCP in dimensions.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
corrplot(var$contrib, is.corr=FALSE)    
dev.off()

png( here::here("figures", "surveys_scale", "surveys_scale", "Surveys in ACP by contry.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
ind.p <- fviz_pca_ind(pca, geom = "point",
                      pointshape=21, 
                      fill.ind = NCP_surveys$SiteCountry)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis of NCPs",
              subtitle = "RLS data set",
              caption = "Source: factoextra",
              xlab = "PC1", ylab = "PC2",
              legend.title = "Country", legend.position = "top",
              ggtheme = theme_gray(), palette = rainbow(38)
)
dev.off()

png( here::here("figures", "surveys_scale", "surveys_scale", "biplot surveys_NCP.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_biplot (pca,
                 geom="point", pointshape=21,
                 col.ind = NCP_surveys$SiteCountry,
                 fill.ind = NCP_surveys$SiteCountry,
                 palette = rainbow(38),
                 addEllipses = F, label = "var",
                 col.var = "black", repel = TRUE,
                 legend.title = "Country")
dev.off()

png( here::here("figures", "surveys_scale", "surveys_scale", "biplot surveys_NCP with Elipse by contry.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_biplot (pca,
                 geom="point", pointshape=21,
                 col.ind = NCP_surveys$SiteCountry,
                 fill.ind = NCP_surveys$SiteCountry,
                 palette = rainbow(38),
                 addEllipses = T, label = "var",
                 col.var = "black", repel = TRUE,
                 legend.title = "Country")
dev.off()

png( here::here("figures", "surveys_scale", "surveys_scale", "biplot surveys_NCP by temperature.png"),  width = 30, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_biplot(pca, 
                geom="point", pointshape=21,
                fill.ind = NCP_surveys$SiteMeanSST,
                col.ind = NCP_surveys$SiteMeanSST,
                gradient.cols = c("blue", "yellow", "red"),
                col.var = "black", repel = TRUE,
                legend.title = "Mean temperature")
dev.off()

