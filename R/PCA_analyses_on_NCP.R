################################################################################
##
## Makes PCA analyses on all NCPs
##
## PCA_analyses_on_NCP.R
##
## 03/01/2023
##
## Ulysse Flandrin
##
################################################################################

#-----------------Loading packages-------------------
pkgs <- c("here", "tidyverse","FactoMineR", "tibble", "questionr", "corrplot",
          "factoextra", "ggpubr", "scico", "RColorBrewer", "plotly", "fishualize", 
          "ggplot2", "patchwork", "colormap", "grDevices", "ggnewscale", "sf",
          "rdacca.hp")
nip <- pkgs[!(pkgs %in% installed.packages())]
nip <- lapply(nip, install.packages, dependencies = TRUE)
ip   <- unlist(lapply(pkgs, require, character.only = TRUE, quietly = TRUE))

rm(list=ls())

##-------------loading data-------------
load(here::here("outputs","all_NCP_site.Rdata"))
#benthic_imputed <- read.csv(here::here("data_raw", "source", "RLS_benthic_imputed.txt"), sep= " ")
load(here::here("biodiversity", "outputs", "occurrence_matrix_family_survey_relative_biomass.Rdata"))
load(here::here("data","metadata_surveys.Rdata"))

NCP_site_clean <- subset(NCP_site, select = -c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                             SiteMeanSST, SiteLatitude, SiteLongitude,
                             HDI, gravtot2, MarineEcosystemDependency,
                             coral_imputation))

##-------------Correlations between NCPs-------------
## With Biomass
plot_correlation <- function(x,y,i){  
  ggplot() +
    geom_point(aes(y = NCP_site[,y][[y]], x = x),
               color = col[i], alpha = 0.6, size = 1) +
    theme_bw() +
    labs(x = x_title, y = y) +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black"), 
          axis.title = element_text(size = 17))
}

x<- NCP_site$Btot
x_title = "total Biomass"
col <- fishualize::fish(n = ncol(NCP_site_clean), option = "Ostracion_whitleyi", begin = 0, end = 0.8)

plots <- lapply( 1:ncol(NCP_site_clean), function(i){
            plot_correlation(x,colnames(NCP_site_clean)[i],i)
        })

plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] + 
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here("outputs", "figures","NCP_correlation_with_biomass.png"), plot, width = 22, height =14 )



## With biodiversity
x<- NCP_site$taxo_richness
x_title = "taxonomic richness"
col <- fish(n = ncol(NCP_site_clean), option = "Ostracion_whitleyi", begin = 0, end = 0.8)

plots <- lapply( 1:ncol(NCP_site_clean), function(i){
            plot_correlation(x,colnames(NCP_site_clean)[i],i)
})

plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] + 
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]]+ plots[[27]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here("outputs", "figures","NCP_correlation_with_biodiversity.png"), plot, width = 22, height =14 )

#### NCPs distribution
col <- fishualize::fish(n = 27, option = "Ostracion_whitleyi", begin = 0, end = 0.8)
names(col) <- colnames(NCP_site_clean)

plots <- lapply(colnames(NCP_site_clean),function(i){
    p <- ggplot(NCP_site_clean) +
          aes(x = NCP_site_clean[,i][[i]]) +
          geom_histogram(bins = 40L,
                         fill = col[i][[1]],
                         col = "black") +
          labs(title = i) +
          xlab("") + ylab("") +
          theme_minimal() +
          theme( legend.position = "none")
  
   p
})

all_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] +
  plots[[8]] + plots[[9]] +plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] +
  plots[[15]] +  plots[[16]] + plots[[17]] + plots[[18]] + plots[[19]] + plots[[20]] + plots[[21]] +
  plots[[22]] +  plots[[23]] + plots[[24]] + plots[[25]] + plots[[26]] + plots[[27]] +
  theme(axis.title.y = element_text(margin = margin(r = -100, unit = "pt"))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here("outputs", "figures","NCP_distribution.png"), all_plot, width = 22, height =14 )


### Correlation test
names <- cor <- p_val <- c()
for( i in colnames(NCP_site_clean)){
  for(j in colnames(NCP_site_clean)){
    correlation <- stats::cor.test(get("NCP_site_clean")[[i]], get("NCP_site_clean")[[j]] )
  names <- c(names, paste0(i,"-",j))
  cor <- c(cor, correlation[["estimate"]][["cor"]])
  p_val <- c(p_val, correlation[["p.value"]])}
}

correlations_between_NCPs <- data.frame(name=names, cor = cor, p_val=p_val)
signif_correlations__NCPs <- dplyr::filter(correlations_between_NCPs, abs(cor) > 0.5)


#### Corr-matrix for all NCPs
M <- cor(NCP_site_clean)
png(filename = here("outputs", "figures","corr_matrix.png"), 
    width= 40, height = 30, units = "cm", res = 1000)
# ## circle + black number
corrplot(M, order = 'AOE', type = 'lower', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r')
corrplot(M, add = TRUE, type = 'upper', method = 'number', order = 'AOE', insig = 'p-value',
           diag = FALSE, tl.pos = 'n', cl.pos = 'n')

# corrplot(M, p.mat = testRes$p, insig = 'p-value')
# corrplot(M, order = 'hclust', type = 'lower', tl.pos = 'tp', tl.srt = 60, cl.pos = 'r')
# corrplot(M, order = 'hclust', add= T, type= 'upper', p.mat = testRes$p, insig = 'p-value', 
#          tl.pos= 'n', cl.pos = 'n')
dev.off() 


##-------------computing pca-------------
rownames(NCP_site_clean) <- NCP_site$SiteCode

NCP_site_for_pca <- subset(NCP_site, select = 
    c( Productivity, funct_distinctiveness, Vitamin_A_C, ED_Mean, aragonite, monohydrocalcite,
    recycling_N, recycling_P, amorphous_carbonate, low_mg_calcite, high_mg_calcite,
    biom_lowTL, biom_mediumTL, biom_highTL, fishery_biomass, iucn_species, 
    elasmobranch_diversity, taxo_richness, phylo_entropy, aesthe_survey, funct_entropy,
    Omega_3_C, Selenium_C, Iron_C, Calcium_C, Zinc_C))
# non_correlated_NCP_site_for_pca <- subset(NCP_site_clean, 
#                             select = c(
#                              #independents NCP
#                              Productivity, funct_distinctiveness, Vitamin_A_C, ED_Mean, aragonite, monohydrocalcite,
#                              Btot, # correlated with biomass: Btot, recycling_N, recycling_P, amorphous_carbonate, low_mg_calcite, 
#                              # high_mg_calcite, biom_lowTL, biom_mediumTL, biom_highTL, fishery_biomass
#                              iucn_species, # iucn_species, elasmobranch_diversity
#                              taxo_richness, # taxo_richness, phylo_entropy, aesthe_survey
#                              funct_entropy,
#                              Omega_3_C, # Omega_3_C, Selenium_C
#                              Iron_C # Iron_C, Calcium_C, Zinc_C
#                            ))

pca <- FactoMineR::PCA(NCP_site_for_pca, scale.unit = T, graph=F, ncp=10)

summary(NCP_site$SiteCountry)

####----------plot -------------

##---- PCA with all NCPs at the community scale ------
### Number of dimension
#### Variance explained by the 10 first dimensions
png(filename = here("outputs", "figures","barplot_variance_explained_all_NCP.png"), 
    width= 15, height = 10, units = "cm", res = 1000)
factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 30), barcolor = "darkred", barfill = "darkred", ncp=15)
dev.off()

#### Contribution of NCP in dimensions
var <- get_pca_var(pca)
png(filename = here("outputs", "figures","contribution_NCP_in_dimensions.png"), 
    width= 12, height = 15, units = "cm", res = 1000)
corrplot(var$contrib, is.corr=FALSE)    
dev.off()


#### PCA in the 2 first dimensions, with representation quality ($cos^{2}$) of each variables
png(filename = here("outputs", "figures","PCA_all_NCP.png"), 
    width= 20, height = 17, units = "cm", res = 1000)
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE 
)
dev.off()

png(filename = here("outputs", "figures","PCA_all_NCP_cos2_0.4.png"), 
    width= 20, height = 17, units = "cm", res = 1000)
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c( "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(cos2 = 0.4)
)
dev.off() 



## Classify variables in Nature for Nature (NN) and Nature for Society (NS)
var <- get_pca_var(pca)
grp <- c(recycling_N="NN", recycling_P="NN",Productivity="NS",taxo_richness="NN", funct_entropy="NN",
         funct_distinctiveness="NN", Selenium_C="NS", Zinc_C="NS", Omega_3_C="NS", Calcium_C="NS",
         Iron_C="NS", Vitamin_A_C="NS", phylo_entropy="NN", ED_Mean="NN", aesthe_survey="NS", iucn_species="NN",
         elasmobranch_diversity="NN", low_mg_calcite="NN", high_mg_calcite="NN", aragonite="NN",
         monohydrocalcite="NN", amorphous_carbonate="NN", biom_lowTL="NN", biom_mediumTL="NN",
         biom_highTL="NN", fishery_biomass="NS") # /!\ the order matter


grp <- as.factor(grp)

png(filename = here("outputs", "figures","PCA_NS_NN_axes1-2.png"), 
    width= 20, height = 17, units = "cm", res = 1000)
fviz_pca_var(pca, col.var = grp, 
             palette = c("forestgreen", "dodgerblue3"),
             legend.title = "Nature Based Contributions")
dev.off()
#### PCA in dimensions 3 and 4 *(for variables with cos2 \> 0.2, and according the NCPs groups)*

# fviz_pca_var(pca, col.var = "cos2",
#              axe=c(3,4),
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE,
#              select.var = list(cos2 = 0.2)
# )

png(filename = here("outputs", "figures","PCA_NS_NN_axes3-4.png"), 
    width= 20, height = 17, units = "cm", res = 1000)
fviz_pca_var(pca, col.var = grp, axe= c(3,4), repel = TRUE,
             palette = c("forestgreen", "dodgerblue3"),
             legend.title = "Nature Based Contributions",
             select.var = list(cos2 = 0))
dev.off() 


#### PCA in dimensions 5 and 6 *(for variables with cos2 \> 0.2, and according the NCPs groups)*
png(filename = here("outputs", "figures","PCA_cos2_0.2_axes5-6.png"), 
    width= 20, height = 17, units = "cm", res = 1000)
fviz_pca_var(pca, col.var = "cos2",
             axe=c(5,6),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(cos2 = 0.2)
)
dev.off()

fviz_pca_var(pca, col.var = grp, axe= c(5,6), repel = TRUE,
             palette = c("forestgreen", "dodgerblue3"),
             legend.title = "Nature Based Contributions",
             select.var = list(cos2 = 0))


#### PCA in dimensions 7 and 8 *(for variables with cos2 \> 0.1)*
fviz_pca_var(pca, col.var = "cos2",
             axe=c(7,8),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(cos2 = 0.1)
)



#------Categories as Diaz 2022-------
### Dimensions 1 and 2
var_3cat <- get_pca_var(pca)
grp_3cat <- c(recycling_N="regulating", recycling_P="regulating",Productivity="material",
              taxo_richness="nature", funct_entropy="regulating", funct_distinctiveness="nature",
              Selenium_C="material", Zinc_C="material", Omega_3_C="material",
              Calcium_C="material",Iron_C="material", Vitamin_A_C="material", 
              phylo_entropy="regulating", ED_Mean="nature", aesthe_survey="non material",
              iucn_species="nature", elasmobranch_diversity="nature", low_mg_calcite="regulating",
              high_mg_calcite="regulating", aragonite="regulating", monohydrocalcite="regulating",
              amorphous_carbonate="regulating", biom_lowTL="regulating", biom_mediumTL="regulating",
              biom_highTL="regulating", fishery_biomass="material") 

grp_3cat <- as.factor(grp_3cat)

png(filename = here("outputs", "figures","PCA_Diaz_categories.png"), 
    width= 23, height = 20, units = "cm", res = 1000)
fviz_pca_var(pca, col.var = grp_3cat,
             palette = c( "dodgerblue3","forestgreen", "darkred", "grey20"),
             legend.title = "Nature Assets")
dev.off()


### study categories separately  
#### Regulating NCP
regulating <- names(grp_3cat)[ grp_3cat == "regulating" ]
NCP_regulating <- NCP_site[,regulating]
rownames(NCP_regulating) <- rownames(NCP_site)

pca_regulating <- FactoMineR::PCA(NCP_regulating, scale.unit = T, graph=F, ncp=10)

hist <- factoextra::fviz_eig(pca_regulating, addlabels = TRUE, ylim = c(0, 45), barfill = "grey20", barcolor = "grey20")

# var <- get_pca_var(pca_regulating)
# corrplot(var$contrib, is.corr=FALSE)    
pca_plot <- fviz_pca_var(pca_regulating, col.var = "grey20",
                         repel = TRUE )

png(filename = here("outputs", "figures","PCA_regulating_NCP.png"), 
    width= 30, height = 17, units = "cm", res = 1000)
gridExtra::grid.arrange(hist, pca_plot, nrow = 1, widths= c(1,2))
dev.off()


#### material NCP
material <- names(grp_3cat)[ grp_3cat == "material" ]
NCP_material <- NCP_site[,material]
rownames(NCP_material) <- rownames(NCP_site)

pca_material <- FactoMineR::PCA(NCP_material, scale.unit = T, graph=F, ncp=10)

hist <- factoextra::fviz_eig(pca_material, addlabels = TRUE, ylim = c(0, 40), barfill = "dodgerblue3", barcolor = "dodgerblue3")

pca_plot <- fviz_pca_var(pca_material, col.var = "dodgerblue3",
                         repel = TRUE)

png(filename = here("outputs", "figures","PCA_regulating_NCP.png"), 
    width= 30, height = 17, units = "cm", res = 1000)
gridExtra::grid.arrange(hist, pca_plot, nrow = 1, widths= c(1,2))
dev.off()


#### Nature intrinsic value
nature <- names(grp_3cat)[ grp_3cat == "nature" ]
NCP_nature <- NCP_site[,nature]
rownames(NCP_nature) <- rownames(NCP_site)
pca_nature <- FactoMineR::PCA(NCP_nature, scale.unit = T, graph=F, ncp=10)

hist <- factoextra::fviz_eig(pca_nature, addlabels = TRUE, ylim = c(0, 35), barfill = "forestgreen", barcolor = "forestgreen")

pca_plot <- fviz_pca_var(pca_nature, col.var = "forestgreen",
                         repel = TRUE)

png(filename = here("outputs", "figures","PCA_nature_Diaz_NCP.png"), 
    width= 30, height = 17, units = "cm", res = 1000)
gridExtra::grid.arrange(hist, pca_plot, nrow = 1, widths= c(1,2))
dev.off()


#### Dimensions 3 and 4
png(filename = here("outputs", "figures","PCA_Diaz_categories_dim3-4.png"), 
    width= 23, height = 20, units = "cm", res = 1000)
fviz_pca_var(pca, col.var = grp_3cat,
             axes=c(3,4),
             palette = c( "dodgerblue3","forestgreen", "darkred", "grey20"),
             legend.title = "Nature Based Contributions")
dev.off()




##----- Sites distribution-----
### By countries
ind.p <- fviz_pca_ind(pca, geom = "point",
                      pointshape=21,
                      fill.ind = NCP_site$SiteCountry)

png(filename = here("outputs", "figures","PCA_countries.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
ggpubr::ggpar(ind.p,
              title = "Principal Component Analysis of NCPs",
              subtitle = "RLS data set",
              caption = "Source: factoextra",
              xlab = "PC1", ylab = "PC2",
              legend.title = "Country", legend.position = "top",
              font.legend = c(8, "plain", "black"),
              ggtheme = theme_gray(), palette = scico(38, palette = "roma"))
dev.off() 


#### By countries, with $cos2 > 0.4$ variables:
png(filename = here("outputs", "figures","PCA_countries_cos2_0.4.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
fviz_pca_biplot (pca,
                 select.var = list(cos2 = 0.4),
                 geom="point", pointshape=21,
                 col.ind = NCP_site$SiteCountry,
                 fill.ind = NCP_site$SiteCountry,
                 palette = scico(38, palette = "roma"),
                 addEllipses = F, label = "var",
                 col.var = "black", repel = TRUE,
                 legend.title = "Country",
                 font.legend = c(7, "plain", "black"))
dev.off()

#### clustering countries, with $cos2 > 0.4$ variables:
fviz_pca_biplot (pca,
                 select.var = list(cos2 = 0.4),
                 geom="point", pointshape=21,
                 col.ind = NCP_site$SiteCountry,
                 fill.ind = NCP_site$SiteCountry,
                 palette = scico(38, palette = "roma"),
                 addEllipses = T, label = "var",
                 col.var = "black", repel = TRUE,
                 legend.title = "Country",
                 font.legend = c(6, "plain", "black")
)



### Temperature pattern
png(filename = here("outputs", "figures","PCA_temperature_pattern.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
fviz_pca_biplot(pca,
                select.var = list(cos2 = 0.4),
                geom="point", pointshape=21,
                fill.ind = NCP_site$SiteMeanSST,
                col.ind = NCP_site$SiteMeanSST,
                gradient.cols = rev(brewer.pal(10,name="RdYlBu")),
                col.var = "black", repel = TRUE,
                legend.title = "Mean temperature")
dev.off()

df <- data.frame(x = NCP_site$SiteMeanSST, y = pca$ind$coord[,"Dim.2"])

png(filename = here("outputs", "figures","PCA_temp_according_Dim2.png"), 
    width= 15, height = 12, units = "cm", res = 1000)
ggplot(df, aes(x , y , alpha = 0.5)) +
  geom_point()+
  geom_smooth() +
  labs(x = "Mean SST", y = "Dimension  in global PCA", title = "Importance of SST in PCA's axe 2") +
  theme_minimal()+
  theme(legend.position = 'none')
dev.off()



### Gravity pattern:  log10(gravtot2)
png(filename = here("outputs", "figures","PCA_gravity_pattern.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
fviz_pca_biplot(pca,
                axes=c(1,2),
                select.var = list(cos2 = 0.4),
                geom="point", pointshape=21,
                fill.ind = log10(NCP_site$gravtot2),
                col.ind = log10(NCP_site$gravtot2),
                gradient.cols = rev(grDevices::colorRampPalette(c("darkred", "white", "forestgreen"))(10)),
                col.var = "black", repel = TRUE,
                invisible= "var",
                legend.title = "gravity impact")
dev.off()

plot(log(NCP_site$gravtot2) ~ pca$ind$coord[,"Dim.1"])

fviz_pca_biplot(pca,
                axes=c(1,3),
                select.var = list(cos2 = 0.4),
                geom="point", pointshape=21,
                fill.ind = log10(NCP_site$gravtot2),
                col.ind = log10(NCP_site$gravtot2),
                gradient.cols = rev(grDevices::colorRampPalette(c("darkred", "white", "forestgreen"))(10)),
                col.var = "black", repel = TRUE,
                invisible= "var",
                legend.title = "gravity impact")





### Country HDI
png(filename = here("outputs", "figures","PCA_HDI_pattern.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
fviz_pca_biplot(pca,
                axes= c(1,2),
                select.var = list(cos2 = 0.4),
                geom="point", pointshape=21,
                fill.ind = NCP_site$HDI,
                col.ind = NCP_site$HDI,
                gradient.cols = rev(grDevices::colorRampPalette(c("darkred", "white", "forestgreen"))(10)),
                col.var = "black", repel = TRUE,
                invisible= "var",
                legend.title = "Country HDI")
dev.off()

### Marine Ecosystem Dependency
png(filename = here("outputs", "figures","PCA_MED_pattern.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
fviz_pca_biplot(pca,
                axes= c(1,2),
                select.var = list(cos2 = 0.4),
                geom="point", pointshape=21,
                fill.ind = NCP_site$MarineEcosystemDependency,
                col.ind = NCP_site$MarineEcosystemDependency,
                gradient.cols = rev(grDevices::colorRampPalette(c("darkred", "white", "forestgreen"))(10)),
                col.var = "black", repel = TRUE,
                invisible= "var",
                legend.title = "MarineEcosystemDependency")
dev.off()


### Imputed coral cover
png(filename = here("outputs", "figures","PCA_coral_cover_pattern.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
fviz_pca_biplot(pca,
                axes= c(1,2),
                select.var = list(cos2 = 0.4),
                geom="point", pointshape=21,
                fill.ind = NCP_site$coral_imputation,
                col.ind = NCP_site$coral_imputation,
                gradient.cols = rev(grDevices::colorRampPalette(c("darkred", "white", "forestgreen"))(10)),
                col.var = "black", repel = TRUE,
                invisible= "var",
                legend.title = "Imputed coral cover")
dev.off()


### Depth pattern
png(filename = here("outputs", "figures","PCA_depth_pattern.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
fviz_pca_biplot(pca,
                geom="point", pointshape=21,
                fill.ind = NCP_site$SurveyDepth,
                col.ind = NCP_site$SurveyDepth,
                gradient.cols = brewer.pal(10,name="RdYlBu"),
                col.var = "black", repel = TRUE,
                legend.title = "Mean depth")
dev.off()


### Lattitude pattern
png(filename = here("outputs", "figures","PCA_lattitude_pattern.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
fviz_pca_biplot(pca,
                geom="point", pointshape=21,
                fill.ind = abs(NCP_site$SiteLatitude),
                col.ind = abs(NCP_site$SiteLatitude),
                gradient.cols = brewer.pal(10,name="RdYlBu"),
                col.var = "black", repel = TRUE,
                legend.title = "abs(latitude)")
dev.off()


### Ecoregions
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

#####biplot sites_NCP with Elipse by ecoregion
fviz_pca_biplot (pca,
                 geom="point", pointshape=21,
                 col.ind = NCP_site$SiteEcoregion,
                 fill.ind = NCP_site$SiteEcoregion,
                 palette = scico(58, palette = "roma"),
                 addEllipses = T, label = "var",
                 col.var = "black", repel = TRUE,
                 legend.title = "Ecoregion")


### Explore with plotly
site_coord_in_pca <- as.data.frame(pca$ind$coord) %>%
  cbind(NCP_site[,c("SiteCountry", "SiteMeanSST", "gravtot2")])

plotly::plot_ly(site_coord_in_pca,
        x= ~Dim.1, y= ~Dim.2, z= ~Dim.3,
        size = 10 ) %>%
  add_markers(color= ~SiteCountry)

plotly::plot_ly(site_coord_in_pca,
        x= ~Dim.1, y= ~Dim.2, z= ~Dim.3,
        size = 10
) %>%
  add_markers(color= ~SiteMeanSST)

plotly::plot_ly(site_coord_in_pca,
        x= ~Dim.1, y= ~Dim.2, z= ~Dim.3,
        size = 10
) %>%
  add_markers(color= ~log(gravtot2))


##------- Study NN and NS separately ------

### Nature contributions to Nature   
NN <- names(grp)[ grp=="NN" ]
NCP_NN <- NCP_site[,NN]
rownames(NCP_NN) <- rownames(NCP_site)

pca <- FactoMineR::PCA(NCP_NN, scale.unit = T, graph=F, ncp=10)

png(filename = here("outputs", "figures","barplot_variance_explained_by_NN.png"), 
    width= 15, height = 10, units = "cm", res = 1000)
factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 35), barfill = "forestgreen", barcolor = "forestgreen")
dev.off()

var <- get_pca_var(pca)
png(filename = here("outputs", "figures","contribution_NN_in_dimensions.png"), 
    width= 17, height = 15, units = "cm", res = 1000)
corrplot(var$contrib, is.corr=FALSE)
dev.off()

#### NN variations in the 2 first dimensions

png(filename = here("outputs", "figures","PCA_NN_only.png"), 
    width= 15, height = 15, units = "cm", res = 1000)
fviz_pca_var(pca, col.var = "forestgreen",
             repel = TRUE)
dev.off()


png(filename = here("outputs", "figures","PCA_countries_NN_only.png"), 
    width= 30, height = 20, units = "cm", res = 1000)
fviz_pca_biplot (pca,
                 select.var = list(cos2 = 0.1),
                 geom="point", pointshape=21,
                 col.ind = "black",
                 fill.ind = NCP_site$SiteCountry,
                 palette = scico(38, palette = "roma"),
                 addEllipses = F, label = "var",
                 col.var = "black", repel = TRUE,
                 legend.title = "Country",
                 font.legend = c(6, "plain", "black")
)
dev.off()

#### NN PCA in dimensions 3 and 4 *(for variables with cos2 \> 0.1)*

fviz_pca_var(pca, col.var = "cos2",
             axe=c(3,4),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(cos2 = 0.1)
)

#### Surveys positions in PC1

NN_PC1_surveys <- pca$ind$coord[,1]
summary(NN_PC1_surveys)

NN_PC1_surveys <- data.frame(NN_PC1=NN_PC1_surveys)
rownames(NN_PC1_surveys) <- rownames(NCP_site)

## NN
histo_NN <- ggplot(NN_PC1_surveys, aes(x = NN_PC1)) +
  geom_histogram( bins = 60, aes(fill=..x..), color="grey50", size=0.3)+
  scale_fill_gradient2(low = "black", mid = "white", high = "forestgreen", guide = "colourbar")+
  labs(x = "Nature for nature PC1 site's evaluation")+
  geom_vline(xintercept = 0, linetype = 2, col= "black")+
  theme_minimal()
histo_NN

### Nature contributions to Society

```{r, echo=FALSE, warning=FALSE}
NS <- names(grp)[ grp=="NS" ]
#NS <- NS[1:8] # enlever fishery biomass ################################
NCP_NS <- NCP_site[,NS]
rownames(NCP_NS) <- rownames(NCP_site)

pca <- FactoMineR::PCA(NCP_NS, scale.unit = T, graph=F, ncp=10)

factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 35), barfill = "dodgerblue3", barcolor = "dodgerblue3")

var <- get_pca_var(pca)
corrplot(var$contrib, is.corr=FALSE)
```

#### Nature for Society variations in the 2 first dimensions 

```{r, echo=TRUE, fig.height=9, fig.width=9}
fviz_pca_var(pca, col.var = "dodgerblue3",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE
)
```

#### NS PCA in dimensions 3 and 4 *(for variables with cos2 \> 0.1)*

```{r, echo=FALSE, fig.height=7, fig.width=7}
fviz_pca_var(pca, col.var = "cos2",
             axe=c(3,4),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(cos2 = 0.1)
)
```

#### Surveys positions in PC1

```{r}
NS_PC1_surveys <- pca$ind$coord[,1]
summary(NS_PC1_surveys)
```

```{r, echo=FALSE}
NN_NS_PC1_surveys <- data.frame(SiteCode = NCP_site$SiteCode ,
                                NN_PC1=NN_PC1_surveys,
                                NS_PC1=NS_PC1_surveys)
rownames(NN_NS_PC1_surveys) <- rownames(NCP_site)

##NS
histo_NS <- ggplot(NN_NS_PC1_surveys, aes(x = NS_PC1)) +
  geom_histogram( bins = 60, aes(fill=..x..), color="grey30", size=0.3)+
  scale_fill_gradient2(low = "black", mid = "white", high = "dodgerblue3", guide = "colourbar")+
  labs(x = "Nature to People PC1 site's evaluation")+
  geom_vline(xintercept = 0, linetype = 2, col= "black")+
  theme_minimal()
histo_NS

```

## NN against NS

```{r, warning=FALSE, message=FALSE}
##NS~NN
#Define colors by quarter

#up right
quarter_up_right <- NN_NS_PC1_surveys |>
  filter(NN_PC1 > 0 & NS_PC1 > 0) |>
  mutate(product_u_r =  pnorm(NN_PC1 * NS_PC1) )

pal <- grDevices::colorRampPalette(c("white", "darkred")) (nrow(quarter_up_right))

quarter_up_right <- quarter_up_right |>
  dplyr::arrange(product_u_r) |>
  dplyr::mutate( cols = pal)

#up left
quarter_up_left <- NN_NS_PC1_surveys |>
  filter(NN_PC1 < 0 & NS_PC1 > 0) |>
  mutate(product_u_l =  pnorm(-NN_PC1 * NS_PC1) )

pal <- grDevices::colorRampPalette(c("white", "dodgerblue3")) (nrow(quarter_up_left))

quarter_up_left <- quarter_up_left |>
  dplyr::arrange(product_u_l) |>
  dplyr::mutate( cols = pal)

#down right
quarter_down_right <- NN_NS_PC1_surveys |>
  filter(NN_PC1 > 0 & NS_PC1 < 0) |>
  mutate(product_d_r =  pnorm(NN_PC1 * -NS_PC1))

pal <- grDevices::colorRampPalette(c("white", "forestgreen")) (nrow(quarter_down_right))

quarter_down_right <- quarter_down_right |>
  dplyr::arrange(product_d_r) |>
  dplyr::mutate( cols = pal)

#down left
quarter_down_left <- NN_NS_PC1_surveys |>
  filter(NN_PC1 < 0 & NS_PC1 < 0) |>
  mutate(product_d_l =  pnorm(NN_PC1 * NS_PC1))

pal <- grDevices::colorRampPalette(c("white", "grey20")) (nrow(quarter_down_left))

quarter_down_left <- quarter_down_left |>
  dplyr::arrange(product_d_l) |>
  dplyr::mutate( cols = pal)


#Merge the 4 quarters
NN_NS_with_col <- quarter_down_left |>
  full_join(quarter_down_right) |>
  full_join(quarter_up_left) |>
  full_join(quarter_up_right)



library(ggplot2)
NN_NS_plot <- ggplot(NN_NS_with_col, aes( y= NS_PC1, x = NN_PC1) ) +
  #up right quarter
  geom_point(data= dplyr::filter(NN_NS_with_col, is.na(product_u_r)==F),
             size = 2,
             aes(colour= product_u_r, alpha = product_u_r)) +
  scale_colour_gradient("product_u_r",
                        low = "white", high="darkred",
                        na.value=NA) +
  
  ggnewscale::new_scale("colour") +
  
  #up left quarter
  geom_point(data= dplyr::filter(NN_NS_with_col, is.na(product_u_l)==F),
             size = 2,
             aes(colour=product_u_l, alpha = product_u_l)) +
  scale_colour_gradient("product_u_l",
                        low = "white", high="dodgerblue3",
                        na.value=NA) +
  geom_point(aes( y= NS_PC1, x = NN_PC1),
             shape = 1, size = 2, stroke = 1,
             color= "black",
             data = filter(NN_NS_with_col ,
                           product_u_l > quantile(NN_NS_with_col$product_u_l, 0.97, na.rm = TRUE))) +
  
  
  ggnewscale::new_scale("colour") +
  
  #down right quarter
  geom_point(data= dplyr::filter(NN_NS_with_col, is.na(product_d_r)==F),
             size = 2,
             aes(colour=product_d_r, alpha = product_d_r)) +
  scale_colour_gradient("product_d_r",
                        low = "white", high="forestgreen",
                        na.value=NA) +
  geom_point(aes( y= NS_PC1, x = NN_PC1),
             shape = 1, size = 2, stroke = 1,
             color= "black",
             data = head(NN_NS_with_col[order(NN_NS_with_col$product_d_r, decreasing = T),], 15)) +
  
  ggnewscale::new_scale("colour") +
  
  #down left quarter
  geom_point(data= dplyr::filter(NN_NS_with_col, is.na(product_d_l)==F),
             size = 2,
             aes(colour= product_d_l, alpha = product_d_l)) +
  scale_colour_gradient("product_d_l",
                        low = "white", high="grey20",
                        na.value=NA) +
  
  scale_alpha_continuous(range = c(0, 0.8)) +
  
  
  #add lines
  geom_vline(xintercept = 0, linetype = 2, col = "black")+
  geom_hline(yintercept = 0, linetype = 2, col = "black")+
  labs( x=  "Nature to Nature", y = "Nature to People")+
  theme_minimal()+
  theme(legend.position = "none")

NN_NS_plot
#ggsave( here::here("outputs", "figures", "Sites in NN and NS PC1.jpg"), plot = NN_NS_plot, width=20, height = 15 )

#plot(NN_NS_with_col$NN_PC1, NN_NS_with_col$NS_PC1, col=NN_NS_with_col$cols, pch=19)

```

## Plot on map

### NS according to NN

```{r, eval=TRUE,echo=FALSE, warning=FALSE, message=FALSE, fig.height=7, fig.width=15}
#world coast shapefile
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))
#data
NN_NS_with_col <- NN_NS_with_col %>%  left_join(NCP_site[,c("SiteCode", "SiteLongitude", "SiteLatitude")])
#plot
map <- ggplot(NN_NS_with_col) +
  geom_sf(data = coast, color = NA, fill = "lightgrey") +
  
  #down left quarter
  geom_point(data= dplyr::filter(NN_NS_with_col, is.na(product_d_l)==F),
             size = 2, alpha = 0.5,
             aes(x = SiteLongitude, y = SiteLatitude,
                 colour= product_d_l, alpha = product_d_l)) +
  scale_colour_gradient("product_d_l",
                        low = "white", high="grey20",
                        na.value=NA) +
  ggnewscale::new_scale("colour") +
  
  #up left quarter
  geom_point(data= dplyr::filter(NN_NS_with_col, is.na(product_u_l)==F),
             size = 2, alpha = 0.5,
             aes(x = SiteLongitude, y = SiteLatitude,
                 colour=product_u_l, alpha = product_u_l)) +
  scale_colour_gradient("product_u_l",
                        low = "white", high="dodgerblue3",
                        na.value=NA) +
  geom_point(aes(x = SiteLongitude, y = SiteLatitude),
             shape = 1, size = 2, stroke = 1,
             color= "black",
             data = head(NN_NS_with_col[order(NN_NS_with_col$product_d_r, decreasing = T),], 8)) +
  
  ggnewscale::new_scale("colour") +
  
  #down right quarter
  geom_point(data= dplyr::filter(NN_NS_with_col, is.na(product_d_r)==F),
             size = 2, alpha = 0.5,
             aes(x = SiteLongitude, y = SiteLatitude,
                 colour=product_d_r, alpha = product_d_r)) +
  scale_colour_gradient("product_d_r",
                        low = "white", high="forestgreen",
                        na.value=NA) +
  geom_point(aes(x = SiteLongitude, y = SiteLatitude),
             shape = 1, size = 2, stroke = 1,
             color= "black",
             data = filter(NN_NS_with_col ,
                           product_d_r >= quantile(NN_NS_with_col$product_d_r, 0.95, na.rm = TRUE))) +
  
  ggnewscale::new_scale("colour") +
  
  #up right quarter
  geom_point(data= dplyr::filter(NN_NS_with_col, is.na(product_u_r)==F),
             size = 2,
             aes(x = SiteLongitude, y = SiteLatitude,
                 colour= product_u_r, alpha = product_u_r)) +
  scale_colour_gradient("product_u_r",
                        low = "white", high="darkred",
                        na.value=NA) +
  #add transparency
  scale_alpha_continuous(range = c(0, 1)) +
  
  
  coord_sf(ylim = c(-36, 31), expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
  theme_minimal()+
  labs(title = "Trade-offs in Nature Based Contribution worldwide",
       x="Longitude", y= "Latitude") +
  theme(legend.position = "none",
        plot.title = element_text(size=10, face="bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
  )
#ggsave( here::here("outputs", "figures", "world map with NN and NS PC1.jpg"), plot = map, width=15, height = 7 )
map
```

### NN worldwide

```{r,echo=FALSE, warning=FALSE, message=FALSE,fig.height=7, fig.width=15}
##NN
map <- ggplot(NN_NS_with_col) +
  geom_sf(data = coast, color = NA, fill = "lightgrey") +
  geom_point(aes(x = SiteLongitude, y = SiteLatitude,
                 color = NN_PC1, alpha= abs(NN_PC1),
                 size = 2)) +
  geom_point(aes(x = SiteLongitude, y = SiteLatitude),
             shape = 1, size = 2, stroke = 0.5,
             color= "black",
             data = head(NN_NS_with_col[order(NN_NS_with_col$product_d_r, decreasing = T),], 15)) +
  
  scale_colour_gradientn(colours = colorRampPalette(rev(c( "forestgreen","white", "grey30")))(1000)) +
  scale_alpha_continuous(range = c(0, 1)) +
  
  coord_sf(ylim = c(-36, 31), expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
  theme_minimal()+
  labs(title = "Nature for Nature",
       x="Longitude", y= "Latitude") +
  theme(
    plot.title = element_text(size=10, face="bold"),
  )
map
```

### NS worldwide

```{r,echo=FALSE, warning=FALSE,fig.height=7, fig.width=15}
##NS
map <- ggplot(NN_NS_with_col) +
  geom_sf(data = coast, color = NA, fill = "lightgrey") +
  geom_point(aes(x = SiteLongitude, y = SiteLatitude,
                 color = NS_PC1, alpha= abs(NS_PC1),
                 size = 2)) +
  geom_point(aes(x = SiteLongitude, y = SiteLatitude),
             shape = 1, size = 2, stroke = 1,
             color= "black",
             data = filter(NN_NS_with_col ,
                           NN_PC1 >= quantile(NN_PC1, 0.99, na.rm = TRUE))) +
  
  scale_colour_gradientn(colours = colorRampPalette(rev(c( "dodgerblue3","white", "grey30")))(1000)) +
  scale_alpha_continuous(range = c(0, 1)) +
  
  coord_sf(ylim = c(-36, 31), expand = FALSE) +
  scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
  theme_minimal()+
  labs(title = "Nature to People",
       x="Longitude", y= "Latitude") +
  theme(
    plot.title = element_text(size=10, face="bold"),
  )

map
```
