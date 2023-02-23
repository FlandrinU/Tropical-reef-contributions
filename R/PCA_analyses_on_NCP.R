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
# load(here::here("outputs","all_NCP_site.Rdata"))
# load(here::here("outputs","NCP_site_coral_reef.Rdata"))
# load(here::here("outputs","NCP_site_SST20.Rdata"))
# load(here::here("outputs","NCP_site_coral_5_imputed.Rdata"))
load(here::here("outputs","NCP_site_wo_australia.Rdata"))

NCP_site <- NCP_site_condition

#benthic_imputed <- read.csv(here::here("data_raw", "source", "RLS_benthic_imputed.txt"), sep= " ")
load(here::here("biodiversity", "outputs", "occurrence_matrix_family_survey_relative_biomass.Rdata"))
coast <- sf::st_read(here::here("data", "ShapeFiles coast", "GSHHS_h_L1.shp"))

##-------------Correlations between NCPs-------------

  ## Clean data
  NCP_site_clean <- subset(NCP_site, select = -c(SiteCode, SiteCountry, SiteEcoregion, SurveyDepth, 
                                                 SiteMeanSST, SiteLatitude, SiteLongitude,
                                                 HDI, MarineEcosystemDependency,
                                                 coral_imputation))
  

  #### NCPs distribution with log correction
  NCP_skewed_distribution <- c("Btot","recycling_N","recycling_P","Productivity",
                               "funct_distinctiveness","Omega_3_C","Calcium_C","Vitamin_A_C",
                               "phylo_entropy","ED_Mean", "iucn_species", "elasmobranch_diversity",
                               "low_mg_calcite", "high_mg_calcite", "aragonite", "monohydrocalcite",
                               "amorphous_carbonate", "biom_lowTL", "biom_mediumTL", "biom_highTL",
                               "fishery_biomass")
  
  NCP_log_transformed <- NCP_site_clean |>
    dplyr::mutate(across(.cols = all_of(NCP_skewed_distribution),
                     .fns = ~ .x +1 , .names = "{.col}")) |>      # Adds 1 to values to log transformed
    dplyr::mutate(across(.cols = all_of(NCP_skewed_distribution),
                         .fns = log10 , .names = "{.col}")) |>       # log(x+1) to avoid negative values
    dplyr::rename_with(.cols = all_of(NCP_skewed_distribution),
                       .fn = ~ paste0("log(", .x, ")"))
  
  
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
  
 


##-------------computing PCA-------------
NCP_to_transform <- c("Btot","recycling_N","recycling_P","Productivity",
                        "funct_distinctiveness","Omega_3_C","Calcium_C","Vitamin_A_C",
                        "phylo_entropy","ED_Mean", "iucn_species", "elasmobranch_diversity",
                        "low_mg_calcite", "high_mg_calcite", "aragonite", "monohydrocalcite",
                        "amorphous_carbonate", "biom_lowTL", "biom_mediumTL", "biom_highTL",
                        "fishery_biomass")

NCP_site <- NCP_site |>
  dplyr::mutate(across(.cols = all_of(NCP_to_transform),
                       .fns = ~ .x +1 , .names = "{.col}")) |>     
  dplyr::mutate(across(.cols = all_of(NCP_to_transform),
                       .fns = log10 , .names = "{.col}")) 


plot_PCA_NCP <- function(NCP_site){
  library(ggplot2)
  NCP_site_selected <- subset(NCP_site, select = 
    c(recycling_N, recycling_P,Productivity,taxo_richness, funct_entropy,
      funct_distinctiveness, Selenium_C, Zinc_C, Omega_3_C, Calcium_C,
      Iron_C, Vitamin_A_C, phylo_entropy, ED_Mean, aesthe_survey, iucn_species,
      elasmobranch_diversity, low_mg_calcite, high_mg_calcite, aragonite,
      monohydrocalcite, amorphous_carbonate, biom_lowTL, biom_mediumTL,
      biom_highTL, fishery_biomass, mean_TL, robustness))
  
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
  
  NCP_site_for_pca <- scale(NCP_site_selected)
  pca <- FactoMineR::PCA(NCP_site_for_pca, scale.unit = FALSE, graph=F, ncp=10) 
  
  summary(NCP_site$SiteCountry)
  
  ####----------plot -------------
  
  ##---- PCA with all NCPs at the community scale ------
  ### Number of dimension
  #### Variance explained by the 15 first dimensions
  png(filename = here::here("outputs", "figures","barplot_variance_explained_all_NCP.png"), 
      width= 15, height = 10, units = "cm", res = 1000)
  print(
    factoextra::fviz_eig(pca, choice = "variance",
                       addlabels = TRUE, ylim = c(0, 35), barcolor = "darkred", 
                       barfill = "darkred", ncp=15) )
  dev.off()
  
  
  eig <- factoextra::get_eig(pca)
  variance_explained <- data.frame(
    Dimensions = c(1:nrow(eig)), 
    cumulative_variance_explained = eig[,3],
    contribution_coefficient = eig[,2]/100)

  cumulative_variance <- ggplot(variance_explained) +
    aes(x = Dimensions, y = cumulative_variance_explained) +
    geom_step(colour = "darkred") +
    geom_hline(yintercept=70, linetype="dashed", color = "black")+
    ylim(0,100) +
    labs(
      x = "Dimensions",
      y = "Variance explained",
      title = "Cumulative variance explained by dimensions (threshold at 70%)"
    ) +
    theme_minimal()
  ggsave(filename = here::here("outputs", "figures","barplot_cumulative_variance_explained_all_NCP.png"), cumulative_variance, width = 15, height =10 )
  
  #### Contribution of NCP in dimensions
  var <- factoextra::get_pca_var(pca)
  contributions <- var$contrib
  for( i in 1:ncol(contributions)){ 
    contributions[,i] <- contributions[,i] * variance_explained$contribution_coefficient[i]}
  
  png(filename = here::here("outputs", "figures","contribution_NCP_in_dimensions.png"), 
      width= 12, height = 15, units = "cm", res = 1000)
  print( corrplot::corrplot(var$contrib, is.corr=FALSE) )
  dev.off()
  
  png(filename = here::here("outputs", "figures","contribution_NCP_in_total_variance.png"), 
      width= 12, height = 15, units = "cm", res = 1000)
  print( corrplot::corrplot(contributions, is.corr=FALSE) )
  dev.off()
  
  #### PCA in the 2 first dimensions, with representation quality ($cos^{2}$) of each variables
  png(filename = here::here("outputs", "figures","PCA_all_NCP.png"), 
      width= 20, height = 17, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE 
  ))
  dev.off()
  
  png(filename = here::here("outputs", "figures","PCA_all_NCP_cos2_0.4.png"), 
      width= 20, height = 17, units = "cm", res = 1000)
  print(factoextra::fviz_pca_var(pca, col.var = "cos2",
               gradient.cols = c( "#E7B800", "#FC4E07"),
               repel = TRUE,
               select.var = list(cos2 = 0.4)
  ))
  dev.off() 
  
  
  
  ## Classify variables in Nature for Nature (NN) and Nature for Society (NS)
  #var <- factoextra::get_pca_var(pca)
  grp <- c(recycling_N="NN", recycling_P="NN",Productivity="NS",taxo_richness="NN", funct_entropy="NN",
           funct_distinctiveness="NN", Selenium_C="NS", Zinc_C="NS", Omega_3_C="NS", Calcium_C="NS",
           Iron_C="NS", Vitamin_A_C="NS", phylo_entropy="NN", ED_Mean="NN", aesthe_survey="NS", iucn_species="NN",
           elasmobranch_diversity="NN", low_mg_calcite="NN", high_mg_calcite="NN", aragonite="NN",
           monohydrocalcite="NN", amorphous_carbonate="NN", biom_lowTL="NN", biom_mediumTL="NN",
           biom_highTL="NN", fishery_biomass="NS", mean_TL = "NN", 
           robustness = "NN") # /!\ the order matter
  
  
  grp <- as.factor(grp)
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes1-2.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp, 
               palette = c("forestgreen", "dodgerblue3"),
               legend.title = "Nature Based Contributions",
               repel = TRUE))
  dev.off()
  
  #### PCA in dimensions 3 and 4 *(for variables with cos2 \> 0.2, and according the NCPs groups)*
  
  # fviz_pca_var(pca, col.var = "cos2",
  #              axe=c(3,4),
  #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  #              repel = TRUE,
  #              select.var = list(cos2 = 0.2)
  # )
  
  png(filename = here::here("outputs", "figures","PCA_NS_NN_axes3-4.png"), 
      width= 30, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp, axe= c(3,4), repel = TRUE,
               palette = c("forestgreen", "dodgerblue3"),
               legend.title = "Nature Based Contributions",
               select.var = list(cos2 = 0)))
  dev.off() 
  
  
  #### PCA in dimensions 5 and 6 *(for variables with cos2 \> 0.2, and according the NCPs groups)*
  png(filename = here::here("outputs", "figures","PCA_cos2_0.2_axes5-6.png"), 
      width= 30, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = "cos2",
               axe=c(5,6),
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE,
               select.var = list(cos2 = 0.2)
  ))
  dev.off()
  
  factoextra::fviz_pca_var(pca, col.var = grp, axe= c(5,6), repel = TRUE,
               palette = c("forestgreen", "dodgerblue3"),
               legend.title = "Nature Based Contributions",
               select.var = list(cos2 = 0))
  
  
  #### PCA in dimensions 7 and 8 *(for variables with cos2 \> 0.1)*
  factoextra::fviz_pca_var(pca, col.var = "cos2",
               axe=c(7,8),
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE,
               select.var = list(cos2 = 0.1)
  )
  
  
  
  #------Categories as Diaz 2022-------
  ### Dimensions 1 and 2
  var_3cat <- factoextra::get_pca_var(pca)
  grp_3cat <- c(recycling_N="regulating", recycling_P="regulating",Productivity="material",
                taxo_richness="nature", funct_entropy="regulating", funct_distinctiveness="nature",
                Selenium_C="material", Zinc_C="material", Omega_3_C="material",
                Calcium_C="material",Iron_C="material", Vitamin_A_C="material", 
                phylo_entropy="regulating", ED_Mean="nature", aesthe_survey="non material",
                iucn_species="nature", elasmobranch_diversity="nature", low_mg_calcite="regulating",
                high_mg_calcite="regulating", aragonite="regulating", monohydrocalcite="regulating",
                amorphous_carbonate="regulating", biom_lowTL="regulating", biom_mediumTL="regulating",
                biom_highTL="regulating", fishery_biomass="material",
                mean_TL = "nature", robustness = "nature") 
  
  grp_3cat <- as.factor(grp_3cat)
  
  png(filename = here::here("outputs", "figures","PCA_Diaz_categories.png"), 
      width= 35, height = 25, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp_3cat,
               palette = c( "dodgerblue3","forestgreen", "darkred", "grey20"),
               legend.title = "Nature Assets",
               repel = TRUE))
  dev.off()
  
  
  ### study categories separately  
  #### Regulating NCP
  regulating <- names(grp_3cat)[ grp_3cat == "regulating" ]
  NCP_regulating <- NCP_site[,regulating]
  rownames(NCP_regulating) <- rownames(NCP_site)
  
  pca_regulating <- FactoMineR::PCA(NCP_regulating, scale.unit = T, graph=F, ncp=10)
  
  hist <- factoextra::fviz_eig(pca_regulating, addlabels = TRUE, ylim = c(0, 45), barfill = "grey20", barcolor = "grey20")
  
  # var <- (pca_regulating)
  # corrplot(var$contrib, is.corr=FALSE)    
  pca_plot <- factoextra::fviz_pca_var(pca_regulating, col.var = "grey20",
                           repel = TRUE )
  
  png(filename = here::here("outputs", "figures","PCA_regulating_NCP.png"), 
      width= 30, height = 17, units = "cm", res = 1000)
  print( gridExtra::grid.arrange(hist, pca_plot, nrow = 1, widths= c(1,2)) )
  dev.off()
  
  
  #### material NCP
  material <- names(grp_3cat)[ grp_3cat == "material" ]
  NCP_material <- NCP_site[,material]
  rownames(NCP_material) <- rownames(NCP_site)
  
  pca_material <- FactoMineR::PCA(NCP_material, scale.unit = T, graph=F, ncp=10)
  
  hist <- factoextra::fviz_eig(pca_material, addlabels = TRUE, ylim = c(0, 40), barfill = "dodgerblue3", barcolor = "dodgerblue3")
  
  pca_plot <- factoextra::fviz_pca_var(pca_material, col.var = "dodgerblue3",
                           repel = TRUE)
  
  png(filename = here::here("outputs", "figures","PCA_material_NCP.png"), 
      width= 30, height = 17, units = "cm", res = 1000)
  print( gridExtra::grid.arrange(hist, pca_plot, nrow = 1, widths= c(1,2)) )
  dev.off()
  
  
  #### Nature intrinsic value
  nature <- names(grp_3cat)[ grp_3cat == "nature" ]
  NCP_nature <- NCP_site[,nature]
  rownames(NCP_nature) <- rownames(NCP_site)
  pca_nature <- FactoMineR::PCA(NCP_nature, scale.unit = T, graph=F, ncp=10)
  
  hist <- factoextra::fviz_eig(pca_nature, addlabels = TRUE, ylim = c(0, 35), barfill = "forestgreen", barcolor = "forestgreen")
  
  pca_plot <- factoextra::fviz_pca_var(pca_nature, col.var = "forestgreen",
                           repel = TRUE)
  
  png(filename = here::here("outputs", "figures","PCA_nature_Diaz_NCP.png"), 
      width= 30, height = 17, units = "cm", res = 1000)
  print( gridExtra::grid.arrange(hist, pca_plot, nrow = 1, widths= c(1,2)) )
  dev.off()
  
  
  #### Dimensions 3 and 4
  png(filename = here::here("outputs", "figures","PCA_Diaz_categories_dim3-4.png"), 
      width= 23, height = 20, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = grp_3cat,
               axes=c(3,4),
               palette = c( "dodgerblue3","forestgreen", "darkred", "grey20"),
               legend.title = "Nature Based Contributions"))
  dev.off()
  
  
  
  
  ##----- Sites distribution-----
  ### By countries
  ind.p <- factoextra::fviz_pca_ind(pca, geom = "point",
                        pointshape=21,
                        fill.ind = NCP_site$SiteCountry)
  
  png(filename = here::here("outputs", "figures","PCA_countries.png"), 
      width= 30, height = 20, units = "cm", res = 1000)
  print( ggpubr::ggpar(ind.p,
                title = "Principal Component Analysis of NCPs",
                subtitle = "RLS data set",
                caption = "Source: factoextra",
                xlab = "PC1", ylab = "PC2",
                legend.title = "Country", legend.position = "top",
                font.legend = c(8, "plain", "black"),
                ggtheme = theme_gray(), palette = scico::scico(38, palette = "roma")))
  dev.off() 
  
  
  #### By countries, with $cos2 > 0.4$ variables:
  png(filename = here::here("outputs", "figures","PCA_countries_cos2_0.4.png"), 
      width= 30, height = 20, units = "cm", res = 1000)
  print( factoextra::fviz_pca_biplot (pca,
                   select.var = list(cos2 = 0.4),
                   geom="point", pointshape=21,
                   col.ind = NCP_site$SiteCountry,
                   fill.ind = NCP_site$SiteCountry,
                   palette = scico::scico(38, palette = "roma"),
                   addEllipses = F, label = "var",
                   col.var = "black", repel = TRUE,
                   legend.title = "Country",
                   font.legend = c(7, "plain", "black")))
  dev.off()
  
  #### clustering countries, with $cos2 > 0.4$ variables:
  factoextra::fviz_pca_biplot (pca,
                   select.var = list(cos2 = 0.4),
                   geom="point", pointshape=21,
                   col.ind = NCP_site$SiteCountry,
                   fill.ind = NCP_site$SiteCountry,
                   palette = scico::scico(38, palette = "roma"),
                   addEllipses = T, label = "var",
                   col.var = "black", repel = TRUE,
                   legend.title = "Country",
                   font.legend = c(6, "plain", "black")
  )
  
  
  
  ### Temperature pattern
  png(filename = here::here("outputs", "figures","PCA_temperature_pattern.png"), 
      width= 30, height = 20, units = "cm", res = 1000)
  print( factoextra::fviz_pca_biplot(pca,
                  select.var = list(cos2 = 0.4),
                  geom="point", pointshape=21,
                  fill.ind = NCP_site$SiteMeanSST,
                  col.ind = NCP_site$SiteMeanSST,
                  gradient.cols = rev(RColorBrewer::brewer.pal(10,name="RdYlBu")),
                  col.var = "black", repel = TRUE,
                  legend.title = "Mean temperature"))
  dev.off()
  
  df <- data.frame(x = NCP_site$SiteMeanSST, y = pca$ind$coord[,"Dim.2"])
  temp <- ggplot(df, aes(x , y , alpha = 0.5)) +
    geom_point()+
    geom_smooth() +
    labs(x = "Mean SST", y = "Dimension  in global PCA", title = "Importance of SST in PCA's axe 2") +
    theme_minimal()+
    theme(legend.position = 'none')
  ggsave(filename = here::here("outputs", "figures", "PCA_temp_according_Dim2.png"), temp, width = 15, height = 12)
  
  
  # ### Gravity pattern:  log10(gravtot2)
  # png(filename = here::here("outputs", "figures","PCA_gravity_pattern.png"), 
  #     width= 30, height = 20, units = "cm", res = 1000)
  # print( factoextra::fviz_pca_biplot(pca,
  #                 axes=c(1,2),
  #                 select.var = list(cos2 = 0.4),
  #                 geom="point", pointshape=21,
  #                 fill.ind = log10(NCP_site$gravtot2),
  #                 col.ind = log10(NCP_site$gravtot2),
  #                 gradient.cols = rev(grDevices::colorRampPalette(c("darkred", "white", "forestgreen"))(10)),
  #                 col.var = "black", repel = TRUE,
  #                 invisible= "var",
  #                 legend.title = "gravity impact"))
  # dev.off()
  # 
  # plot(log(NCP_site$gravtot2) ~ pca$ind$coord[,"Dim.1"])
  # 
  # factoextra::fviz_pca_biplot(pca,
  #                 axes=c(1,3),
  #                 select.var = list(cos2 = 0.4),
  #                 geom="point", pointshape=21,
  #                 fill.ind = log10(NCP_site$gravtot2),
  #                 col.ind = log10(NCP_site$gravtot2),
  #                 gradient.cols = rev(grDevices::colorRampPalette(c("darkred", "white", "forestgreen"))(10)),
  #                 col.var = "black", repel = TRUE,
  #                 invisible= "var",
  #                 legend.title = "gravity impact")
  # 
  
  
  
  
  ### Country HDI
  png(filename = here::here("outputs", "figures","PCA_HDI_pattern.png"), 
      width= 30, height = 20, units = "cm", res = 1000)
  print( factoextra::fviz_pca_biplot(pca,
                  axes= c(1,2),
                  select.var = list(cos2 = 0.4),
                  geom="point", pointshape=21,
                  fill.ind = NCP_site$HDI,
                  col.ind = NCP_site$HDI,
                  gradient.cols = rev(grDevices::colorRampPalette(c("darkred", "white", "forestgreen"))(10)),
                  col.var = "black", repel = TRUE,
                  invisible= "var",
                  legend.title = "Country HDI"))
  dev.off()
  
  ### Marine Ecosystem Dependency
  png(filename = here::here("outputs", "figures","PCA_MED_pattern.png"), 
      width= 30, height = 20, units = "cm", res = 1000)
  print( factoextra::fviz_pca_biplot(pca,
                  axes= c(1,2),
                  select.var = list(cos2 = 0.4),
                  geom="point", pointshape=21,
                  fill.ind = NCP_site$MarineEcosystemDependency,
                  col.ind = NCP_site$MarineEcosystemDependency,
                  gradient.cols = rev(grDevices::colorRampPalette(c("darkred", "white", "forestgreen"))(10)),
                  col.var = "black", repel = TRUE,
                  invisible= "var",
                  legend.title = "MarineEcosystemDependency"))
  dev.off()
  
  
  ### Imputed coral cover
  png(filename = here::here("outputs", "figures","PCA_coral_cover_pattern.png"), 
      width= 30, height = 20, units = "cm", res = 1000)
  print( factoextra::fviz_pca_biplot(pca,
                  axes= c(1,2),
                  select.var = list(cos2 = 0.4),
                  geom="point", pointshape=21,
                  fill.ind = NCP_site$coral_imputation,
                  col.ind = NCP_site$coral_imputation,
                  gradient.cols = rev(grDevices::colorRampPalette(c("darkred", "white", "forestgreen"))(10)),
                  col.var = "black", repel = TRUE,
                  invisible= "var",
                  legend.title = "Imputed coral cover"))
  dev.off()
  
  
  ### Depth pattern
  png(filename = here::here("outputs", "figures","PCA_depth_pattern.png"), 
      width= 30, height = 20, units = "cm", res = 1000)
  print( factoextra::fviz_pca_biplot(pca,
                  geom="point", pointshape=21,
                  fill.ind = NCP_site$SurveyDepth,
                  col.ind = NCP_site$SurveyDepth,
                  gradient.cols = RColorBrewer::brewer.pal(10,name="RdYlBu"),
                  col.var = "black", repel = TRUE,
                  legend.title = "Mean depth"))
  dev.off()
  
  
  ### Lattitude pattern
  png(filename = here::here("outputs", "figures","PCA_lattitude_pattern.png"), 
      width= 30, height = 20, units = "cm", res = 1000)
  print( factoextra::fviz_pca_biplot(pca,
                  geom="point", pointshape=21,
                  fill.ind = abs(NCP_site$SiteLatitude),
                  col.ind = abs(NCP_site$SiteLatitude),
                  gradient.cols = RColorBrewer::brewer.pal(10,name="RdYlBu"),
                  col.var = "black", repel = TRUE,
                  legend.title = "abs(latitude)") )
  dev.off()
  
  
  ### Ecoregions
  ind.p <- factoextra::fviz_pca_ind(pca, geom = "point",
                        pointshape=21,
                        fill.ind = NCP_site$SiteEcoregion)
  ggpubr::ggpar(ind.p,
                title = "Principal Component Analysis of NCPs",
                subtitle = "RLS data set",
                caption = "Source: factoextra",
                xlab = "PC1", ylab = "PC2",
                legend.title = "Ecoregion", legend.position = "top",
                ggtheme = theme_gray(), palette = scico::scico(58, palette = "roma")
  )
  
  #####biplot sites_NCP with Elipse by ecoregion
  factoextra::fviz_pca_biplot (pca,
                   geom="point", pointshape=21,
                   col.ind = NCP_site$SiteEcoregion,
                   fill.ind = NCP_site$SiteEcoregion,
                   palette = scico::scico(58, palette = "roma"),
                   addEllipses = T, label = "var",
                   col.var = "black", repel = TRUE,
                   legend.title = "Ecoregion")
  
  
  # ### Explore with plotly
  # site_coord_in_pca <- as.data.frame(pca$ind$coord) |>
  #   cbind(NCP_site[,c("SiteCountry", "SiteMeanSST", "gravtot2")])
  # 
  # plotly::plot_ly(site_coord_in_pca,
  #         x= ~Dim.1, y= ~Dim.2, z= ~Dim.3,
  #         size = 10 ) |>
  #   add_markers(color= ~SiteCountry)
  # 
  # plotly::plot_ly(site_coord_in_pca,
  #         x= ~Dim.1, y= ~Dim.2, z= ~Dim.3,
  #         size = 10
  # ) |>
  #   add_markers(color= ~SiteMeanSST)
  # 
  # plotly::plot_ly(site_coord_in_pca,
  #         x= ~Dim.1, y= ~Dim.2, z= ~Dim.3,
  #         size = 10
  # ) |>
  #   add_markers(color= ~log(gravtot2))
  
  
  ##------- Study NN and NS separately ------
  grp <- as.factor( c(recycling_N="NN", recycling_P="NN",Productivity="NS",taxo_richness="NN", funct_entropy="NN",
           funct_distinctiveness="NN", Selenium_C="NS", Zinc_C="NS", Omega_3_C="NS", Calcium_C="NS",
           Iron_C="NS", Vitamin_A_C="NS", phylo_entropy="NN", ED_Mean="NN", aesthe_survey="NS", iucn_species="NN",
           elasmobranch_diversity="NN", low_mg_calcite="NN", high_mg_calcite="NN", aragonite="NN",
           monohydrocalcite="NN", amorphous_carbonate="NN", biom_lowTL="NN", biom_mediumTL="NN",
           biom_highTL="NN", fishery_biomass="NS", mean_TL = "NN",
           robustness = "NN")) # /!\ the order matter
  
  library(ggplot2)
  
  ###------- NN: Nature contributions to Nature -------   
  NN <- names(grp)[ grp=="NN" ]
  NCP_NN <- NCP_site[,NN]
  rownames(NCP_NN) <- rownames(NCP_site)
  
  pca <- FactoMineR::PCA(NCP_NN, scale.unit = T, graph=F, ncp=10)
  
  png(filename = here::here("outputs", "figures","barplot_variance_explained_by_NN.png"), 
      width= 15, height = 10, units = "cm", res = 1000)
  print( factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 40),
                              barfill = "forestgreen", barcolor = "forestgreen") )
  dev.off()
  
  var <- factoextra::get_pca_var(pca)
  png(filename = here::here("outputs", "figures","contribution_NN_in_dimensions.png"), 
      width= 17, height = 15, units = "cm", res = 1000)
  print( corrplot::corrplot(var$contrib, is.corr=FALSE) )
  dev.off()
  
  #### NN variations in the 2 first dimensions
  
  png(filename = here::here("outputs", "figures","PCA_NN_only.png"), 
      width= 15, height = 15, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = "forestgreen",
               repel = TRUE))
  dev.off()
  
  
  png(filename = here::here("outputs", "figures","PCA_countries_NN_only.png"), 
      width= 30, height = 20, units = "cm", res = 1000)
  print( factoextra::fviz_pca_biplot (pca,
                   select.var = list(cos2 = 0.1),
                   geom="point", pointshape=21,
                   col.ind = "black",
                   fill.ind = NCP_site$SiteCountry,
                   palette = scico::scico(38, palette = "roma"),
                   addEllipses = F, label = "var",
                   col.var = "black", repel = TRUE,
                   legend.title = "Country",
                   font.legend = c(6, "plain", "black")
  ))
  dev.off()
  
  #### NN PCA in dimensions 3 and 4 *(for variables with cos2 \> 0.1)*
  
  factoextra::fviz_pca_var(pca, col.var = "cos2",
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
    geom_histogram( bins = 60, aes(fill=after_stat(x)), color="grey50", size=0.3)+
    scale_fill_gradient2(low = "black", mid = "white", high = "forestgreen", guide = "colourbar")+
    labs(x = "Nature for nature PC1 site's evaluation")+
    geom_vline(xintercept = 0, linetype = 2, col= "black")+
    theme_minimal()
  histo_NN
  ggsave(filename = here::here("outputs", "figures","hist_NN_in_PC1.png"), histo_NN, width = 8, height =6 )
  
  
  
  
  ###------- NS: Nature contributions to Society -------
  NS <- names(grp)[ grp=="NS" ]
  NCP_NS <- NCP_site[,NS]
  rownames(NCP_NS) <- rownames(NCP_site)
  
  pca <- FactoMineR::PCA(NCP_NS, scale.unit = T, graph=F, ncp=10)
  
  png(filename = here::here("outputs", "figures","barplot_variance_explained_by_NS.png"), 
      width= 15, height = 10, units = "cm", res = 1000)
  print( factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 35), 
                       barfill = "dodgerblue3", barcolor = "dodgerblue3"))
  dev.off()
  
  
  png(filename = here::here("outputs", "figures","contribution_NS_in_dimensions.png"), 
      width= 17, height = 15, units = "cm", res = 1000)
  var <- factoextra::get_pca_var(pca)
  print( corrplot::corrplot(var$contrib, is.corr=FALSE) )
  dev.off()
  
  #### Nature for Society variations in the 2 first dimensions 
  png(filename = here::here("outputs", "figures","PCA_NS_only.png"), 
      width= 15, height = 15, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = "dodgerblue3",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE))
  dev.off()
  
  #### NS PCA in dimensions 3 and 4 *(for variables with cos2 \> 0.1)*
  png(filename = here::here("outputs", "figures","PCA_NS_only_dim3&4.png"), 
      width= 15, height = 15, units = "cm", res = 1000)
  print( factoextra::fviz_pca_var(pca, col.var = "cos2",
               axe=c(3,4),
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE,
               select.var = list(cos2 = 0.1)))
  dev.off()

  
  #### Surveys positions in PC1
  NS_PC1_surveys <- pca$ind$coord[,1]
  summary(NS_PC1_surveys)

  NN_NS_PC1_surveys <- data.frame(SiteCode = NCP_site$SiteCode ,
                                  NN_PC1=NN_PC1_surveys,
                                  NS_PC1=NS_PC1_surveys)
  rownames(NN_NS_PC1_surveys) <- rownames(NCP_site)
  
  save(NN_NS_PC1_surveys, file = here::here("outputs", "NN_NS_score_PC1.Rdata"))
  
  ##NS
  histo_NS <- ggplot(NN_NS_PC1_surveys, aes(x = NS_PC1)) +
    geom_histogram( bins = 60, aes(fill=..x..), color="grey30", size=0.3)+
    scale_fill_gradient2(low = "black", mid = "white", high = "dodgerblue3", guide = "colourbar")+
    labs(x = "Nature to People PC1 site's evaluation")+
    geom_vline(xintercept = 0, linetype = 2, col= "black")+
    theme_minimal()
  histo_NS
  ggsave(filename = here::here("outputs", "figures","hist_NS_in_PC1.png"), histo_NS, width = 8, height =6 )
  

  ## ------- NN against NS -------
    #Define colors by quarter
    #up right
  quarter_up_right <- NN_NS_PC1_surveys |>
    dplyr::filter(NN_PC1 > 0 & NS_PC1 > 0) |>
    dplyr::mutate(product_u_r =  pnorm(NN_PC1 * NS_PC1) )
  
  pal <- grDevices::colorRampPalette(c("white", "darkred")) (nrow(quarter_up_right))
  
  quarter_up_right <- quarter_up_right |>
    dplyr::arrange(product_u_r) |>
    dplyr::mutate( cols = pal)
  
  #up left
  quarter_up_left <- NN_NS_PC1_surveys |>
    dplyr::filter(NN_PC1 < 0 & NS_PC1 > 0) |>
    dplyr::mutate(product_u_l =  pnorm(-NN_PC1 * NS_PC1) )
  
  pal <- grDevices::colorRampPalette(c("white", "dodgerblue3")) (nrow(quarter_up_left))
  
  quarter_up_left <- quarter_up_left |>
    dplyr::arrange(product_u_l) |>
    dplyr::mutate( cols = pal)
  
  #down right
  quarter_down_right <- NN_NS_PC1_surveys |>
    dplyr::filter(NN_PC1 > 0 & NS_PC1 < 0) |>
    dplyr::mutate(product_d_r =  pnorm(NN_PC1 * -NS_PC1))
  
  pal <- grDevices::colorRampPalette(c("white", "forestgreen")) (nrow(quarter_down_right))
  
  quarter_down_right <- quarter_down_right |>
    dplyr::arrange(product_d_r) |>
    dplyr::mutate( cols = pal)
  
  #down left
  quarter_down_left <- NN_NS_PC1_surveys |>
    dplyr::filter(NN_PC1 < 0 & NS_PC1 < 0) |>
    dplyr::mutate(product_d_l =  pnorm(NN_PC1 * NS_PC1))
  
  pal <- grDevices::colorRampPalette(c("white", "grey20")) (nrow(quarter_down_left))
  
  quarter_down_left <- quarter_down_left |>
    dplyr::arrange(product_d_l) |>
    dplyr::mutate( cols = pal)
  
  
  #Merge the 4 quarters
  NN_NS_with_col <- quarter_down_left |>
    dplyr::full_join(quarter_down_right) |>
    dplyr::full_join(quarter_up_left) |>
    dplyr:: full_join(quarter_up_right)
  
  
  
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
               data = head(NN_NS_with_col[order(NN_NS_with_col$product_u_l, decreasing = T),], 15)) +
    
    
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
  ggsave( here::here("outputs", "figures", "Sites in NN and NS PC1.jpg"), plot = NN_NS_plot, width=10, height = 8 )
  

  
  
  ##------- Plot on map -------
  ### NS according to NN
  
  #data
  NN_NS_with_col <- NN_NS_with_col |>  dplyr::left_join(NCP_site[,c("SiteCode", "SiteLongitude", "SiteLatitude")])
  #plot
  function_NN_NS_on_map <- function(coord_NN_NS = NN_NS_with_col, ylim = c(-36, 31),
                                    xlim= c(-180,180)){
    ggplot(coord_NN_NS) +
    geom_sf(data = coast, color = NA, fill = "lightgrey") +
    
    #down left quarter
    geom_point(data= dplyr::filter(coord_NN_NS, is.na(product_d_l)==F),
               size = 2, alpha = 0.5,
               aes(x = SiteLongitude, y = SiteLatitude,
                   colour= product_d_l, alpha = product_d_l)) +
    scale_colour_gradient("product_d_l",
                          low = "white", high="grey20",
                          na.value=NA) +
    ggnewscale::new_scale("colour") +
    
    #up left quarter
    geom_point(data= dplyr::filter(coord_NN_NS, is.na(product_u_l)==F),
               size = 2, alpha = 0.5,
               aes(x = SiteLongitude, y = SiteLatitude,
                   colour=product_u_l, alpha = product_u_l)) +
    scale_colour_gradient("product_u_l",
                          low = "white", high="dodgerblue3",
                          na.value=NA) +
    geom_point(aes(x = SiteLongitude, y = SiteLatitude),
               shape = 1, size = 2, stroke = 1,
               color= "black",
               data = head(coord_NN_NS[order(coord_NN_NS$product_u_l, decreasing = T),], 15)) +
    
    ggnewscale::new_scale("colour") +
    
    #down right quarter
    geom_point(data= dplyr::filter(coord_NN_NS, is.na(product_d_r)==F),
               size = 2, alpha = 0.5,
               aes(x = SiteLongitude, y = SiteLatitude,
                   colour=product_d_r, alpha = product_d_r)) +
    scale_colour_gradient("product_d_r",
                          low = "white", high="forestgreen",
                          na.value=NA) +
    geom_point(aes(x = SiteLongitude, y = SiteLatitude),
               shape = 1, size = 2, stroke = 1,
               color= "black",
               data = head(coord_NN_NS[order(coord_NN_NS$product_d_r, decreasing = T),], 15)) +
    
    ggnewscale::new_scale("colour") +
    
    #up right quarter
    geom_point(data= dplyr::filter(coord_NN_NS, is.na(product_u_r)==F),
               size = 2,
               aes(x = SiteLongitude, y = SiteLatitude,
                   colour= product_u_r, alpha = product_u_r)) +
    scale_colour_gradient("product_u_r",
                          low = "white", high="darkred",
                          na.value=NA) +
    #add transparency
    scale_alpha_continuous(range = c(0, 1)) +
    
    
    coord_sf(xlim= xlim, ylim = ylim, expand = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    scale_size_continuous(range = c(0.5, 4), guide = FALSE) +
    theme_minimal()+
    labs(title = "Trade-offs in Nature Based Contribution",
         x="Longitude", y= "Latitude") +
    theme(legend.position = "none",
          plot.title = element_text(size=10, face="bold"),
          # axis.text.x = element_blank(),
          # axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm")
    )
  }
  
  world_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_col, 
                                        ylim = c(-36, 31),
                                        xlim= c(-180,180))
  ggsave( here::here("outputs", "figures", "world map with NN and NS PC1.jpg"), plot = world_map_NN_NS, width=10, height = 6 )
  
  
  australia_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_col, 
                                              ylim = c(-39, 0),
                                              xlim= c(100,180))
  ggsave( here::here("outputs", "figures", "australia map with NN and NS PC1.jpg"), plot = australia_map_NN_NS, width=10, height = 6 )
  
  
  polynesia_map_NN_NS <- function_NN_NS_on_map( NN_NS_with_col, 
                                                ylim = c(-25, -10),
                                                xlim= c(-180,-130))
  ggsave( here::here("outputs", "figures", "polynesia map with NN and NS PC1.jpg"), plot = polynesia_map_NN_NS, width=10, height = 6 )
  
  ### NN worldwide
  map_NN <- ggplot(NN_NS_with_col) +
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
    scale_size_continuous(range = c(0.5, 4), guide = "none") +
    theme_minimal()+
    labs(title = "Nature for Nature",
         x="Longitude", y= "Latitude") +
    theme(
      legend.position = "none",
      plot.title = element_text(size=10, face="bold"),
    )
  ggsave( here::here("outputs", "figures", "world map with NN PC1.jpg"), plot = map_NN, width=10, height = 6 )
  
  
  
  
  ### NS worldwide
  map_NS <- ggplot(NN_NS_with_col) +
    geom_sf(data = coast, color = NA, fill = "lightgrey") +
    geom_point(aes(x = SiteLongitude, y = SiteLatitude,
                   color = NS_PC1, alpha= abs(NS_PC1),
                   size = 2)) +
    geom_point(aes(x = SiteLongitude, y = SiteLatitude),
               shape = 1, size = 2, stroke = 1,
               color= "black",
               data = head(NN_NS_with_col[order(NN_NS_with_col$product_u_l, decreasing = T),], 15)) +
    
    scale_colour_gradientn(colours = colorRampPalette(rev(c( "dodgerblue3","white", "grey30")))(1000)) +
    scale_alpha_continuous(range = c(0, 1)) +
    
    coord_sf(ylim = c(-36, 31), expand = FALSE) +
    scale_size_continuous(range = c(0.5, 4), guide = "none") +
    theme_minimal()+
    labs(title = "Nature to People",
         x="Longitude", y= "Latitude") +
    theme(
      legend.position = "none",
      plot.title = element_text(size=10, face="bold"),
    )
  
  ggsave( here::here("outputs", "figures", "world map with NS PC1.jpg"), plot = map_NS, width=10, height = 6 )
  
  
} # End of function plot_ACP_NCP

#run plot_ACP_NCP function
plot_PCA_NCP(NCP_site)
