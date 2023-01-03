################################################################################
##
## Analyse correlations among phylogenetic indices, ane between phylogenetic
## indices and taxonomic richness
##
## phylogenetic_indices_correlations.R
##
## 21/10/2022
##
## Ulysse Flandrin
##
################################################################################

## cleaning memory
rm(list=ls())

##-------------loading data-------------
load(here::here("biodiversity", "outputs", "phylogenetic_indices_surveys_without_outliers.Rdata"))
load(here::here("biodiversity", "outputs", "surveys_biodiversity.Rdata"))


# Filter
phylo_indices <-phylo_indices_surveys %>%
  dplyr::left_join(surveys_biodiversity_without_outliers)%>%
  dplyr::select(SurveyID, taxo_richness, ED_Mean, PD_Mean, SES_PD_Mean, residuals_PD_richness, PE_Mean, phylo_entropy )



##-------------correlation among indices-------------

### overview with ACP
phylo_indices <- questionr::na.rm(phylo_indices)
pca <- FactoMineR::PCA(phylo_indices[,-c(1,2)], scale.unit = T)
print(pca)
summary(pca)

factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))


png( here::here("biodiversity", "figures", "PCA_phylogenetic_indices_axes1-2.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Évite le chevauchement de texte 
) #give the quality of representatin of each variables in the 2 first dimensions
dev.off()

png( here::here("biodiversity", "figures", "PCA_phylogenetic_indices_axes3-4.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
fviz_pca_var(pca, col.var = "cos2",
             axe=c(4,4),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)
dev.off()



png( here::here("biodiversity", "figures", "Correlation ED ~ SES.PD.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
plot(phylo_indices$ED_Mean, phylo_indices$SES_PD_Mean)
dev.off()
png( here::here("biodiversity", "figures", "Correlation ED ~ PD.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
plot(phylo_indices$ED_Mean, phylo_indices$SES_PD_Mean)
dev.off()
png( here::here("biodiversity", "figures", "Correlation PE ~ SES.PD.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
plot(phylo_indices$PE_Mean, phylo_indices$SES_PD_Mean)
dev.off()


##-------------correlation with taxonomic richness-------------
png( here::here("biodiversity", "figures", "Correlation ED ~ nb species.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
plot(phylo_indices$ED_Mean ~ phylo_indices$taxo_richness)
dev.off()

png( here::here("biodiversity", "figures", "Correlation phylogenetic entropy ~ nb species.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
plot(phylo_indices$phylo_entropy ~ phylo_indices$taxo_richness)
dev.off()

png( here::here("biodiversity", "figures", "Correlation SES.PD ~ nb species.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
plot(phylo_indices$SES_PD_Mean ~ phylo_indices$taxo_richness)
dev.off()

png( here::here("biodiversity", "figures", "Correlation residuals PD ~ nb species.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
plot(phylo_indices$residuals_PD_richness ~ phylo_indices$taxo_richness)
dev.off()

png( here::here("biodiversity", "figures", "Correlation PE ~ nb species.png"),  width = 15, height = 15, units = "cm", pointsize = 12, res = 500)
plot(phylo_indices$PE_Mean ~ phylo_indices$taxo_richness)
dev.off()


