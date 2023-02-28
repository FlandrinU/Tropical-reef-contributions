################################################################################
##
## Aggregate cultural contributions estimated on species at the surveys scale
##
## public_and_scientific_interest_surveys.R
##
## 28/02/2023
##
## Ulysse Flandrin
##
################################################################################

rm(list=ls())

##-------------loading data------------
cultural_contrib <- read.csv(here::here("cultural_contributions", "data", 
                      "cultural_contributions.csv"), sep = ";", dec= ",")
load(here::here("biodiversity", "outputs", "occurrence_matrix_sp_survey_01.Rdata"))
load( here::here("data", "data_species.Rdata"))

##-------------resume data at species scale------------
cultural <- cultural_contrib |>
  dplyr::mutate(public_interest = (Wiki_views + Flickr)/2,
                scientific_interest = (NCBI + Sci_litt)/2 ,
                species_corrected = gsub(" ", "_", Scientific)) 

data_species <- tibble::column_to_rownames(data_species, var = "species")
colnames(surveys_sp_occ) <- data_species[colnames(surveys_sp_occ), "species_corrected"]


##-------------aggregates cultural scores at survey scale------------
colnames(surveys_sp_occ)[which(is.element(colnames(surveys_sp_occ),cultural$species_corrected) == F)]
# "Cirrhilabrus_laboutei" is lacking

cultural <- tibble::column_to_rownames(cultural, var ="species_corrected")
cultural_survey <- lapply( rownames(surveys_sp_occ), function(id){
  sp <- names(surveys_sp_occ[id, which(surveys_sp_occ[id,] >0)])
  if( length(sp) > 0){
    scientific_interest <- mean(cultural[sp, "scientific_interest"])
    public_interest <- mean(cultural[sp, "public_interest"])
    wiki_verna <- mean(cultural[sp, "Wiki_verna"])
    c(id, scientific_interest, public_interest, wiki_verna)
  }else{c(id, NA, NA, NA)}
})

cultural_contribution_surveys <- data.frame(do.call(rbind, cultural_survey)) |>
  dplyr::rename(SurveyID=X1, scientific_interest=X2, public_interest=X3,
                wiki_verna = X4) |>
  dplyr::mutate(scientific_interest = as.numeric(scientific_interest),
                public_interest = as.numeric(public_interest),
                wiki_verna = as.numeric(wiki_verna))

save(cultural_contribution_surveys, file = here::here("cultural_contributions",
                    "outputs", "cultural_contributions_surveys.Rdata"))

##-------------plot data------------
plot(cultural_contribution_surveys$scientific_interest ~ cultural_contribution_surveys$public_interest)
##study correlations
rownames(cultural_contribution_surveys) <- cultural_contribution_surveys$SurveyID
pca <- FactoMineR::PCA(questionr::na.rm(cultural_contribution_surveys[,-1]), 
                       scale.unit = T, graph=T, ncp=10)
factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 80))
var <- factoextra::get_pca_var(pca)
corrplot::corrplot(var$contrib, is.corr=FALSE)  
factoextra::fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE )
site_coord_in_pca <- as.data.frame(pca$ind$coord) 
plotly::plot_ly(site_coord_in_pca, 
        x= ~Dim.1, y= ~Dim.2, z= ~Dim.3,
        size = 10)
