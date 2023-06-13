################################################################################
##
## Aggregate cultural contributions estimated on species at the surveys scale
##
## public_and_academic_knowledge_surveys.R
##
## 28/02/2023
##
## Ulysse Flandrin, data from Nicolas Mouquet
##
################################################################################

rm(list=ls())

##-------------loading data------------
# cultural_contrib <- read.csv(here::here("cultural_contributions", "data", 
#                       "cultural_contributions.csv"), sep = ";", dec= ",")
cultural_contrib <- read.csv(here::here("cultural_contributions", "data", 
                                        "05_Human_Interest_final_table.csv"), 
                             sep = ";", dec= ",")
load(here::here("biodiversity", "outputs", "occurrence_matrix_sp_survey_01.Rdata"))
load( here::here("data", "data_species.Rdata"))

##-------------resume data at species scale------------
cultural <- cultural_contrib |>
  dplyr::rename(public_interest = public,
                academic_knowledge = acad,
                species_corrected = fb_sci_name) |> 
  dplyr::right_join(
    dplyr::select(data_species, species, species_corrected, spec_code))

## check list of species
uncommon_species <- cultural$species_corrected[ which(is.na(cultural$rls_sci_name))] #15 species are lacking in interest estimation
for( i in uncommon_species){
  cat(i, ": " )
  cat(table(surveys_sp_occ[,i]), "\n")
}   # quite rare species
##

##-------------aggregates cultural scores at survey scale------------
cultural <- tibble::column_to_rownames(cultural, var ="species")
cultural_survey <- lapply( rownames(surveys_sp_occ), function(id){
  cat(id, "\n")
  sp <- names(surveys_sp_occ[id, which(surveys_sp_occ[id,] >0)])
  if( length(sp) > 0){
    mean_academic_knowledge <- mean(cultural[sp, "academic_knowledge"], na.rm = T)
    mean_public_interest <- mean(cultural[sp, "public_interest"], na.rm = T)

    academic_knowledge <- quantile(cultural[sp, "academic_knowledge"], 0.75, na.rm = T)
    public_interest <- quantile(cultural[sp, "public_interest"], 0.75, na.rm = T)

    cult <- as.numeric(c(id, academic_knowledge, public_interest,
                         mean_academic_knowledge, mean_public_interest))
  }else{cult <- c(id, NA, NA, NA, NA)}
  
  names(cult) <- c("SurveyID", "academic_knowledge", "public_interest", 
                   "mean_academic_knowledge", "mean_public_interest")
  cult
})

cultural_contribution_surveys <- data.frame(do.call(rbind, cultural_survey)) |>
  dplyr::mutate( across(everything(), ~as.numeric(.))) |>
  dplyr::mutate(SurveyID = as.character(SurveyID))

save(cultural_contribution_surveys, file = here::here("cultural_contributions",
                    "outputs", "cultural_contributions_surveys.Rdata"))

##-------------plot data------------
plot(cultural_contribution_surveys$academic_knowledge ~ cultural_contribution_surveys$public_interest)
plot(cultural_contribution_surveys$mean_public_interest ~ cultural_contribution_surveys$public_interest)

##study correlations
rownames(cultural_contribution_surveys) <- cultural_contribution_surveys$SurveyID
pca <- FactoMineR::PCA(questionr::na.rm(cultural_contribution_surveys[,-1]), 
                       scale.unit = T, graph=T, ncp=10)
factoextra::fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100))
var <- factoextra::get_pca_var(pca)
corrplot::corrplot(var$contrib, is.corr=FALSE)  
factoextra::fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE )
site_coord_in_pca <- as.data.frame(pca$ind$coord) 
plotly::plot_ly(site_coord_in_pca, 
        x= ~Dim.1, y= ~Dim.2, z= ~Dim.3,
        size = 10)
