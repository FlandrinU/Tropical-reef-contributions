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
cultural_contrib <- read.csv(here::here("cultural_contributions", "data", 
                      "cultural_contributions.csv"), sep = ";", dec= ",")
load(here::here("biodiversity", "outputs", "occurrence_matrix_sp_survey_01.Rdata"))
load( here::here("data", "data_species.Rdata"))

##-------------resume data at species scale------------
cultural <- cultural_contrib |>
  dplyr::mutate(public_interest = (Wiki_views + Flickr)/2,
                academic_knowledge = (NCBI + Sci_litt)/2 ,
                species_corrected = gsub(" ", "_", Scientific)) 

##-------------deal with names ------------
for( i in 1:nrow(cultural)){
  if(cultural$species_corrected[i] %in% data_species$species){
    cultural$species_corrected[i] <- data_species$species_corrected[
      which(data_species$species == cultural$species_corrected[i])]
  }
}

data_species <- tibble::column_to_rownames(data_species, var = "species")
colnames(surveys_sp_occ) <- data_species[colnames(surveys_sp_occ), "species_corrected"] # corrected names in surveys_sp_occ

#check uncommon species
colnames(surveys_sp_occ)[which(is.element(colnames(surveys_sp_occ),cultural$species_corrected) == F)] # Cirrhilabrus_laboutei is lacking


##-------------aggregates cultural scores at survey scale------------

cultural <- tibble::column_to_rownames(cultural, var ="species_corrected")
cultural_survey <- lapply( rownames(surveys_sp_occ), function(id){
  cat(id, "\n")
  sp <- names(surveys_sp_occ[id, which(surveys_sp_occ[id,] >0)])
  if( length(sp) > 0){
    mean_academic_knowledge <- mean(cultural[sp, "academic_knowledge"], na.rm = T)
    mean_public_interest <- mean(cultural[sp, "public_interest"], na.rm = T)
    mean_wiki_verna <- mean(cultural[sp, "Wiki_verna"], na.rm = T)
    
    academic_knowledge <- quantile(cultural[sp, "academic_knowledge"], 0.75, na.rm = T)
    public_interest <- quantile(cultural[sp, "public_interest"], 0.75, na.rm = T)
    wiki_verna <- quantile(cultural[sp, "Wiki_verna"], 0.75, na.rm = T)
    
    cult <- as.numeric(c(id, academic_knowledge, public_interest, wiki_verna, 
      mean_academic_knowledge, mean_public_interest, mean_wiki_verna))
  }else{cult <- c(id, NA, NA, NA,  NA, NA, NA)}
  names(cult) <- c("SurveyID", "academic_knowledge", "public_interest", "wiki_verna", 
                   "mean_academic_knowledge", "mean_public_interest", "mean_wiki_verna")
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
