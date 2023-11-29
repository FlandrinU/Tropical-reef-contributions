#################################################################################
#'
#'This script determines the total biomass of fishery species in each transect,
#' using the list of families which are fished in the world 
#' (expert opinion, Cinner et al. 2016)
#'
#'@author Ulysse Flandrin, \email{ulysse.flandrin@@gmail.com}
#'
#' @date 2022/11/24
#'
################################################################################

#-----------------Loading data---------------------
load(here::here("data", "data_surveys.Rdata"))
load(here::here("data", "data_species.Rdata"))
load(here::here("recycling", "outputs","flux_final_data_surveys.Rdata"))

#----------------- fishery species ---------------------
data_surveys_fishery <- data_surveys |>
  dplyr::left_join(data_species[,c("species", "family", "Size")]) |>
  dplyr::rename(max_size = Size)
table(data_surveys_fishery$family)

#families in which all species are targeted in the family: (Cinner et al. 2020 + other expert opinion )
all_sp <- c("Acanthuridae", "Caesionidae", "Carangidae", "Ephippidae", "Haemulidae", "Kyphosidae",
            "Labridae", "Lethrinidae", "Lutjanidae", "Mullidae", "Nemipteridae", "Scaridae",
            "Scombridae", "Serranidae", "Siganidae", "Sparidae", "Sphyraenidae",
            "Mugilidae", "Bothidae", "Scorpaenidae", "Sciaenidae")

#families with species targeted if larger than 20cm: (Cinner et al. 2020 + other expert opinion )
target_larger_20cm <- c("Balistidae", "Holocentridae", "Pomacanthidae", "Priacanthidae")

#non-targeted families: (Cinner et al. 2020 + other expert opinion )
untargeted_fam <- c("Chaetodontidae","Cirrhitidae", "Diodontidae", "Grammatidae", 
                   "Monacanthidae", "Pempheridae", "Pinguipedidae", "Pseudochromidae", 
                   "Synodontidae", "Tetraodontidae", "Zanclidae","Fistulariidae",
                   "Ostraciidae", "Pomacentridae")

unclassified_fam <- setdiff(unique(data_surveys_fishery$family), c(all_sp, target_larger_20cm, untargeted_fam))
unclassified_fam #OK


#keep only targeted fishes
data_surveys_fishery <- data_surveys_fishery |> 
  dplyr::filter(family %in% all_sp | (family %in% target_larger_20cm & max_size > 20))

# -------------------biomass of fishery species in each survey -------------------
surveys_fishery_sp_biom <- data_surveys_fishery |> 
  dplyr::select(SurveyID, species, biomass) |>
  dplyr::group_by(SurveyID, species) |>
  dplyr:: summarize( fishery_biomass=sum(biomass) ) |>
  tidyr::pivot_wider(names_from = species, values_from = fishery_biomass, values_fill = 0) |>
  tibble::column_to_rownames(var="SurveyID") 


surveys_fishery_biom <- as.data.frame(rowSums(surveys_fishery_sp_biom))
surveys_fishery_biom <- tibble::rownames_to_column(surveys_fishery_biom)
colnames(surveys_fishery_biom) <- c("SurveyID", "fishery_biomass")


# -------------------save data -------------------
save(surveys_fishery_biom, file = here::here("nutrients", "outputs", "fishery_tot_biomass.Rdata"))
