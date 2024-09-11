# CORAL REEF CONTRIBUTIONS

Research compendium to reproduce analyses and figures of the article: 
_Contributions of the world’s tropical reef fish communities to nature and people_ 
by Flandrin _et al._ published in One Earth.


## General

This repository is structured as follow:

- `data_raw/`: contains raw data from RLS protocol, don't delete this file
- `data/`: contains cleaned data from 'data_raw' required to reproduce figures and tables
- `R/`: contains R functions developed for this project
- `analyses/`: contains files organized by theme. Each files contains R scripts
to run specific analysis
- `outputs/`: contains intermediate results, final table and all figures 



## Workflow

### Contributions assessement

Most of the contributions used in this work have been assessed independently by 
researcher of the 'REEF FUTURE' project

#### Recycling
The script 'recycling/analyses/recycling_analyses.R' 
  - uses the work of **Sebastien Villéger** and **Nina Schiettekatte** (see 
  **Schiettekatte et al. 2022** for the methods), to assess from metabolic 
  parameters and 'rfishflux' package, the recycling of N and P by fishes.
  
  - produces the table 'recycling/outputs/flux_final_data_surveys.Rdata' that 
  contains the amount of nitrogen and phosphorous recycled by fishes, aggregated 
  at the survey level
  
> :boom: WARNING: The recycling data are under embargo, and will be published soon.
The final table aggregating excretion data at the survey scale is available to
enable the rest of the analysis. Please contact **Sebastien Villéger** for further
informations.

  

#### Biodiversity
The script 'biodiversity/analyses/biodiversity_analysis.R' 
  - produces the occurrence matrices of species, and calculates functional 
  diversity indices at the survey scale.
  
  - produces phylogenetic diversity indices of each surveys, from occurrence
  matrices and the phylogeny of fishes from **Siqueira et al. 2020**
  
  - assesses the richness of surveys in elasmobranch and endangered species
  according to the IUCN.
  
  
#### Productivity
The script 'productivity/analyses/productivity_analysis.R' 
  - uses the work of **Raphael Seguin** and **Nicolas Loiseau** (**Seguin et al. 2023**) 
  to assess the biomass production and biomass turnover (=productivity) of all 
  observed fish in each survey.
  
  - produces the 'productivity/outputs/RLS_prod_transect.Rdata' file containing
  production and productivity measures at the survey scale.


#### Nutrient content 
The script 'nutrients/analyses/nutrients_analysis.R' 
  - uses the work of **Eva Maire et al. 2021** by taking the observed and inferred
  nutrient content of tropical fishes in the file 
  'nutrients/outputs/nutrient_sp_data.RData'
  
  - aggregates all nutrients at the surveys scale, and produces the table of the 
  amount of nutrients available in 100g fresh portion of a mean fish 
  ('nutrients/outputs/nutrient_concentration_surveys.Rdata')
  

#### Carbonates excretion
The script 'carbonates/analyses/carbonates_analysis.R' 
  - uses the work of **Mattia Ghilardi** and **Sonia Bejarano**  
  (**Ghilardi et al. 2023**) by taking the observed and inferred
  carbonates excretions of tropical fishes in the file 
  'carbonates/outputs/survey_caco3_production.rds'
  
  - aggregates the excretions at the surveys scale in the table
  'carbonates/outputs/caco3_per_day.Rdata'
  
  
#### Trophic web
The script 'trophic_web/trophic_web_analysis.R' infers a trophic metaweb among
all the species observed in the RLS surveys, based on the work of 
**Camille Albouy** and **Dominique Gravel** (**Albouy et al. 2019**) using a 
 allometric niche model of trophic interactions.
  - the niche model is calibrated with observed data from **Barnes et al. 2008**
  - Local trophic web are extracted form the global metaweb using the occurrence 
  matrix.
  -Local trophic indicators are calculated and saved in the file 
  'trophic_indicators_surveys.Rdata'.


#### cultural
The script 'cultural_contributions/analyses/cultural_analyses.R'
  - uses the aesthetic values of each surveys, assessed by **Nicolas Mouquet** 
  and **Matthew McClean** in **McClean et al. 2023** (see **Langlois et al. 2022** 
  for the method of measuring species' aesthetic)
  
  - uses the cultural value of each species, assessed by **Nicolas Mouquet et al. 2023**
  and gathers the information at the survey scale.


### Contributions analysis

After assessing all the contributions, the script 'make.R' merges all of them and 
study their distribution, dimensionality and correlations. This script calls both 
scripts '2a_make_fig_1.R' and '2b_make_fig_2.R' to produce  both Figure 1 and 2 
of the paper Flandrin at al. 





## File with all results  

The file `outputs/tropical_reef_contributions_final_table.csv` contains all the information 
used and produced in this study at the reef level for the 1,237 sampled sites 
concerned. You are welcome to use it by citing properly our work, but even more 
welcome to contact us (ulysse.flandrin@gmail.com) if you want to collaborate :smiley:




## Figures and Tables

Figures will be stored in `outputs/figures/`.

The following Figures and Tables can be reproduced with the script indicated in 
brackets:
    
- Figure 1 (`R/2a_make_fig_1.R`)

- Figure 2 (`R/2b_make_fig_2.R`)

- Figure 3 has been produced by other mean.

- Figure S1 (`R/1b_Plot_contributions.R`)

- Figure S2 (`R/1b_Plot_contributions.R`)

- Figure S3 (`R/1c_PCA_analyses_on_contributions.R`)

- Figure S4 (`R/1c_PCA_analyses_on_contributions.R`)

- Figure S5 (`R/1b_Plot_contributions.R`)

- Figure S6 (`R/1b_Plot_contributions.R`)

- Figure S7 (`R/1d_weighted_mean_NP_NN_score.R`)

- Figure S8 (`R/1d_weighted_mean_NP_NN_score.R`)

- Figure S9 (`R/1e_spatial_autocorrelation.R`)

- Figure S10 (`R/1d_weighted_mean_NP_NN_score.R`)

- Figure S11 (`R/1d_weighted_mean_NP_NN_score.R`)

- Figure S12 (`R/1b_Plot_contributions.R`)

- Figure S13 (`R/1f_test_composite_scores_NP_NN.R`)




## Usage

Clone the repository and run this command in R/RStudio:

```r
source("make.R")
```

All required packages will be installed (if necessary) and loaded.

> :boom: WARNING: running `make.R` calls all the scripts and takes hours so if 
you want to work on one or a few scripts, you should run lines 17-45 of 
`make.R` and then go to the other script. Note that the lines 47-72 let to assess
each contributions in each surveys; you can run the lines 76-108 to quickly 
reproduce all the figures of the paper Flandrin et al. 

Enjoy!