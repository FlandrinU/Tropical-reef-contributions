README file for James Cook University's Tropical Data Hub repository, created on April 8, 2020

-------------------------------------------------
Citation
-------------------------------------------------

Trophic innovations fuel reef fish diversification - Nature Communications
Authors: Siqueira, Alexandre C; Morais, Renato A; Bellwood, David R; Cowman, Peter F
Corresponding author: Siqueira AC (alexandre.siqueira@my.jcu.edu.au)


-------------------------------------------------
Description of files
-------------------------------------------------
(1) Data folder

This folder contains the curated dataset ("Data_final_Siqueira_etal.csv") with classified geographical and ecological variables per reef fish species, as described in the manuscript.
Additionally, this dataset contains the median rate estimates (lambda, mu and net_div) derived from the BAMM analysis.

-------------------------------------------------
(2) BAMM folder

Within this folder, all the results from the BAMM analysis are stored in different sub-folders. Each sub-folder contains the outputs of a BAMM run in one (out of 100) TACT tree.
Also, this folder contains the R script ("BAMM_script_Siqueira_etal.R") used to extract the results from all BAMM runs and calculate the median rate estimates (provided in the final dataset).

-------------------------------------------------
(3) BRT folder

This folder refers to the Gradient Boosted Regression Tree analysis. 
It contains the R script ("Xgboost_script_Siqueira_etal.R") to reproduce the analysis and plot the respective results.

-------------------------------------------------
(4) MuSSE folder

This folder contains the R script ("MuSSE_script_Siqueira_etal.R") to reproduce the multistate speciation and extinction (MuSSE) model and plot the respective results.
Please, be aware that this analysis takes quite a long time to run. 

-------------------------------------------------
(5) Simmap folder

This folder contains the R script ("Simmap_script_Siqueira_etal.R") to perform the stochastic character mappings (simmaps) and plot the respective results.
The customized R function ("make.simmap.musse.R") to perform simmaps using MuSSE reconstructions (see methods) is also provided. It is internally sourced within the simmap R script.
Please, be aware that this analysis takes quite a long time to run.

-------------------------------------------------
(6) TACT folder

This folder contains the outputs of 100 runs of the stochastic polytomy resolution algorithm (Taxonomic Addition for Complete Trees [TACT]).
It also contains the file "Reef_fish_all_combined.trees" that consists of all trees combined in a single file (used as an input in the MuSSE and Simmap R scripts).

-------------------------------------------------

** Notes: 
* All R scripts require setting the working directory in which all these folders are saved
* All R package versions are provided within the respective scripts
* The "pre-cooked" results to generate the figures are provided along with the article in the file "Source Data"
* The scripts were written in R version 3.5.3 using a Mac OS 10.13.6 

Please contact corresponding author for any queries
