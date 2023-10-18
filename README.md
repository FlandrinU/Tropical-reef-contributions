# CORAL REEF CONTRIBUTIONS

Research compendium to reproduce analyses and figures of the article: 
_..._ 
by Flandrin _et al._ 2023 published in .


## General

This repository is structured as follow:

- `data_raw/`: contains raw data from RLS protocol, don't delete this file
- `data/`: contains cleaned data from 'data_raw' required to reproduce figures and tables
- `R/`: contains R functions developed for this project
- `analyses/`: contains folders organized by theme. Each folder contains R scripts
to run specific analysis
- `results/`: contains intermediate results and all figures



## File with all results  

The file `outputs/tropical_reef_contributions_final_table.csv` contains all the information 
used and produced in this study at the reef level for the 1,237 sampled sites 
concerned. You are welcome to use it by citing properly our work, but even more 
welcome to contact us (ulysse.flandrin@gmail.com) if you want to collaborate :smiley:



## Figures and Tables

Figures will be stored in `outputs/figures/`.

The following Figures and Tables can be reproduced with the script indicated in 
brackets:
    
- Figure 1b (`deep/06_prediction_performances.R`)

- Figure S1 A (`features/cluster.R`)



## Usage

Clone the repository and run this command in R/RStudio:

```r
source("make.R")
```


All required packages will be installed (if necessary) and loaded.

> :boom: WARNING: running `make.R` calls all the scripts and takes days so if 
you want to work on one or a few scripts, you should run lines 17-42 of 
`make.R` and then go to the other script.

Enjoy!