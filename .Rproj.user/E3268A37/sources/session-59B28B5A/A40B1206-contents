# Note: 
# A c++ toolchain is needed
# See here for details on how to install: https://mc-stan.org/docs/2_26/cmdstan-guide/cmdstan-installation.html#installing-the-c-toolchain

# packages that are needed
packages <- c(
  "purrr",
  "tidybayes",
  "ggplot2",
  "fishflux",
  "Rphylopars",
  "dplyr",
  "tidyr",
  "drake",
  "e1071",
  "brms",
  "utils",
  "patchwork",
  "fishualize",
  "ggnewscale",
  "tidytree",
  "ggtree",
  "forcats",
  "ggrepel",
  "tidyr",
  "readr",
  "Rcpp",
  "broom",
  "flextable",
  "officer",
  "drake",
  "rnaturalearth",
  "rnaturalearthdata")

# install packages that are missing
to_install <- packages[!packages %in% installed.packages()]

install.packages(to_install)

# we recommend running this is a fresh R session or restarting your current session
install.packages("cmdstanr", 
                 repos = c("https://mc-stan.org/r-packages/", 
                           getOption("repos")))

# check toolchain
cmdstanr::check_cmdstan_toolchain()

# cmdstan installation
cmdstanr::install_cmdstan()
