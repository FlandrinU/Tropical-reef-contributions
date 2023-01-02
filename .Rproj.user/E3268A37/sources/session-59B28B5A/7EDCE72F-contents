source("R/packages.R")  # loads packages
source("R/functions_wrangling.R") # defines the create_plot() function
source("R/functions_analysis.R") # defines the create_plot() function
source("R/functions_plots.R") # defines the create_plot() function
source("R/functions_tables.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan


future::plan(future.callr::callr)

config <- drake_config(plan,
                       lock_envir = FALSE,
                       garbage_collection = TRUE,
                       memory_strategy = "autoclean",
                       history = FALSE)
