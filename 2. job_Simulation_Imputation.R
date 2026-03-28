############################################################
## JOB SCRIPT FOR SIMULATION + IMPUTATION
############################################################

library(future)
library(future.apply)
library(progressr)
library(tidyverse)
library(here)

# Source unified simulation-imputation functions
source(here("1. functions_Simulation_Imputation.R"))

############################################################
## PARAMETERS
############################################################

# Species and site info
species <- "Notomys.alexis"
SiteName <- "Tobermorey East"
included_years <- 1998:2011 |> as.character()

# Mode type
mode_type = "sequential" # or 'sequential'

# Simulation years and missingness
if(mode_type == "fixed"){
  sim_years_vec <- c(60, 90, 900, 9000)
  miss_ids <- c(8, 9, 10, 12, 19, 24, 39, 43, 48, 55)
  
} else {
  sim_years_vec <- c(90)
  miss_ids <- c(8, 9, 10, 12, 19, 24, 38, 39, 43, 55, 56, 58,
                65, 68, 75, 80, 82, 83, 84, 85, 89, 90)
}

# Seeds
set.seed(69)
total_seeds <- 50
species_seeds <- c(1996,2023,sample(9999, (total_seeds - 2), replace = F))
rain_seeds <- c(165,1568,sample(9999, (total_seeds - 2), replace = F))
seed_combs <- paste(species_seeds,rain_seeds, sep = "_")


############################################################
## CREATE COMBINATIONS
############################################################

# Expand grid for all combinations
combs_tbl <- expand_grid(
  sim_years = sim_years_vec,
  seed_combs = seed_combs,
  species = species,
  SiteName = SiteName
) |>
  tibble::tibble() |> 
  tidyr::separate("seed_combs",sep="_",into = c("species_seed", "rain_seed")) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.character))


############################################################
## RUN SIMULATIONS
############################################################

## Parallel setup
plan(future.callr::callr, workers = 6)

with_progress({
  
  p <- progressor(steps = nrow(combs_tbl))
  
  results <- future_lapply(seq_len(nrow(combs_tbl)), function(i) {
    
    row <- combs_tbl[i, ]
    
    p()
    
    run_simulation_imputation(
      mode = mode_type,
      species = row$species,
      SiteName = row$SiteName,
      sim_years = row$sim_years,
      species_seed = row$species_seed,
      rain_seed = row$rain_seed,
      miss_ids = miss_ids,
      included_years = included_years
    )
    
  }, future.seed = TRUE)
})

############################################################
## SAVE RESULTS
############################################################

# For fixed mode (CSV)
fixed_results <- results %>% 
  keep(~ inherits(.x, "tbl_df") && all(c("mice","mean","ma","kf") %in% names(.x)))

if(length(fixed_results) > 0) {
  write_csv(
    bind_rows(fixed_results),
    here("results",
         "Fixed_imputation_TS_lengths_results.csv")
  )
}


# For sequential mode (RDS)
sequential_results <- results %>% 
  keep(~ inherits(.x, "tbl_df") && "sp_1_mice" %in% names(.x))  # heuristic: sequential results have sp_* columns

if(length(sequential_results) > 0) {
  write_csv(
    bind_rows(sequential_results),
    here("results",
         "Sequential_Imputation_results.csv")
  )
}

