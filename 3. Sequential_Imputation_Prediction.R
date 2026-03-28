############################################################
## MODEL EVALUATION SCRIPT
## Combines imputed and no-imputation model fits
############################################################

library(tidyverse)
library(furrr)
library(future)
library(R2jags)
library(MLmetrics)
library(here)

############################################################
## JAGS MODEL
############################################################
mars_model <- function() {
  x[1] ~ dnorm(0,1)
  B ~ dnorm(0,1)
  R ~ dnorm(0,1)
  tauQ ~ dgamma(0.01,0.01)
  tauR ~ dgamma(0.01,0.01)
  
  mu_x[1] ~ dnorm(0,1)
  for(i in 2:n_years) {
    mu_x[i] <- B * x[i-1] + R * rain[i-1]
    x[i] ~ dnorm(mu_x[i], tauQ)
    y[i] ~ dnorm(x[i], tauR)
  }
}

############################################################
## ERROR CALCULATION FUNCTION
############################################################
get_error_for_model_seed <- function(model, species_seed_, rain_seed_, 
                                     input_data, n_imputations = NA) {
  
  model_data <- input_data |>
    filter(species_seed == species_seed_, rain_seed == rain_seed_) |>
    select(species_sim, !!sym(model), rain_sim)
  
  # Split train/test
  train_indices <- seq(1, round(nrow(model_data) * 0.8))
  train_data <- model_data[train_indices, ]
  test_data <- model_data[-train_indices, ]
  
  # JAGS parameters
  jags.data <- list(
    y = train_data[[model]],
    n_years = nrow(train_data),
    rain = train_data$rain_sim
  )
  
  jags.params <- c("x","B","tauQ","tauR","R")
  
  jags_model <- jags(
    jags.data,
    inits = NULL,
    parameters.to.save = jags.params,
    model.file = mars_model,
    n.chains = 3,
    n.thin = 10,
    n.burnin = 6000,
    n.iter = 10000,
    DIC = FALSE,
    jags.seed = 919
  )
  
  means <- jags_model$BUGSoutput$mean
  
  # Preidict on test data
  x <- numeric(nrow(test_data))
  x[1] <- train_data$species_sim[nrow(train_data)]
  set.seed(123)
  
  for(t in 2:nrow(test_data)) {
    mean_x <- means$B * x[t-1] + means$R * test_data$rain_sim[t-1]
    sd_x <- 1 / sqrt(means$tauQ)
    x[t] <- rnorm(1, mean = mean_x, sd = sd_x)
  }
  
  tibble(
    MAE = MLmetrics::MAE(x, test_data$species_sim),
    nrows = nrow(model_data),
    n_imputations = n_imputations
  )
}

############################################################
## LOAD SIMULATION-IMPUTED DATA
############################################################
sim_seq_imp_results <- read_csv(here("R/Manuscript - Missing Value Imputation/results/Sequential_Imputation_results.csv"))

miss_ids <- c(8, 9, 10, 12, 19, 24, 38, 39, 43, 55, 56, 58, 65, 68, 75, 80, 82, 83, 84, 85, 89, 90)

############################################################
## MODEL FIT FOR IMPUTED DATA
############################################################
seed_model_combs <- expand_grid(
  sim_seq_imp_results |> distinct(species_seed, rain_seed),
  model = expand_grid("sp", as.character(seq(1,14)), c("kf","mean","ma","mice")) %>%
    as.matrix() %>%
    apply(1, paste0, collapse = "_")
)

plan(multisession, workers = 4)

# Imputed models
model_fit_imputed <- split(seed_model_combs, seq_len(nrow(seed_model_combs))) |>
  future_map_dfr(function(row) {
    # Count number of imputed points from step
    step_num <- as.numeric(strsplit(row$model, "_")[[1]][2])
    bind_cols(row,
              get_error_for_model_seed(row$model, row$species_seed, row$rain_seed,
                                       input_data = sim_seq_imp_results,
                                       n_imputations = step_num))
  })

############################################################
## MODEL FIT - NO IMPUTATION (MISSING RECORDS ARE IGNORED)
############################################################
seed_model_combs_noimp <- expand_grid(
  sim_seq_imp_results |> distinct(species_seed, rain_seed),
  model = paste0("sp_", seq(1,14), "_priori")
)

noimp_test <- imap_dfc(miss_ids[1:14], function(i, k) {
  vec <- sim_seq_imp_results$species_sim
  vec[miss_ids[1:k]] <- NA
  vec
}) |>
  setNames(paste0("sp_", seq(1,14), "_priori")) |>
  bind_cols(sim_seq_imp_results |> select(species_sim, rain_seed, species_seed, rain_sim))

model_fit_noimp <- split(seed_model_combs_noimp, seq_len(nrow(seed_model_combs_noimp))) |>
  future_map_dfr(function(row) {
    # Count number of missing points
    step_num <- as.numeric(strsplit(row$model, "_")[[1]][2])
    bind_cols(row,
              get_error_for_model_seed(row$model, row$species_seed, row$rain_seed,
                                       input_data = noimp_test,
                                       n_imputations = step_num))
  })

############################################################
## COMBINE RESULTS
############################################################
model_fit_final <- bind_rows(
  model_fit_imputed |> mutate(imputation_type = "imputed"),
  model_fit_noimp  |> mutate(imputation_type = "no_imputation")
)

saveRDS(model_fit_final,
        here("results/model_fit_final.rds"))
