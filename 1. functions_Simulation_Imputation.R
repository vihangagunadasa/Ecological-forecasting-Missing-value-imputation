############################################################
## SIMULATION + IMPUTATION SCRIPT
## Modes:
##   "fixed"     -> fixed missing indices (STEP 1)
##   "sequential" -> progressive missingness (STEP 2)
############################################################

# Load required libraries
library(tidyverse)
library(mice)
library(forecast)
library(tseries)
library(imputeTS)
library(R2jags)
library(cli)
library(glue)
library(here)

############################################################
## HELPER FUNCTIONS
############################################################

# Filter dataset for a given species and site
# `included_years` is the years for which we have continuous species abundance and rain data
make_dataset <- function(species, site_name, abun_df, included_years) {
  
  df <- abun_df %>%
    select(Year, SiteName,
           abun_species = !!sym(species),
           rain = RainBom1) %>%
    filter(SiteName == site_name)
  
  # Add rows for missing years to ensure a continuous time series
  df <- rbind(
    df,
    tibble(
      Year = unique(abun_df$Year)[!unique(abun_df$Year) %in% unique(df$Year)],
      SiteName = site_name,
      abun_species = NA,
      rain = NA
    )
  ) %>% 
    arrange(Year) %>% 
    filter(Year %in% included_years)
  
  return(df)
}

############################################################
## SIMULATION MODELS 
#### STATE-SPACE MODEL FOR ABUNDANCE
#### ARIMA FOR RAIN
############################################################

## Simulation model for abundance ##
# Define Bayesian state-space model (latent process + observation model)
model_func <- function() {
  
  # Prior for autoregressive coefficient
  B ~ dnorm(0,0.01)
  
  # Process variance (state evolution)
  tauQ ~ dgamma(0.01, 0.01)
  
  # Observation variance (measurement error)
  tauR ~ dgamma(0.01, 0.01)
  
  # Initial latent state
  x[1] ~ dnorm(0, 0.1)
  
  # State process (latent dynamics)
  for(t in 2:n_years) {
    x[t] ~ dnorm(B * x[(t-1)], tauQ)
  }
  
  # Observation model (linking latent state to observed data)
  for(t in 1:n_years) {
    y[t] ~ dnorm(x[t],tauR)
  }
}

# Fit state-space model using JAGS
fit_ssm <- function(std_log_y) {
  
  # Prepare data for JAGS
  jags.data <- list("y" = std_log_y, "n_years" = length(std_log_y))
  
  # Run MCMC sampling
  jags(
    jags.data,
    inits = NULL,
    parameters.to.save = c("x","B","tauQ","tauR"),
    model.file = model_func,
    n.chains = 3,
    n.thin = 10,
    n.burnin = 6000,
    n.iter = 10000,
    DIC = FALSE,
    jags.seed = 919
  )
}
##

## Simulation model for abundance ##
# Simulate values from a fitted ARIMA model
# Extracts AR and MA coefficients and simulate series
rain_sim_function <- function(arima_model, n_steps) {
  
  # Extract AR and MA coefficients
  ar_coefs <- arima_model$coef[grep("ar",names(arima_model$coef))]
  ma_coefs <- arima_model$coef[grep("ma",names(arima_model$coef))]
  sd_val <- sqrt(arima_model$sigma2)
  
  # Determine burn-in length based on model order
  if(length(ar_coefs) == 0 && length(ma_coefs) == 0){
    # random walk arima(0,0,0)
    n_start <- 1
  } else {
    # burn-in 'n.start' must be as long as 'ar + ma'
    n_start <- length(ar_coefs) + length(ma_coefs)
  }
  
  # Generate simulated time series
  arima.sim(
    list(ar = ar_coefs, ma = ma_coefs),
    n = n_steps,
    n.start = n_start,
    start.innov = rep(0, times = n_start),
    sd = sd_val
  )
}
##

############################################################
## SIMULATION
############################################################

# Simulate species abundance from fitted state-space model
# Converts latent standardized values back to original scale
simulate_abundance <- function(means, std_log_y, mean_log_y, 
                               sd_log_y, sim_years, species_seed) {
  
  set.seed(species_seed)
  
  # Initialize latent state (x) and observation (y)
  x <- y <- vector(mode = "numeric", length = sim_years)
  x[1] <- std_log_y[1]
  
  # Simulate latent state process
  for(t in 2:sim_years) {
    mean_x <- means$B * x[t-1]
    sd_x <- (1 / sqrt(means$tauQ))
    x[t] <- rnorm(n = 1, mean = mean_x, sd = sd_x)
  }
  
  # Simulate observation process
  for(t in 1:sim_years) {
    mean_y <- 1 * x[t]
    sd_y <- (1/sqrt(means$tauR))
    y[t] <- rnorm(n = 1, mean = mean_y,sd = sd_y)
  }
  
  # Back-transform to original abundance scale
  abs(exp((y * sd_log_y) + mean_log_y) - 1)
}


# Simulate rainfall time series using ARIMA model
simulate_rain <- function(rain_vec, sim_years, rain_seed, save_path) {
  
  # If simulation not already saved, generate new one
  if(!file.exists(save_path)) {
    print("Generating fresh data")
    
    # Standardize rainfall
    mean_rain <- mean(rain_vec, na.rm = TRUE)
    sd_rain   <- sd(rain_vec)
    
    std_rain <- (rain_vec - mean_rain)/sd_rain
    
    # Fit ARIMA model
    model <- auto.arima(std_rain)
    
    # Simulate rainfall series
    set.seed(rain_seed)
    rain_sim <- rain_sim_function(model, sim_years)
    
    # Convert back to original scale
    abs_rain_sim <- abs((rain_sim * sd_rain) + mean_rain)
    
    # Save model and simulated data for reproducibility
    saveRDS(list(model = model, sim_abs = abs_rain_sim, sim = rain_sim), save_path)
    
  } else {
    print("Reading saved data")
    abs_rain_sim <- readRDS(save_path)$sim_abs
  }
  
  return(abs_rain_sim)
}

############################################################
## IMPUTATION
############################################################

# Apply different imputation methods to missing species data
run_imputation <- function(df, method) {
  
  if(method == "mice") {
    # Multiple Imputation by Chained Equations (predictive mean matching)
    set.seed(10)
    imp <- mice(df, m=5, maxit=50, method='pmm', seed=500, printFlag=FALSE)
    
    # Average imputed values across multiple datasets
    vals <- apply(as.matrix(imp$imp$species_sim),1,mean)
    df$species_sim[as.numeric(rownames(imp$imp$species_sim))] <- vals
    
  } else if(method == "mean") {
    # Mean imputation
    df$species_sim <- na_mean(df$species_sim)
    
  } else if(method == "ma") {
    # Moving average imputation
    df$species_sim <- na_ma(df$species_sim)
    
  } else if(method == "kf") {
    # Kalman filtering-based imputation
    df$species_sim <- na_kalman(df$species_sim, type="level")
  }
  
  return(df$species_sim)
}

############################################################
## MASTER FUNCTION
############################################################

# Main workflow: simulate data, introduce missingness, and apply imputation
run_simulation_imputation <- function(
    mode,
    species,
    SiteName,
    sim_years,
    species_seed,
    rain_seed,
    miss_ids,
    included_years
) {
  
  cli::cli_h1(glue("{species}, {SiteName}"))
  
  # Load raw dataset
  abun_df <- read.csv(
    here("R/Manuscript - Missing Value Imputation/data/qryMammalSpVegRainYEAR.csv")
  ) |> as_tibble()
  
  # For fixed mode, remove first row (assumed preprocessing step)
  if(mode == "fixed") abun_df <- abun_df[-1,]
  
  # Prepare dataset for selected species and site
  dataset <- make_dataset(
    species, SiteName, abun_df,
    included_years
  )
  
  # Extract species abundance
  target_y <- dataset$abun_species
  
  # Log-transform and standardize data
  log_y <- log(target_y + 1)
  mean_log_y <- mean(log_y)
  sd_log_y <- sd(log_y)
  std_log_y <- (log_y - mean_log_y)/sd_log_y
  
  # Fit state-space model
  jags_model <- fit_ssm(std_log_y)
  means <- jags_model$BUGSoutput$mean
  
  # Simulate species abundance
  abs_sim_species <- simulate_abundance(
    means, std_log_y, mean_log_y, sd_log_y,
    sim_years, species_seed
  )
  
  # Simulate rainfall data
  abs_rain_sim <- simulate_rain(
    dataset$rain,
    sim_years,
    rain_seed,
    here("R/Manuscript - Missing Value Imputation/results/fixed_imputation_rain",
         glue("{SiteName}_{sim_years}_{rain_seed}.rds"))
  )
  
  # Combine into working dataframe
  df <- tibble(species_sim = abs_sim_species, rain_sim = abs_rain_sim)
  
  ##########################################################
  ## MODE SWITCH
  ##########################################################
  
  if(mode == "fixed") {
    
    # Store true values before introducing missingness
    orig <- df$species_sim[miss_ids]
    df$species_sim[miss_ids] <- NA
    
    # Apply all imputation methods
    mice_imp  <- run_imputation(df, "mice")
    mean_imp  <- run_imputation(df, "mean")
    ma_imp    <- run_imputation(df, "ma")
    kf_imp    <- run_imputation(df, "kf")
    
    results <- tibble(
      year = miss_ids,
      original = orig,
      mice = mice_imp[miss_ids],
      mean = mean_imp[miss_ids],
      ma   = ma_imp[miss_ids],
      kf   = kf_imp[miss_ids]
    ) |> 
      mutate(species_name = species, site = SiteName, 
             sim_years = sim_years, species_seed = species_seed,
             rain_seed = rain_seed)
    
    write_csv(
      results,
      here("results/fixed_imputation_sp_abundance",
           glue("{species}_{SiteName}_{sim_years}_{species_seed}_{rain_seed}.csv"))
    )
    
  } else {
    
    # we need to have the simulated abundnace data iwth where missingness occurred
    df <- df |> mutate(species_miss_full = replace(species_sim, miss_ids, NA))
    
    # Sequentially increase missingness and evaluate each method
    results <- expand_grid(
      miss_index = seq_along(miss_ids),
      method = c("mice","mean","ma","kf")
    ) %>%
      pmap_dfc(function(miss_index, method) {
        
        # Introduce increasing missing values
        ids <- miss_ids[1:miss_index]
        df$species_sim[ids] <- NA
        
        # Apply imputation
        imputed <- run_imputation(df, method)
        
        tibble(!!glue("sp_{miss_index}_{method}") := imputed)
      })
    
    # Combine with working dataframe
    final <- bind_cols(df, results) %>%
      mutate(species_seed = species_seed,
             rain_seed = rain_seed)
    
    write_csv(
      final,
      here("results/sequential_imputation_sp_abundance_90",
           glue("{species}_{SiteName}_{sim_years}_{species_seed}_{rain_seed}.csv"))
    )
  }
}
