# Load libraries
library(tidyverse)
library(here)
library(patchwork)

# ===============================
# Load data
# ===============================
original_df           <- read.csv(here("data/qryMammalSpVegRainYEAR.csv"))
sim_seq_imp_results   <- read_csv(here("results/Sequential_Imputation_results.csv"))
sim_fixed_imp_results <- read_csv(here("results/Fixed_imputation_TS_lengths_results.csv"))
model_fit_final       <- readRDS(here("results/model_fit_final.rds"))

miss_ids <- c(8, 9, 10, 12, 19, 24, 38, 39, 43, 55, 56, 58, 65, 68, 75, 80, 82, 83, 84, 85, 89, 90)

# ===============================
# Helper: Centralized save function
# ===============================
save_plot <- function(plot, filename, width_mm = 500, height_mm = 500) { # size may vary from published in the manuscript
  ggsave(filename = filename, path = here("plots/"),
         plot = plot, width = width_mm, height = height_mm, units = "mm", limitsize = FALSE)
}

# ===============================
# Figure 1: Site map
# ===============================
# Site map was created through ArcGIS application.


# ===============================
# Figure 2: Time series of Species Abundance and Rain covariate
# ===============================
species   <- "Notomys.alexis"
SiteName  <- "Tobermorey East"
covariate <- "RainBom1"
title_cov <- "Rain"
years <- 1990:2012
  
species_site_df <- original_df %>%
  filter(SiteName == !!SiteName) %>%
  select(Year, SiteName, species = !!sym(species), rain = !!sym(covariate))
  
miss_years_species <- years[!years %in% species_site_df$Year]

cov_site_df <- species_site_df %>% tidyr::drop_na(rain)
miss_years_cov <- years[!years %in% cov_site_df$Year]
  
plot_species <- ggplot(species_site_df, aes(x = Year, y = species)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = miss_years_species, color = "red", linetype = "longdash") +
    labs(y = "Capture per 100 trap nights") +
    theme_light() +
    theme(strip.text.x = element_blank(),
        axis.text.x = element_blank(), # Remove x-axis text for alignment
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18),  
        axis.title.y = element_text(size = 14))

plot_rain <- ggplot(cov_site_df, aes(x = Year, y = rain)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = miss_years_cov, color = "red", linetype = "longdash") +
    labs(y = "Previous year Rainfall (mm)") +
    theme_light() +
    theme(strip.text.x = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),  
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 14))

final_plot <- (plot_species + labs(tag = "(A)")) /
  (plot_rain + labs(tag = "(B)")) +
  plot_layout(ncol = 1, heights = c(1, 1.2)) & 
  theme(plot.tag = element_text(size = 14, face = "bold"))  # Adjust label size


save_plot(final_plot, "ts_plot_species_rain.png")



# ===============================
# Figure 3: Sequential Arrangement of Missing Values
# ===============================
coords <- expand_grid(x = seq(8, 90, 1), y = seq(8, 90, 1))

coords_color <- map_dfr(seq_along(miss_ids), function(n) {
  tibble(x = miss_ids[n], y = miss_ids[1:n])
}) |> mutate(color = TRUE)

fig3 <- coords |>
  left_join(coords_color, by = c("x", "y")) |>
  replace_na(list(color = FALSE)) |>
  ggplot(aes(x = factor(x), y = factor(y), fill = color)) +
  geom_tile(color = "grey") +
  scale_fill_manual(values = c("white", "#349DCA")) +
  scale_x_discrete(breaks = c(1, miss_ids)) +
  scale_y_discrete(breaks = c(1, miss_ids)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 20, angle = 90, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20, vjust = -0.15),
    axis.title.y = element_text(size = 20, vjust = 2)
  ) +
  xlab("Position in the time series of the cumulative missing values") +
  ylab("Position in the time series used to cumulate missing values")

save_plot(fig3, "sequential_arrangement.png")



# ===============================
# Figure 4: MAE Across Simulation Lengths for Imputation Methods
# ===============================
df_long <- sim_fixed_imp_results |>
  pivot_longer(
    cols = c("mice", "mean", "ma", "kf"),
    names_to = "method",
    values_to = "imps"
  ) |>
  mutate(error = abs(original - imps))

method_levels <- c("Mean", "MA", "KF", "MICE")

line_plot_w_CI <- df_long |>
  group_by(site, sim_years, method) |>
  summarise(
    error_mean = mean(error),
    error_sd   = sd(error),
    n          = n(),
    .groups = "drop"
  ) |>
  mutate(
    method   = factor(method, labels = method_levels),
    ci_lower = error_mean - 1.96 * (error_sd / sqrt(n)),
    ci_upper = error_mean + 1.96 * (error_sd / sqrt(n))
  ) |>
  ggplot(aes(x = sim_years, y = error_mean)) +
  geom_point(aes(color = method), size = 3) +
  geom_line(aes(group = method, color = method), linewidth = 1) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = method),
              alpha = 0.2, color = NA) +
  scale_x_log10() +
  facet_wrap(~ site) +
  geom_vline(xintercept = c(60, 90, 900, 9000),
             color = "grey", linetype = "dashed") +
  xlab("Length of simulation") +
  ylab("Mean Absolute Error") +
  scale_color_discrete(name = "Method", breaks = method_levels) +
  scale_fill_discrete(name = "Method", breaks = method_levels) +
  theme_light() +
  theme(
    legend.position = "bottom",
    strip.text.x = element_blank(),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 1.5)),
    fill  = guide_legend(override.aes = list(alpha = 0.3))
  )

save_plot(line_plot_w_CI, "optimal_length.png")



# ===============================
# Figure 5: MAE by Missing Value Position for Sequential Imputation
# ===============================
violin_data <- map_dfr(seq_along(miss_ids), function(index) {
  sim_seq_imp_results |>
    select(starts_with(paste0("sp_", index, "_")), species_sim, species_miss_full, species_seed, rain_seed) |>
    pivot_longer(cols = starts_with(paste0("sp_", index, "_"))) |>
    separate(name, into = c("name", "n_missing", "method"), sep = "_") |>
    select(-name) |>
    group_by(rain_seed, species_seed, method) |>
    mutate(id = row_number()) |>
    filter(is.na(species_miss_full)) |>
    group_modify(function(data, grouping) {
      head(data, as.numeric(data$n_missing[1])) |>
        summarise(mae = mean(abs(species_sim - value)))
    }) |>
    mutate(n_missing = index)
}) |>
  ungroup() |>
  mutate(method = factor(method, levels = c("mean", "ma", "kf", "mice"),
                         labels = c("Mean", "MA", "KF", "MICE")))

method_colors <- c("Mean" = "#00BFC4", "MA" = "#7CAE00", "KF" = "#F8766D", "MICE" = "#C77CFF")

fig5 <- violin_data |>
  ggplot(aes(x = factor(n_missing), y = mae, fill = method)) +
  geom_boxplot() +
  scale_x_discrete(labels = factor(miss_ids)) +
  xlab("Position in the time series of the cumulative missing values") +
  ylab("Mean Absolute Error") +
  scale_fill_manual(name = "Method", values = method_colors) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15)
  )

save_plot(fig5, "imputation_error_comparison.png")



# ===============================
# Figure 6: Species Abundance Prediction through Sequential Imputation
# ===============================
bplot_full_data <- model_fit_final |>
  mutate(
    step_method = str_replace(model, "sp_", ""),
    step = as.numeric(str_extract(step_method, "^[0-9]+")),
    method = str_replace(step_method, "^[0-9]+_", ""),
    method = case_when(
      method == "mice" ~ "MICE",
      method == "mean" ~ "Mean",
      method == "ma"   ~ "MA",
      method == "kf"   ~ "KF",
      method == "priori" ~ "A Priori",
      TRUE ~ method
    )
  )

method_order <- c("A Priori", "Mean", "MA", "KF", "MICE")
method_colors <- c("A Priori" = "#D3D3D3", "Mean" = "#00BFC4", "MA" = "#7CAE00", "KF" = "#F8766D", "MICE" = "#C77CFF")

fig6 <- bplot_full_data |>
  ggplot(aes(x = factor(step), y = MAE, fill = factor(method, levels = method_order))) +
  geom_boxplot() +
  scale_x_discrete(labels = miss_ids) +
  xlab("Position in the time series of the cumulative missing values") +
  ylab("Mean Absolute Error") +
  scale_fill_manual(name = "Method", values = method_colors, labels = method_order) +
  theme_light() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 18),
    axis.title.x = element_text(size = 20, vjust = -0.75),
    axis.title.y = element_text(size = 20, vjust = 1.5),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15)
  )

save_plot(fig6, "imputation_prediction.png")
