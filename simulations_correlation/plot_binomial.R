source("simulations_correlation/utils_plot_metrics.R") 
library(BayesGlm)

base <- "simulations_correlation/irls2_normal_start_zero/binomial_filtrato/Nmcmc_1000_alpha_0.2_a0_0.05_y0_0.5_0.5_n500"

config <- list(
  name = "binomial",
  base_path = base,
  figure_dir = paste0(base, "/figures")
)

source("simulations_correlation/compare_results_seed.R")


 