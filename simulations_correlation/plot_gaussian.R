source("simulations_correlation/utils_plot_metrics.R") 
library(BayesGlm)

base <- "simulations_correlation/irls2_normal_start_zero/gaussian/Nmcmc_1000_alpha_1_a0_1e-04_y0_1_4_n100"

config <- list(
  name = "gaussian",
  base_path = base,
  figure_dir = paste0(base, "/figures")
)

source("simulations_correlation/compare_results_seed.R")


 