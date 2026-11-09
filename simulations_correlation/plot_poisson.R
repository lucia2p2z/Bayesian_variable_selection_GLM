source("simulations_correlation/utils_plot_metrics.R") 
library(BayesGlm)

base <- "simulations_correlation/irls2_normal_start_zero/poisson/Nmcmc_1000_alpha_1_a0_0.001_y0_0_1_n100"

config <- list(
  name = "poisson",
  base_path = base,
  figure_dir = paste0(base, "/figures")
)

source("simulations_correlation/compare_results_seed.R")


 