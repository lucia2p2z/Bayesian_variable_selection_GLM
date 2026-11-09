source("utils_simulations_baselines.R") 
library(BayesGlm)

base <- "simulations_with_baselines/irls2_start_glm/gaussian_cd/Nmcmc_5000_a0_1e-04_y0_1_4"

config <- list(
  name = "gaussian_cd",
  base_path = base,
  figure_dir =  paste0(base, "/figures_with_baseline_new2")
)

source("simulations_with_baselines/run_plot_seed_cd.R")