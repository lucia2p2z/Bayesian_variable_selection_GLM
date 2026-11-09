source("utils_simulations_baselines.R") 
library(BayesGlm)

base <- "simulations_with_baselines/irls2_start_glm/poisson_cd/Nmcmc_5000_a0_0.001_y0_0_1"


config <- list(
  name = "poisson_cd",
  base_path = base,
  figure_dir =  paste0(base, "/figures_with_baseline_new3")
)

source("simulations_with_baselines/run_plot_seed_cd.R")
#source("simulations_with_baselines/run_plot_poisson_specfic_cd.R")