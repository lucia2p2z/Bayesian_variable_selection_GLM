source("utils_simulations_baselines.R") 
library(BayesGlm)

base <- "simulations_with_baselines/optim_start_glm/binomial_cd/Nmcmc_5000_a0_0.01_y0_0_1_smallbeta"

config <- list(
  name = "binomial_cd",
  base_path = base,
  figure_dir =  paste0(base, "/figures_baselines_new3")
)

source("simulations_with_baselines/run_plot_seed_cd.R")
