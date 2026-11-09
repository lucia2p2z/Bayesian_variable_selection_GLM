library(BoomSpikeSlab)
library(bayesreg)
library(BayesGlm)
library(bestglm)
library(ggplot2)
source('experiments_pollution_data/utils_pollution_data.R')
source('experiments_pollution_data/run_experiments_pollution.R')
source('bayes_glmvs_oracle_class_s4.R')
source('utils_plot_experiments.R')

alpha <- 15
a0 <- 0.1

folder_figures <- "experiments_pollution_data/figures_data"
folder_results <- paste0("experiments_pollution_data/results_noHC_alpha",
                         alpha,"_a0_", a0, "_", "20260612") #format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(folder_results))  dir.create(folder_results)
if (!dir.exists(folder_figures)) dir.create(folder_figures)

data("mcdonald")
dataset <- mcdonald
df <- preprocess(dataset, folder_figures)

num_mcmc <- 10000
burn.in <- 1000
seed <- 1234
y0_mean <- 0
y0_std <- 1
set.seed(seed)
hyperparam = list(
  alpha = alpha,
  paramtau = list(shape=3, rate =1),
  y0_mean = y0_mean,
  y0_std = y0_std,
  y0 =  rnorm(nrow(df), y0_mean, y0_std), 
  a0 = a0,
  method_sir='irls2', 
  proposal='normal', 
  num_mcmc = num_mcmc + burn.in,
  burn.in = burn.in,
  beta0 = 'zero',
  seed=seed,
  folder_results = folder_results
)

hyperparam_baselines = list(
  priorspike = NULL,
  num_mcmc = num_mcmc,
  num_mcmc_hs = num_mcmc,
  burnin_hs=1e4,
  seed=seed,
  folder_results = folder_results
)

run_only_gibbs(df, hyperparam)
run_only_gibbs_baselines(df, hyperparam_baselines)


folder_figures <- paste0(folder_results, "/figures_model_paper")
if (!dir.exists(folder_figures)) dir.create(folder_figures)

mod <- readRDS(file.path(folder_results, "mod.Rds"))
mod_z1 <- readRDS(file.path(folder_results, "mod_z1.Rds"))
mod_glm <- readRDS(file.path(folder_results, "mod_glm.Rds"))
mod_spike <- readRDS(file.path(folder_results, "mod_spike.Rds"))
mod_hs <- readRDS(file.path(folder_results, "mod_hs.Rds"))
names_variables <- readRDS(file.path(folder_results, "names_variables.Rds"))

pdf(paste0(folder_figures,'/betas_times_z_with_baselines.pdf'), width = 10)
p <- plot_beta_z_boxplot_with_baselines(mod, mod_z1, mod_glm, 
                                        mod_spike, mod_hs, names_variables)
print(p)
dev.off()

pdf(paste0(folder_figures,'/betas_times_z.pdf'))
p <- plot_beta_z_boxplot_full(mod, mod_z1, mod_glm, names_variables)
print(p)
dev.off()