library(BoomSpikeSlab)
library(bayesreg)
library(BayesGlm)
library(bestglm)
library(ggplot2)
source('experiments_heart_data/utils_heart_dataset.R')
source('experiments_heart_data/run_experiments_heart.R')
source('bayes_glmvs_oracle_class_s4.R')
source('utils_plot_experiments.R')


a0 <- 0.01
alpha <- 1

folder_figures <- "experiments_heart_data/figures_data"
folder_results <- paste0("experiments_heart_data/results_tot_normal_optim_a0_", a0, 
                         "_alpha", alpha, "_", "20260618")  #format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(folder_results))  dir.create(folder_results)
if (!dir.exists(folder_figures)) dir.create(folder_figures)


dataset <- read.csv(file = "experiments_heart_data/dataset_heart_disease_uci.csv", sep = ",")
#summary(dataset)

df <- preprocess(dataset, folder_figures)


num_mcmc <- 2000
seed <- 1234
y0_mean <- 0.5
y0_std <- 0.5
set.seed(seed)
hyperparam = list(
  alpha = alpha,
  paramtau =1, 
  y0_mean = y0_mean,
  y0_std = y0_std,
  y0 = rbinom(nrow(df), 1, 0.5), 
  a0 = a0,
  method_sir='optim',
  proposal='normal', 
  num_mcmc = num_mcmc,
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

#run_only_gibbs(df, hyperparam)
#run_only_gibbs_baselines(df, hyperparam_baselines)


folder_figures <- paste0(folder_results, "/figures_model_paper")
if (!dir.exists(folder_figures)) dir.create(folder_figures)

mod <- readRDS(file.path(folder_results, "mod.Rds"))
mod_z1 <- readRDS(file.path(folder_results, "mod_z1.Rds"))
mod_glm <- readRDS(file.path(folder_results, "mod_glm.Rds"))
mod_spike <- readRDS(file.path(folder_results, "mod_spike.Rds"))
mod_hs <- readRDS(file.path(folder_results, "mod_hs.Rds"))
names_variables <- readRDS(file.path(folder_results, "names_variables.Rds"))

pdf(paste0(folder_figures,'/betas_times_z_with_baselines.pdf'), 
    width = 10, height = 7)
p <- plot_beta_z_boxplot_with_baselines(mod, mod_z1, mod_glm, 
                                        mod_spike, mod_hs, names_variables)
print(p)
dev.off()


pdf(paste0(folder_figures,'/betas_times_z.pdf'))
p <- plot_beta_z_boxplot_full(mod, mod_z1, mod_glm, names_variables)
print(p)
dev.off()

pdf(paste0(folder_figures,'/betas_only_baselines.pdf'))
p <- plot_beta_boxplot_only_baselines(mod_spike, mod_hs, names_variables)
print(p)
dev.off()
