library(BoomSpikeSlab)
library(bayesreg)
library(BayesGlm)
library(glm2)
library(ggplot2)

source('experiments_crabs_data/utils_crabs_dataset.R')
source('experiments_crabs_data/run_experiments_crabs.R')
source('bayes_glmvs_oracle_class_s4.R')
source('utils_plot_experiments.R')


folder_results <- paste0("experiments_crabs_data/results_", "20260612") #format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(folder_results))  dir.create(folder_results)

data(crabs)
dataset <- crabs

df <- preprocess(dataset, NULL)
seed <- 1234
num_mcmc <- 3000
hyperparam_baselines = list(
  priorspike ='zellner',
  num_mcmc = num_mcmc,
  num_mcmc_hs = num_mcmc,
  burnin_hs=1e4,
  seed=seed,
  folder_results = folder_results
)

y0_mean <- 0
y0_std <- 1
set.seed(seed)
hyperparam = list(
  alpha = 1,
  paramtau =1, 
  y0_mean = y0_mean,
  y0_std = y0_std,
  y0 =  rep(3, nrow(df)), 
  a0 = 0.01,
  method_sir='irls2', 
  proposal='normal', 
  num_mcmc = num_mcmc,
  beta0 = 'zero',
  seed=seed,
  folder_results = folder_results
)

run_only_gibbs_baselines(df, hyperparam_baselines)
run_only_gibbs(df, hyperparam)

folder_figures <- paste0(folder_results, "/figures_model_paper")
if (!dir.exists(folder_figures))  dir.create(folder_figures)
mod <- readRDS(file.path(folder_results, "mod.Rds"))
mod_z1 <- readRDS(file.path(folder_results, "mod_z1.Rds"))
mod_glm <- readRDS(file.path(folder_results, "mod_glm.Rds"))
mod_spike <- readRDS(file.path(folder_results, "mod_spike.Rds"))
mod_hs <- readRDS(file.path(folder_results, "mod_hs.Rds"))
names_variables <- readRDS(file.path(folder_results, "names_variables.Rds"))

pdf(paste0(folder_figures,'/betas_times_z_with_baselines.pdf'), width=10, height=8)
p <- plot_beta_z_boxplot_with_baselines(mod, mod_z1, mod_glm, 
                                        mod_spike, mod_hs, names_variables)
print(p)
dev.off()

pdf(paste0(folder_figures,'/betas_only_baselines.pdf'))
p <- plot_beta_boxplot_only_baselines(mod_spike, mod_hs, names_variables)
print(p)
dev.off()

pdf(paste0(folder_figures,'/betas_times_z.pdf'))
p <- plot_beta_z_boxplot_full(mod, mod_z1, mod_glm, names_variables)
print(p)
dev.off()


# posterior inclusion computed as:
#post.incl.prob1 <- colMeans(model1$beta != 0)


