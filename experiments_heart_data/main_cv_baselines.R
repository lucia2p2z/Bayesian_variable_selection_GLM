library(BoomSpikeSlab)
library(bayesreg)
library(BayesGlm)
library(ggplot2)
library(caret)
library(parallel)
source('experiments_heart_data/utils_heart_dataset.R')
source('experiments_heart_data/run_experiments_heart.R')
source('bayes_glmvs_oracle_class_s4.R')
source('utils_plot_experiments.R')

NUMSIM <- 10
vect_seed <- 1234 + 1:NUMSIM

a0 <- 0.01
alpha <- 1

folder_figures <- "experiments_heart_data/figures_data"
folder_results <- paste0("experiments_heart_data/results_cv_normal_optim_a0_", a0, "_alpha", alpha, "_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(folder_results))  dir.create(folder_results)
if (!dir.exists(folder_figures)) dir.create(folder_figures)

dataset <- read.csv(file = "experiments_heart_data/dataset_heart_disease_uci.csv", sep = ",")
df <- preprocess(dataset, folder_figures)

f <- function(s) {
  set.seed(s)
  sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.90,0.10)) 
  train  <- df[sample, ]
  test   <- df[!sample, ]
  
  num_mcmc <- 2000
  seed <- 1234
  set.seed(seed)
  folder_results <- paste0(folder_results, "/seed_", s)
  hyperparam = list(
    seed_data = s ,
    priorspike = NULL,
    num_mcmc = num_mcmc,
    num_mcmc_hs = num_mcmc,
    burnin_hs=1e4,
    seed=seed,
    folder_results = folder_results
  )
  run_train_predict_baselines(df, train, test, hyperparam)
}

mclapply(X = vect_seed, FUN=f, mc.cores = NUMSIM)

