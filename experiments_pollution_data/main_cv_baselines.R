library(BoomSpikeSlab)
library(bayesreg)
library(BayesGlm)
library(ggplot2)
library(parallel)
library(bestglm)
source('experiments_pollution_data/utils_pollution_data.R')
source('experiments_pollution_data/run_experiments_pollution.R')
source('bayes_glmvs_oracle_class_s4.R')
source('utils_plot_experiments.R')

NUMSIM <- 30
vect_seed <- 1234 + 1:NUMSIM

alpha <- 15
a0 <- 0.1

folder_figures <- "experiments_pollution_data/figures_data"
folder_results <- paste0("experiments_pollution_data/results_noHC_cv5k_alpha",
                         alpha,"_a0_", a0, "_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(folder_results))  dir.create(folder_results)
if (!dir.exists(folder_figures)) dir.create(folder_figures)

data("mcdonald")
dataset <- mcdonald
df <- preprocess(dataset, folder_figures)

f <- function(s) {
  set.seed(s)
  sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.90,0.10)) 
  train  <- df[sample, ]
  test   <- df[!sample, ]
  
  num_mcmc <- 5000
  seed <- 1234
  set.seed(seed)
  folder_results <- paste0(folder_results, "/seed_", s)
  hyperparam = list(
    seed_data = s ,
    priorspike =NULL,
    num_mcmc = num_mcmc,
    num_mcmc_hs = num_mcmc,
    burnin_hs=1e4,
    seed=seed,
    folder_results = folder_results
  )
  run_train_predict_baselines(df, train, test, hyperparam)
}


#f(1235)
# f(1236)
# f(1237)
# f(1238)
# f(1239)
# f(1240)
# f(1241)
#f(1242)
f(1243)
f(1244)
f(1245)
f(1246)
f(1247)
f(1248)
f(1249)
f(1250)
f(1251)
f(1252)
f(1253)
f(1254)
f(1255)
f(1256)
f(1257)
f(1258)
f(1259)
f(1260)
f(1261)
f(1262)
f(1263)
f(1264)
#mclapply(X = vect_seed, FUN=f, mc.cores = NUMSIM)


