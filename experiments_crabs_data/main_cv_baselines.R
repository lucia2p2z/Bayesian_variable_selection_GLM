library(BoomSpikeSlab)
library(bayesreg)
library(ggplot2)
library(caret)
library(glm2)
library(parallel)
source('experiments_crabs_data/utils_crabs_dataset.R')
source('experiments_crabs_data/run_experiments_crabs.R')

NUMSIM <- 10
vect_seed <- 1234 + 1:NUMSIM


folder_figures <- "experiments_crabs_data/figures_data"
folder_results <- paste0("experiments_crabs_data/resultsGLM_cv_90_10_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(folder_results))  dir.create(folder_results)
if (!dir.exists(folder_figures)) dir.create(folder_figures)

data(crabs)
dataset <- crabs
df <- preprocess(dataset)
saveRDS(df, file = file.path(folder_results, "df.Rds"))

# run training and prediction for the k-th fold 
f <- function(s) {
  set.seed(s)
  sample <- sample(c(TRUE, FALSE), nrow(df), replace=TRUE, prob=c(0.9,0.1)) 
  train  <- df[sample, ]
  test   <- df[!sample, ]
  
  num_mcmc <- 3000
  seed <- 1234
  set.seed(seed)
  folder_results <- paste0(folder_results, "/seed_", s)
  hyperparam = list(
    seed_data = s ,
    priorspike ='zellner',
    num_mcmc = num_mcmc,
    num_mcmc_hs = num_mcmc,
    burnin_hs=1e4,
    seed=seed,
    folder_results = folder_results
  )
  run_train_predict_baselines(df, train, test, hyperparam)
}

#mclapply(X = vect_seed, FUN=f, mc.cores = NUMSIM)

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

