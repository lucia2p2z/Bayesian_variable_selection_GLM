
run_only_gibbs <- function(df, hyperparam) {
  
  num_mcmc <- hyperparam$num_mcmc
  folder_results <- hyperparam$folder_results
  if (!dir.exists(folder_results))  dir.create(folder_results)
  
  folder_figures <- paste0(folder_results, "/figures_model")
  if (!dir.exists(folder_figures)) dir.create(folder_figures)
  
  # Standardize numerical features and perform one-hot encoding ==================
  num_vars <- c('Width') #, 'Rep1', 'Rep2')
  df[num_vars] <- scale(df[num_vars])

  x <- df[, !(names(df)%in%c("target"))]
  y <- df[, (names(df)%in%c("target"))]
  x_df <- prepare_dataset_for_BayesGlm(x)
  x_mt <- as.matrix(x_df)
  
  saveRDS(colnames(x_df), file.path(folder_results, "names_variables.Rds"))
  
  
  ### GLM model ================================================================== 
  data_train <- cbind(y, x_df)
  mod_glm <- glm(y ~.-1, family = poisson, data=data_train)

  saveRDS(mod_glm, file = file.path(folder_results, "mod_glm.Rds"))
  cat('Saved `mod_glm` into ', folder_results, "\n")
  
  
  ### BayesVS-GLM  ================================================================== 
  if (hyperparam$beta0 == 'glm_coefficients'){
    beta0 <- mod_glm$coefficients
  } else if (hyperparam$beta0 == 'zero'){
    beta0 <- rep(0, ncol(x))
  } else {
    stop('Unknow starting point of the chain. Please check: `hyperparam$beta0`. ')
  }
  
  #set.seed(seed)
  mod <- bayes_glmvs(x_mt, y, family=poisson.bayes(), control = hyperparam)
  mod <- fit(mod, num_mcmc, beta0=beta0, burn.in =floor(num_mcmc/10), asmcmc =F)
  saveRDS(mod, file = file.path(folder_results, "mod.Rds"))
  saveRDS(hyperparam, file = file.path(folder_results, "hyperparam.Rds"))
  cat('Saved `mod` and `hyperparam` into ', folder_results, "\n")
  
  ### Bayes-GLM  ===============================================================
  mod_z1 <- bayes_glmvs_oracle(x_mt, y, family=poisson.bayes(), 
                               control = hyperparam, z_true=rep(1, ncol(x)))
  mod_z1 <- fit(mod_z1, num_mcmc, beta0=beta0, burn.in =floor(num_mcmc/10), asmcmc =FALSE)
  saveRDS(mod_z1, file = file.path(folder_results, "mod_z1.Rds"))
  cat('Saved `mod_z1`  into ', folder_results, "\n")
  
  
  ### PLOTS:  ================================================================== 
  
  plot(mod, data_x=x_mt, figure_folder = folder_figures)
  plot(mod, data_x=x_mt, asmcmc=TRUE, figure_folder= folder_figures)
  
  # pdf(paste0(folder_figures,'/betas_times_z_bayesvsGLM.pdf'))
  # p <- plot_beta_z_boxplot(mod, mod_glm, x_df)
  # print(p)
  # dev.off()
  
  pdf(paste0(folder_figures,'/betas_times_z.pdf'))
  p <- plot_beta_z_boxplot_full(mod, mod_z1, mod_glm, colnames(x_df))
  print(p)
  dev.off()
    
  pdf(paste0(folder_figures,'/betas.pdf'))
  boxplot(mod@fit$beta, main = 'beta')
  dev.off()
  
}


run_only_gibbs_baselines <- function(df, hyperparam){
  
  num_mcmc <- hyperparam$num_mcmc
  num_mcmc_hs <- hyperparam$num_mcmc_hs
  burnin_hs <- hyperparam$burnin_hs
  
  folder_results <- hyperparam$folder_results
  if (!dir.exists(folder_results))  dir.create(folder_results)
  
  folder_figures <- paste0(folder_results, "/figures_model")
  if (!dir.exists(folder_figures)) dir.create(folder_figures)
  
  
  # Standardize numerical features and perform one-hot encoding ==================
  num_vars <- c('Width') 
  df[num_vars] <- scale(df[num_vars])
  
  x <- df[, !(names(df)%in%c("target"))]
  y <- df[, (names(df)%in%c("target"))]
  x_df <- prepare_dataset_for_BayesGlm(x)
  x_mt <- as.matrix(x_df)
  data_train <- cbind(y, x_df)
  
  saveRDS(colnames(x_df), file.path(folder_results, "names_variables.Rds"))
  
  ### SPIKE AND SLAB model =====================================================
  if (hyperparam$priorspike == 'zellner') {
    prior <- PoissonZellnerPrior(predictors = x_mt, counts = y)
  } else {
    prior <- NULL
  }
  mod_spike <- poisson.spike(y~.-1, data = data_train, niter=num_mcmc, prior = prior)
  saveRDS(mod_spike, file = file.path(folder_results, "mod_spike.Rds"))
  cat('Saved `mod_spike` into ', folder_results, "\n")
  
  
  ### Shrinkage HORSESHOE model ================================================
  mod_hs <- bayesreg(y~., data=data_train, model ="poisson", prior = "horseshoe", 
                     n.samples=num_mcmc_hs, burnin=burnin_hs)
  saveRDS(mod_hs, file = file.path(folder_results, "mod_hs.Rds"))
  cat('Saved `mod_hs` into ', folder_results, "\n")
  
  
  pdf(paste0(folder_figures,'/betas_only_baselines.pdf'))
  p <- plot_beta_boxplot_only_baselines(mod_spike, mod_hs, colnames(x_df))
  print(p)
  dev.off()
  
}


run_train_predict <- function(df, train, test, hyperparam) {
  
  num_mcmc <- hyperparam$num_mcmc
  folder_results <- hyperparam$folder_results
  if (!dir.exists(folder_results))  dir.create(folder_results)
  
  folder_figures <- paste0(folder_results, "/figures_model")
  if (!dir.exists(folder_figures)) dir.create(folder_figures)
  
  # Standardize numerical features and perform one-hot encoding ==================
  num_vars <- c('Width') #, 'Rep1', 'Rep2')
  train[num_vars] <- scale(train[num_vars])
  test[num_vars] <- scale(test[num_vars])
  
  xtrain <- train[, !(names(df)%in%c("target"))]
  ytrain <- train[, (names(df)%in%c("target"))]
  xtrain_df <- prepare_dataset_for_BayesGlm(xtrain)
  xtrain_mt <- as.matrix(xtrain_df)
  
  xtest <- test[, !(names(df)%in%c("target"))]
  ytest <- test[, (names(df)%in%c("target"))]
  xtest_df <- prepare_dataset_for_BayesGlm(xtest)
  xtest_mt <- as.matrix(xtest_df)
  
  saveRDS(colnames(xtrain_df), file.path(folder_results, "names_variables.Rds"))
  
  
  ### GLM model ================================================================
  data_train <- cbind(ytrain, xtrain_df)
  mod_glm <- glm(ytrain ~.-1, family = poisson, data=data_train)
  
  summary(mod_glm)
  
  
  saveRDS(mod_glm, file = file.path(folder_results, "mod_glm.Rds"))
  cat('Saved `mod_glm` into ', folder_results, "\n")
  
  
  ### BayesVS-GLM  =============================================================
  if (hyperparam$beta0 == 'glm_coefficients'){
    beta0 <- mod_glm$coefficients
  } else if (hyperparam$beta0 == 'zero'){
    beta0 <- rep(0, ncol(xtrain))
  } else {
    stop('Unknow starting point of the chain. Please check: `hyperparam$betA0`. ')
  }
  
  mod <- bayes_glmvs(xtrain_mt, ytrain, family=poisson.bayes(), 
                     control = hyperparam)
  mod <- fit(mod, num_mcmc, beta0=beta0, burn.in =floor(num_mcmc/10), asmcmc =F)
  saveRDS(mod, file = file.path(folder_results, "mod.Rds"))
  saveRDS(hyperparam, file = file.path(folder_results, "hyperparam.Rds"))
  cat('Saved `mod` and `hyperparam` into ', folder_results, "\n")
  
  ### Bayes-GLM  ===============================================================
  mod_z1 <- bayes_glmvs_oracle(xtrain_mt, ytrain, family=poisson.bayes(), 
                               control = hyperparam, z_true=rep(1, ncol(xtrain_mt)))
  mod_z1 <- fit(mod_z1, num_mcmc, beta0=beta0, burn.in =floor(num_mcmc/10), asmcmc =FALSE)
  saveRDS(mod_z1, file = file.path(folder_results, "mod_z1.Rds"))
  cat('Saved `mod_z1`  into ', folder_results, "\n")
  
  # PREDICTION with GLM and BayesGlm  ============================================
  pred_prob <- predict(mod, newdata=xtest_mt, type = "response")
  ypred <- rpois(length(pred_prob), lambda = pred_prob)
  
  pred_prob_z1 <- predict(mod_z1, newdata=xtest_mt, type = "response")
  ypred_z1 <- rpois(length(pred_prob_z1), lambda = pred_prob_z1)
  
  pred_prob_glm <- predict.glm(mod_glm, newdata = xtest_df, type = "response")
  ypred_glm <- rpois(length(pred_prob_glm), lambda = pred_prob_glm)
  
  mae_glm <- mean(abs(ytest - ypred_glm))
  rmse_glm <- sqrt(mean((ytest - ypred_glm)^2))
  ll_glm <- sum(dpois(ytest, lambda = pred_prob_glm, log = TRUE))
  
  mae <- mean(abs(ytest - ypred))
  rmse <- sqrt(mean((ytest - ypred)^2))
  ll <- sum(dpois(ytest, lambda = pred_prob, log = TRUE))
  
  mae_z1 <- mean(abs(ytest - ypred_z1))
  rmse_z1 <- sqrt(mean((ytest - ypred_z1)^2))
  ll_z1 <- sum(dpois(ytest, lambda = pred_prob_z1, log = TRUE))
  
  
  # Combine into a metrics matrix
  metrics <- rbind(
    MAE  = c(BayesVSGLM = mae, BayesGLM = mae_z1,  MLE = mae_glm),
    RMSE = c(BayesVSGLM = rmse, BayesGLM = rmse_z1, MLE = rmse_glm),
    LL = c(BayesVSGLM = ll, BayesGLM = ll_z1, MLE = ll_glm)
  )
  colnames(metrics) <- c("BayesVS-GLM", "Bayes-GLM", "MLE")
  saveRDS(metrics, file = file.path(folder_results, "prediction_metrics.Rds"))
  
  
  ### PLOTS:  ================================================================== 
  
  plot(mod, data_x=xtrain_mt, figure_folder= folder_figures)
  plot(mod, data_x=xtrain_mt, asmcmc=TRUE, figure_folder= folder_figures)
  
  pdf(paste0(folder_figures,'/betas_times_z.pdf'))
  p <- plot_beta_z_boxplot_full(mod, mod_z1, mod_glm, colnames(xtrain_df))
  print(p)
  dev.off()
  
  pdf(paste0(folder_figures,'/betas.pdf'))
  boxplot(mod@fit$beta, main = 'beta')
  dev.off()
  
}


run_train_predict_baselines <- function(df, train, test, hyperparam) {
  
  num_mcmc <- hyperparam$num_mcmc
  num_mcmc_hs <- hyperparam$num_mcmc_hs
  burnin_hs <- hyperparam$burnin_hs
  folder_results <- hyperparam$folder_results
  if (!dir.exists(folder_results))  dir.create(folder_results)
  
  folder_figures <- paste0(folder_results, "/figures_model")
  if (!dir.exists(folder_figures)) dir.create(folder_figures)
  
  # Standardize numerical features and perform one-hot encoding ==================
  num_vars <- c('Width') #, 'Rep1', 'Rep2')
  train[num_vars] <- scale(train[num_vars])
  test[num_vars] <- scale(test[num_vars])
  
  xtrain <- train[, !(names(df)%in%c("target"))]
  ytrain <- train[, (names(df)%in%c("target"))]
  xtrain_df <- prepare_dataset_for_BayesGlm(xtrain)
  xtrain_mt <- as.matrix(xtrain_df)
  
  xtest <- test[, !(names(df)%in%c("target"))]
  ytest <- test[, (names(df)%in%c("target"))]
  xtest_df <- prepare_dataset_for_BayesGlm(xtest)
  xtest_mt <- as.matrix(xtest_df)
  
  data_train <- cbind(ytrain, xtrain_df)
  
  saveRDS(colnames(xtrain_df), file.path(folder_results, "names_variables.Rds"))
  
  
  ### SPIKE AND SLAB model =====================================================
  if (hyperparam$priorspike == 'zellner') {
    prior <- PoissonZellnerPrior(predictors = xtrain_mt, counts = ytrain)
  } else {
    prior <- NULL
  }
  mod_spike <- poisson.spike(ytrain~.-1, data =data_train, niter=num_mcmc, prior = prior)
  saveRDS(mod_spike, file = file.path(folder_results, "mod_spike.Rds"))
  cat('Saved `mod_spike` into ', folder_results, "\n")
  
  
  ### Shrinkage HORSESHOE model ================================================
  mod_hs <- bayesreg(ytrain~., data=data_train, model ="poisson", prior = "horseshoe", 
                     n.samples=num_mcmc_hs, burnin=burnin_hs)
  saveRDS(mod_hs, file = file.path(folder_results, "mod_hs.Rds"))
  cat('Saved `mod_hs` into ', folder_results, "\n")
  
  ### Shrinkage HORSESHOE+ model ================================================
  mod_hsplus <- bayesreg(ytrain~., data=data_train, model ="poisson", prior = "horseshoe+", 
                       n.samples=num_mcmc_hs, burnin=burnin_hs)
  saveRDS(mod_hsplus, file = file.path(folder_results, "mod_hsplus.Rds"))
  cat('Saved `mod_hs` into ', folder_results, "\n")
  
  
  # PREDICTION   ============================================
  pred_prob_spike <- BoomSpikeSlab::predict.poisson.spike(mod_spike, newdata=xtest_mt, type = "response")
  pred_prob_spike <- apply(pred_prob_spike, 1, mean)
  ypred_spike <- rpois(length(pred_prob_spike), lambda = pred_prob_spike)
  
  pred_prob_hs <- bayesreg:::predict.bayesreg(mod_hs, newdata=xtest_df, type = "response")
  ypred_hs <- rpois(length(pred_prob_hs), lambda = pred_prob_hs)
  
  pred_prob_hsplus <- bayesreg:::predict.bayesreg(mod_hsplus, newdata=xtest_df, type = "response")
  ypred_hsplus <- rpois(length(pred_prob_hsplus), lambda = pred_prob_hsplus)
  
  
  mae_spike <- mean(abs(ytest - ypred_spike))
  rmse_spike <- sqrt(mean((ytest - ypred_spike)^2))
  ll_spike <- sum(dpois(ytest, lambda = pred_prob_spike, log = TRUE))
  
  mae_hs <- mean(abs(ytest - ypred_hs))
  rmse_hs <- sqrt(mean((ytest - ypred_hs)^2))
  ll_hs <- sum(dpois(ytest, lambda = pred_prob_hs, log = TRUE))
  
  mae_hsplus <- mean(abs(ytest - ypred_hsplus))
  rmse_hsplus <- sqrt(mean((ytest - ypred_hsplus)^2))
  ll_hsplus <- sum(dpois(ytest, lambda = pred_prob_hsplus, log = TRUE))
  
  
  # Combine into a metrics matrix
  metrics <- rbind(
    MAE  = c(SpikeSlab = mae_spike, Horseshoe = mae_hs,  HorseshoePlus = mae_hsplus),
    RMSE = c(SpikeSlab = rmse_spike, Horseshoe = rmse_hs, HorseshoePlus = rmse_hsplus),
    LL = c(SpikeSlab = ll_spike, Horseshoe = ll_hs, HorseshoePlus = ll_hsplus)
  )
  colnames(metrics) <- c("SpikeSlab", "Horseshoe", "Horseshoe+")
  saveRDS(metrics, file = file.path(folder_results, "prediction_metrics_baselines.Rds"))
  
  
  ### PLOTS:  ================================================================== 
  
  pdf(paste0(folder_figures,'/beta_spike.pdf'))
  boxplot(mod_spike$beta, main='Posterior Coefficients (SpikeSlab)')
  dev.off()
  
  pdf(paste0(folder_figures,'/beta_horseshoe.pdf'))
  boxplot(t(mod_hs$beta), main='Posterior Coefficients (Horseshoe)')
  dev.off()
  
  pdf(paste0(folder_figures,'/beta_horseshoeplus.pdf'))
  boxplot(t(mod_hsplus$beta), main='Posterior Coefficients (Horseshoe+)')
  dev.off()
  
}
