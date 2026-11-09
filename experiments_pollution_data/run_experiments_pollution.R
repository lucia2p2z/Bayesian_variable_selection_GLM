run_only_gibbs <- function(df, hyperparam) {
  
  num_mcmc <- hyperparam$num_mcmc
  folder_results <- hyperparam$folder_results
  if (!dir.exists(folder_results))  dir.create(folder_results)
  
  folder_figures <- paste0(folder_results, "/figures_model")
  if (!dir.exists(folder_figures)) dir.create(folder_figures)
  
  x <- df[, !(names(df)%in%c("target"))]
  y <- df[, (names(df)%in%c("target"))]
  x_df <- x # prepare_dataset_for_BayesGlm(x)
  x_mt <- as.matrix(x_df)
  
  saveRDS(colnames(x_df), file.path(folder_results, "names_variables.Rds"))
  
  burn.in <- ifelse(is.null(hyperparam$burn.in), floor(num_mcmc/10), hyperparam$burn.in) 
  
  
  ### GLM model ================================================================== 
  data_train <- cbind(y, x_df)
  mod_glm <- glm(y ~.-1, family = gaussian, data=data_train)

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
  mod <- bayes_glmvs(x_mt, y, family=gaussian.bayes(), control = hyperparam)
  mod <- fit(mod, num_mcmc, beta0=beta0, burn.in =burn.in, asmcmc =F)
  saveRDS(mod, file = file.path(folder_results, "mod.Rds"))
  saveRDS(hyperparam, file = file.path(folder_results, "hyperparam.Rds"))
  cat('Saved `mod` and `hyperparam` into ', folder_results, "\n")
  
  Z <- mod@fit$z
  ZZ <- apply(Z, 1, function(x) paste(x, collapse = ""))
  tZ <- table(ZZ)
  stZ <- sort(tZ, decreasing = TRUE)
  saveRDS(stZ, file=file.path(folder_results, "table_z.Rds"))
  
  ### Bayes-GLM  ===============================================================
  mod_z1 <- bayes_glmvs_oracle(x_mt, y, family=gaussian.bayes(), 
                               control = hyperparam, z_true=rep(1, ncol(x)))
  mod_z1 <- fit(mod_z1, num_mcmc, beta0=beta0, burn.in =burn.in, asmcmc =FALSE)
  saveRDS(mod_z1, file = file.path(folder_results, "mod_z1.Rds"))
  cat('Saved `mod_z1`  into ', folder_results, "\n")
  
  
  ### PLOTS:  ================================================================== 
  
  plot(mod, data_x=x_mt, figure_folder= folder_figures)
  plot(mod, data_x=x_mt, asmcmc=TRUE, figure_folder= folder_figures)
  
  
  pdf(paste0(folder_figures,'/betas_times_z_comparison.pdf'))
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
  num_mcmc <- hyperparam$num_mcmc
  folder_results <- hyperparam$folder_results
  if (!dir.exists(folder_results))  dir.create(folder_results)
  
  folder_figures <- paste0(folder_results, "/figures_model")
  if (!dir.exists(folder_figures)) dir.create(folder_figures)
  
  x <- df[, !(names(df)%in%c("target"))]
  y <- df[, (names(df)%in%c("target"))]
  x_df <- x # prepare_dataset_for_BayesGlm(x)
  x_mt <- as.matrix(x_df)
  data_train <- cbind(y, x_df)
  
  saveRDS(colnames(x_df), file.path(folder_results, "names_variables.Rds"))
  
  burn.in <- ifelse(is.null(hyperparam$burn.in), floor(num_mcmc/10), hyperparam$burn.in) 
  
  
  ### SPIKE AND SLAB model =====================================================
  mod_spike <- lm.spike(y~.-1, data = data_train, niter=num_mcmc, 
                        error.distribution = "gaussian", 
                        prior = hyperparam$priorspike)
  saveRDS(mod_spike, file = file.path(folder_results, "mod_spike.Rds"))
  cat('Saved `mod_spike` into ', folder_results, "\n")
  
  
  ### Shrinkage HORSESHOE model ================================================
  mod_hs <- bayesreg(y~., data=data_train, model ="gaussian", prior = "horseshoe", 
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
  
  burn.in <- ifelse(is.null(hyperparam$burn.in), floor(num_mcmc/10), hyperparam$burn.in) 
  
  # Standardize numerical features and perform one-hot encoding ==================
  xtrain <- train[, !(names(df)%in%c("target"))]
  ytrain <- train[, (names(df)%in%c("target"))]
  xtrain_df <- xtrain #prepare_dataset_for_BayesGlm(xtrain)
  xtrain_mt <- as.matrix(xtrain_df)
  
  xtest <- test[, !(names(df)%in%c("target"))]
  ytest <- test[, (names(df)%in%c("target"))]
  xtest_df <- xtest #prepare_dataset_for_BayesGlm(xtest)
  xtest_mt <- as.matrix(xtest_df)
  
  saveRDS(colnames(xtrain_df), file.path(folder_results, "names_variables.Rds"))
  
  
  ### GLM model ================================================================
  data_train <- cbind(ytrain, xtrain_df)
  mod_glm <- glm(ytrain ~.-1, family = gaussian, data=data_train)
  
  summary(mod_glm)
  
  
  saveRDS(mod_glm, file = file.path(folder_results, "mod_glm.Rds"))
  cat('Saved `mod_glm` into ', folder_results, "\n")
  
  
  ### BayesVS-GLM  =============================================================
  if (hyperparam$beta0 == 'glm_coefficients'){
    beta0 <- mod_glm$coefficients
  } else if (hyperparam$beta0 == 'zero'){
    beta0 <- rep(0, ncol(xtrain))
  } else {
    stop('Unknow starting point of the chain. Please check: `hyperparam$beta0`. ')
  }
  
  mod <- bayes_glmvs(xtrain_mt, ytrain, family=gaussian.bayes(), 
                     control = hyperparam)
  mod <- fit(mod, num_mcmc, beta0=beta0, burn.in = burn.in, asmcmc =F)
  saveRDS(mod, file = file.path(folder_results, "mod.Rds"))
  saveRDS(hyperparam, file = file.path(folder_results, "hyperparam.Rds"))
  cat('Saved `mod` and `hyperparam` into ', folder_results, "\n")
  
  ### Bayes-GLM  ===============================================================
  mod_z1 <- bayes_glmvs_oracle(xtrain_mt, ytrain, family=gaussian.bayes(), 
                               control = hyperparam, z_true=rep(1, ncol(xtrain_mt)))
  mod_z1 <- fit(mod_z1, num_mcmc, beta0=beta0, burn.in =burn.in, asmcmc =FALSE)
  saveRDS(mod_z1, file = file.path(folder_results, "mod_z1.Rds"))
  cat('Saved `mod_z1`  into ', folder_results, "\n")
  
  # PREDICTION with GLM and BayesGlm  ============================================
  ypred <- predict(mod, newdata=xtest_mt, type = "response")
  ypred_z1 <- predict(mod_z1, newdata=xtest_mt, type = "response")
  ypred_glm <- predict.glm(mod_glm, newdata = xtest_df, type = "response")
  
  n <- length(ytest)
  p <- ncol(xtest_mt)
  tss <- sum((ytest - mean(ytest))^2)       # Total sum of squares
  
  mae_glm <- mean(abs(ytest - ypred_glm))
  rmse_glm <- sqrt(mean((ytest - ypred_glm)^2))
  rss_glm <- sum((ytest - ypred_glm)^2)         # Residual sum of squares
  r2_glm <- 1 - rss_glm/tss
  adj_r2_glm <- 1 - (1 - r2_glm) * (n - 1)/(n - p - 1)
  
  mae <- mean(abs(ytest - ypred))
  rmse <- sqrt(mean((ytest - ypred)^2))
  rss <- sum((ytest - ypred)^2)         # Residual sum of squares
  r2 <- 1 - rss/tss
  adj_r2 <- 1 - (1 - r2) * (n - 1)/(n - p - 1)
  
  mae_z1 <- mean(abs(ytest - ypred_z1))
  rmse_z1 <- sqrt(mean((ytest - ypred_z1)^2))
  rss_z1 <- sum((ytest - ypred_z1)^2)         # Residual sum of squares
  r2_z1 <- 1 - rss_z1/tss
  adj_r2_z1 <- 1 - (1 - r2_z1) * (n - 1)/(n - p - 1)
  
  
  # Combine into a metrics matrix
  metrics <- rbind(
    MAE  = c(BayesVSGLM = mae, BayesGLM = mae_z1,  MLE = mae_glm),
    RMSE = c(BayesVSGLM = rmse, BayesGLM = rmse_z1, MLE = rmse_glm),
    RSS = c(BayesVSGLM = rss, BayesGLM = rss_z1, MLE = rss_glm),
    R2 = c(BayesVSGLM = r2, BayesGLM = r2_z1, MLE = r2_glm),
    AdjR2 = c(BayesVSGLM = adj_r2, BayesGLM = adj_r2_z1, MLE = adj_r2_glm)
  )
  colnames(metrics) <- c("BayesVS-GLM", "Bayes-GLM", "MLE")
  saveRDS(metrics, file = file.path(folder_results, "prediction_metrics.Rds"))
  
  
  ### PLOTS:  ================================================================== 
  
  plot(mod, data_x=xtrain_mt, figure_folder= folder_figures)
  plot(mod, data_x=xtrain_mt, asmcmc=TRUE, figure_folder= folder_figures)
  
  pdf(paste0(folder_figures,'/betas_times_z_comparison.pdf'))
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
  
  burn.in <- ifelse(is.null(hyperparam$burn.in), floor(num_mcmc/10), hyperparam$burn.in) 
  
  # Standardize numerical features and perform one-hot encoding ==================
  xtrain <- train[, !(names(df)%in%c("target"))]
  ytrain <- train[, (names(df)%in%c("target"))]
  xtrain_df <- xtrain #prepare_dataset_for_BayesGlm(xtrain)
  xtrain_mt <- as.matrix(xtrain_df)
  
  xtest <- test[, !(names(df)%in%c("target"))]
  ytest <- test[, (names(df)%in%c("target"))]
  xtest_df <- xtest #prepare_dataset_for_BayesGlm(xtest)
  xtest_mt <- as.matrix(xtest_df)
  
  saveRDS(colnames(xtrain_df), file.path(folder_results, "names_variables.Rds"))
  
  data_train <- cbind(ytrain, xtrain_df)

  
  ### SPIKE AND SLAB model =====================================================
  mod_spike <- lm.spike(ytrain~.-1, data =data_train, niter=num_mcmc, 
                        prior = hyperparam$priorspike,
                        error.distribution = "gaussian"
                        )
  saveRDS(mod_spike, file = file.path(folder_results, "mod_spike.Rds"))
  cat('Saved `mod_spike` into ', folder_results, "\n")
  
  
  ### Shrinkage HORSESHOE model ================================================
  mod_hs <- bayesreg(ytrain~., data=data_train, model ="gaussian", prior = "horseshoe", 
                     n.samples=num_mcmc_hs, burnin=burnin_hs)
  saveRDS(mod_hs, file = file.path(folder_results, "mod_hs.Rds"))
  cat('Saved `mod_hs` into ', folder_results, "\n")

  
  # PREDICTION   ============================================
  ypred_spike <- BoomSpikeSlab:::predict.lm.spike(mod_spike, newdata=xtest_mt, mean.only=TRUE)
  ypred_hs <- bayesreg:::predict.bayesreg(mod_hs, newdata=xtest_df, type = "response")

  n <- length(ytest)
  p <- ncol(xtest_mt)
  tss <- sum((ytest - mean(ytest))^2)       # Total sum of squares
  
  mae_spike <- mean(abs(ytest - ypred_spike))
  rmse_spike <- sqrt(mean((ytest - ypred_spike)^2))
  rss_spike <- sum((ytest - ypred_spike)^2)         # Residual sum of squares
  r2_spike <- 1 - rss_spike/tss
  adj_r2_spike <- 1 - (1 - r2_spike) * (n - 1)/(n - p - 1)
  
  mae_hs <- mean(abs(ytest - ypred_hs))
  rmse_hs <- sqrt(mean((ytest - ypred_hs)^2))
  rss_hs <- sum((ytest - ypred_hs)^2)         # Residual sum of squares
  r2_hs <- 1 - rss_hs/tss
  adj_r2_hs <- 1 - (1 - r2_hs) * (n - 1)/(n - p - 1)
  
  # Combine into a metrics matrix
  metrics <- rbind(
    MAE  = c(SpikeSlab = mae_spike, Horseshoe = mae_hs),
    RMSE = c(SpikeSlab = rmse_spike, Horseshoe = rmse_hs),
    RSS = c(SpikeSlab = rss_spike, Horseshoe = rss_hs),
    R2 = c(SpikeSlab = r2_spike, Horseshoe = r2_hs),
    AdjR2 = c(SpikeSlab = adj_r2_spike, Horseshoe = adj_r2_hs)
  )
  colnames(metrics) <- c("SpikeSlab", "Horseshoe")
  saveRDS(metrics, file = file.path(folder_results, "prediction_metrics_baselines.Rds"))
  
  
  
  ### PLOTS:  ================================================================== 
  
  pdf(paste0(folder_figures,'/beta_spike.pdf'))
  boxplot(mod_spike$beta, main='Posterior Coefficients (SpikeSlab)')
  dev.off()
  
  pdf(paste0(folder_figures,'/beta_horseshoe.pdf'))
  boxplot(t(mod_hs$beta), main='Posterior Coefficients (Horseshoe)')
  dev.off()
  
  
}

