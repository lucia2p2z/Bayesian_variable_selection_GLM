run_only_gibbs <- function(df, hyperparam) {
  
  num_mcmc <- hyperparam$num_mcmc
  folder_results <- hyperparam$folder_results
  if (!dir.exists(folder_results))  dir.create(folder_results)
  
  folder_figures <- paste0(folder_results, "/figures_model")
  if (!dir.exists(folder_figures)) dir.create(folder_figures)
  
  # Standardize numerical features and perform one-hot encoding ==================
  num_vars <- sapply(df, is.numeric)
  df[num_vars] <- scale(df[num_vars])
  
  x <- df[, !(names(df)%in%c("target"))]
  y <- df[, (names(df)%in%c("target"))]
  y <- (ifelse(y=="Disease", 1, 0))
  x_df <- prepare_dataset_for_BayesGlm(x)
  x_mt <- as.matrix(x_df)
  
  saveRDS(colnames(x_df), file.path(folder_results, "names_variables.Rds"))
  
  ### GLM model ================================================================== 
  data_train <- cbind(y, x_df)
  mod_glm <- glm(y ~.-1, family = binomial, data=data_train)
  
  saveRDS(mod_glm, file = file.path(folder_results, "mod_glm.Rds"))
  cat('Saved `mod_glm` into ', folder_results, "\n")
  
  ### BayesVS-GLM  =============================================================
  if (hyperparam$beta0 == 'glm_coefficients'){
    beta0 <- mod_glm$coefficients
  } else if (hyperparam$beta0 == 'zero'){
    beta0 <- rep(0, ncol(x_mt))
  } else {
    stop('Unknow starting point of the chain. Please check: `hyperparam$betA0`. ')
  }
  
  #set.seed(seed)
  mod <- bayes_glmvs(x_mt, y, family=binomial.bayes(), control = hyperparam)
  mod <- fit(mod, num_mcmc, beta0=beta0, burn.in =floor(num_mcmc/10), asmcmc =FALSE)
  saveRDS(mod, file = file.path(folder_results, "mod.Rds"))
  saveRDS(hyperparam, file = file.path(folder_results, "hyperparam.Rds"))
  cat('Saved `mod` and `hyperparam` into ', folder_results, "\n")
  
  ### Bayes-GLM  ===============================================================
  mod_z1 <- bayes_glmvs_oracle(x_mt, y, family=binomial.bayes(), 
                               control = hyperparam, z_true=rep(1, ncol(x_mt)))
  mod_z1 <- fit(mod_z1, num_mcmc, beta0=beta0, burn.in =floor(num_mcmc/10), asmcmc =FALSE)
  saveRDS(mod_z1, file = file.path(folder_results, "mod_z1.Rds"))
  cat('Saved `mod_z1`  into ', folder_results, "\n")
  
  
  ### PLOTS:  ================================================================== 
  
  plot(mod, data_x=x_mt, figure_folder= folder_figures)
  plot(mod, data_x=x_mt, asmcmc=TRUE, figure_folder= folder_figures)
  
  
  pdf(paste0(folder_figures,'/betas_times_z_bayesvsglm.pdf'))
  p <- plot_beta_z_boxplot(mod, mod_glm, x_df)
  print(p)
  dev.off()
  
  pdf(paste0(folder_figures,'/betas_times_z_comparison.pdf'))
  p <- plot_beta_z_boxplot_full(mod, mod_z1, mod_glm, colnames(x_df))
  print(p)
  dev.off()
  
  pdf(paste0(folder_figures,'/betas_bayesvsglm.pdf'))
  boxplot(mod@fit$beta, main = 'beta')
  dev.off()
  
  pdf(paste0(folder_figures,'/betas_bayesglm.pdf'))
  boxplot(mod_z1@fit$beta, main = 'beta')
  dev.off()
  
  #### Table Z =================================================================
  
  Z <- mod@fit$z
  ZZ <- apply(Z, 1, function(x) paste(x, collapse = ""))
  tZ <- table(ZZ)
  stZ <- sort(tZ, decreasing = TRUE)
  saveRDS(stZ, file=file.path(folder_results, "table_z.Rds"))
  
}


run_only_gibbs_baselines <- function(df, hyperparam) {
  
  num_mcmc <- hyperparam$num_mcmc
  num_mcmc_hs <- hyperparam$num_mcmc_hs
  burnin_hs <- hyperparam$burnin_hs
  
  folder_results <- hyperparam$folder_results
  if (!dir.exists(folder_results))  dir.create(folder_results)
  
  folder_figures <- paste0(folder_results, "/figures_model")
  if (!dir.exists(folder_figures)) dir.create(folder_figures)
  
  # Standardize numerical features and perform one-hot encoding ==================
  num_vars <- sapply(df, is.numeric)
  df[num_vars] <- scale(df[num_vars])
  
  x <- df[, !(names(df)%in%c("target"))]
  y <- df[, (names(df)%in%c("target"))]
  y <- (ifelse(y=="Disease", 1, 0))
  x_df <- prepare_dataset_for_BayesGlm(x)
  x_mt <- as.matrix(x_df)
  
  saveRDS(colnames(x_df), file.path(folder_results, "names_variables.Rds"))
  
  data_train <- cbind(y, x_df)
  
  ### SPIKE AND SLAB model =====================================================
  mod_spike <- logit.spike(y~.-1, data = data_train, niter=num_mcmc, 
                        prior = hyperparam$priorspike)
  saveRDS(mod_spike, file = file.path(folder_results, "mod_spike.Rds"))
  cat('Saved `mod_spike` into ', folder_results, "\n")
  
  data_train$y <- as.factor(y)
  
  ### Shrinkage HORSESHOE model ================================================
  mod_hs <- bayesreg(y~., data=data_train, model ="binomial", prior = "horseshoe", 
                     n.samples=num_mcmc_hs, burnin=burnin_hs)
  saveRDS(mod_hs, file = file.path(folder_results, "mod_hs.Rds"))
  cat('Saved `mod_hs` into ', folder_results, "\n")
  
  
  ### PLOTS:  ================================================================== 
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
  num_vars <- sapply(df, is.numeric)
  train[num_vars] <- scale(train[num_vars])
  test[num_vars] <- scale(test[num_vars])
  
  xtrain <- train[, !(names(df)%in%c("target"))]
  ytrain <- train[, (names(df)%in%c("target"))]
  ytrain <- (ifelse(ytrain=="Disease", 1, 0))
  xtrain_df <- prepare_dataset_for_BayesGlm(xtrain)
  xtrain_mt <- as.matrix(xtrain_df)
  
  xtest <- test[, !(names(df)%in%c("target"))]
  ytestfactor <- test[, (names(df)%in%c("target"))]
  ytest <- (ifelse(ytestfactor=="Disease", 1, 0))
  xtest_df <- prepare_dataset_for_BayesGlm(xtest)
  xtest_mt <- as.matrix(xtest_df)
  
  saveRDS(colnames(xtrain_df), file.path(folder_results, "names_variables.Rds"))
  
  ### GLM model ================================================================== 
  data_train <- cbind(ytrain, xtrain_df)
  mod_glm <- glm(ytrain ~.-1, family = binomial, data=data_train)
  
  saveRDS(mod_glm, file = file.path(folder_results, "mod_glm.Rds"))
  cat('Saved `mod_glm` into ', folder_results, "\n")
  
  ### BayesVS-GLM  =============================================================
  if (hyperparam$beta0 == 'glm'){
    beta0 <- mod_glm$coefficients
  } else if (hyperparam$beta0 == 'zero'){
    beta0 <- rep(0, ncol(xtrain_mt))
  } else {
    stop('Unknow starting point of the chain. Please check: `hyperparam$beta0`. ')
  }
  
  #set.seed(seed)
  mod <- bayes_glmvs(xtrain_mt, ytrain, family=binomial.bayes(), control = hyperparam)
  mod <- fit(mod, num_mcmc, beta0=beta0, burn.in =floor(num_mcmc/10), asmcmc =FALSE)
  
  if(is.null(mod@fit$stopped)) {
    saveRDS(mod, file = file.path(folder_results, "mod.Rds"))
    saveRDS(hyperparam, file = file.path(folder_results, "hyperparam.Rds"))
    cat('Saved `mod` and `hyperparam` into ', folder_results, "\n")
    
    ### Bayes-GLM  ===============================================================
    mod_z1 <- bayes_glmvs_oracle(xtrain_mt, ytrain, family=binomial.bayes(), 
                                 control = hyperparam, z_true=rep(1, ncol(xtrain_mt)))
    mod_z1 <- fit(mod_z1, num_mcmc, beta0=beta0, burn.in =floor(num_mcmc/10), asmcmc =FALSE)
    saveRDS(mod_z1, file = file.path(folder_results, "mod_z1.Rds"))
    cat('Saved `mod_z1`  into ', folder_results, "\n")
    
    # PREDICTION with GLM and BayesGlm  ============================================
    pred_prob <- predict(mod, newdata=xtest_mt, type = "response")
    ypred <- (ifelse(pred_prob>=0.5, 1, 0))
    ypred_factor <- factor(ypred, levels = c(1, 0), 
                           labels=c("Disease", "NoDisease"))
    
    pred_prob_z1 <- predict(mod_z1, newdata=xtest_mt, type = "response")
    ypred_z1 <- (ifelse(pred_prob_z1>=0.5, 1, 0))
    ypred_factor_z1 <- factor(ypred_z1, levels = c(1, 0), 
                              labels=c("Disease", "NoDisease"))
    
    
    pred_prob_glm <- predict.glm(mod_glm, newdata = xtest_df, type = "response")
    ypred_glm <- (ifelse(pred_prob_glm>=0.5, 1, 0)) 
    ypred_glm_factor <- factor(ypred_glm, levels = c(1, 0), 
                               labels=c("Disease", "NoDisease"))
    
    
    c <- confusionMatrix(data=ypred_factor, reference = ytestfactor)
    c_z1 <- confusionMatrix(data=ypred_factor_z1, reference = ytestfactor)
    c_glm <- confusionMatrix(data=ypred_glm_factor, reference = ytestfactor)
    
    metrics <- cbind(as.matrix(c, "classes"), 
                     as.matrix(c_z1, "classes"),
                     as.matrix(c_glm, "classes")) 
    colnames(metrics) <- c("BayesVS-GLM", "Bayes-GLM", "MLE")
    saveRDS(metrics, file = file.path(folder_results, "prediction_metrics.Rds"))
    
    ### PLOTS:  ================================================================== 
    
    plot(mod, data_x=xtrain_mt, figure_folder= folder_figures)
    plot(mod, data_x=xtrain_mt, asmcmc=TRUE, figure_folder= folder_figures)
  
    
    pdf(paste0(folder_figures,'/betas_times_z_bayesvsGLM.pdf'))
    p <- plot_beta_z_boxplot(mod, mod_glm, xtrain_df)
    print(p)
    dev.off()
    
    pdf(paste0(folder_figures,'/betas_times_z_comparison.pdf'))
    p <- plot_beta_z_boxplot_full(mod, mod_z1, mod_glm, colnames(xtrain_df))
    print(p)
    dev.off()
  
    pdf(paste0(folder_figures,'/betas_bayesvsglm.pdf'))
    boxplot(mod@fit$beta, main = 'beta')
    dev.off()
    
    pdf(paste0(folder_figures,'/betas_bayesglm.pdf'))
    boxplot(mod_z1@fit$beta, main = 'beta')
    dev.off()
    
    #### Table Z =================================================================
    
    Z <- mod@fit$z
    ZZ <- apply(Z, 1, function(x) paste(x, collapse = ""))
    tZ <- table(ZZ)
    stZ <- sort(tZ, decreasing = TRUE)
    saveRDS(stZ, file=file.path(folder_results, "table_z.Rds"))
    
  }
  
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
  num_vars <- sapply(df, is.numeric)
  train[num_vars] <- scale(train[num_vars])
  test[num_vars] <- scale(test[num_vars])
  
  xtrain <- train[, !(names(df)%in%c("target"))]
  ytrain <- train[, (names(df)%in%c("target"))]
  ytrain <- (ifelse(ytrain=="Disease", 1, 0))
  ytrain_factor <- as.factor(ytrain)
  xtrain_df <- prepare_dataset_for_BayesGlm(xtrain)
  xtrain_mt <- as.matrix(xtrain_df)
  
  xtest <- test[, !(names(df)%in%c("target"))]
  ytestfactor <- test[, (names(df)%in%c("target"))]
  ytest <- (ifelse(ytestfactor=="Disease", 1, 0))
  xtest_df <- prepare_dataset_for_BayesGlm(xtest)
  xtest_mt <- as.matrix(xtest_df)
  
  saveRDS(colnames(xtrain_df), file.path(folder_results, "names_variables.Rds"))
  
  data_train <- cbind(ytrain, xtrain_df)
  
  
  ### SPIKE AND SLAB model =====================================================
  mod_spike <- logit.spike(ytrain~.-1, data =data_train, niter=num_mcmc, 
                        prior = hyperparam$priorspike)
  
  saveRDS(mod_spike, file = file.path(folder_results, "mod_spike.Rds"))
  cat('Saved `mod_spike` into ', folder_results, "\n")
  
  
  ### Shrinkage HORSESHOE model ================================================
  ytrain <- as.factor(ytrain)
  data_train_factor <- cbind(ytrain, xtrain_df)
  mod_hs <- bayesreg(ytrain~., data=data_train_factor, model ="binomial", prior = "horseshoe", 
                     n.samples=num_mcmc_hs, burnin=burnin_hs)
  saveRDS(mod_hs, file = file.path(folder_results, "mod_hs.Rds"))
  cat('Saved `mod_hs` into ', folder_results, "\n")
  
    
  # PREDICTION   ===============================================================
  pred_prob_spike <- BoomSpikeSlab::predict.logit.spike(mod_spike, newdata=xtest_mt, 
                                                  type='prob', mean.only=TRUE)
  pred_prob_spike <- rowMeans(pred_prob_spike)
  ypred_spike <- (ifelse(pred_prob_spike>=0.5, 1, 0))
  ypred_spike_factor <- factor(ypred_spike, levels = c(1, 0), 
                               labels=c("Disease", "NoDisease"))
  
  pred_prob_hs <- bayesreg:::predict.bayesreg(mod_hs, newdata=xtest_df, type = "response")
  ypred_hs <- (ifelse(pred_prob_hs>=0.5, 1, 0))
  ypred_hs_factor <- factor(ypred_hs, levels = c(1, 0), 
                            labels=c("Disease", "NoDisease"))
  
  
  c_spike <- confusionMatrix(data=ypred_spike_factor, reference = ytestfactor)
  c_hs <- confusionMatrix(data=ypred_hs_factor, reference = ytestfactor)

  metrics <- cbind(as.matrix(c_spike, "classes"), 
                   as.matrix(c_hs, "classes")) 
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



