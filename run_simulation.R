simula <- function(config){

  set.seed(config$seed)
  
  start_chain <- config$start_chain
  str_n <- paste0("Nmcmc_", config$num_mcmc)
  str_a0 <- paste0("a0_", config$hyperparam$a0)
  str_y0 <- paste0("y0_", config$hyperparam$y0_mean, "_", config$hyperparam$y0_std)
  alpha <- config$hyperparam$alpha
  
  str_param <- paste0("p", config$p)
  if (config$c != 0) str_param <- paste0(str_param, "_c", config$c)
  if (config$d != 0) str_param <- paste0(str_param, "_d", config$d)
  str_param <- paste0(str_param, "_n", config$n)

  #str_param <- paste0(str_n, "_", str_a0, "_", str_y0, "_NOburnin_alpha", alpha, "/", str_param, "/seed_", config$seed)
  #str_param <- paste0(str_n, "_", str_a0, "_", str_y0, "_smallbeta", "/", str_param, "/seed_", config$seed)
  str_param <- paste0(str_n, "_", str_a0, "_", str_y0, "/", str_param, "/seed_", config$seed)
  
  # Create folder output simulations_synthetic/.../p..._n.../seed..
  outdir <- file.path(paste0("simulations_with_baselines/", config$hyperparam$method_sir, "_start_", config$start_chain), 
                      config$name, str_param)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cat("Results will be saved into: \n", outdir, "\n")
  
  # Generate data
  data <- gen_data_corr(n=config$n, p=config$p, k=config$k, c=config$c, d=config$d,
                        family=config$family,
                        true_beta=config$true_beta, true_tau=config$true_tau)
  
  # Define starting point of the Gibbs chain
  if(start_chain=='glm'){
    g <- glm(data$y ~ data$X -1 , family=config$family)
    beta0 <- g$coefficients
    
  } else if(start_chain=='zero') {
    beta0 <- rep(0, config$p) 
  }
  
  cat('\nGibbs for BayesVS-GLM: \n')
  mod <- bayes_glmvs(x = data$X, y = data$y, family = config$family, control = config$hyperparam)
  mod <- fit(mod, num_mcmc = config$num_mcmc, burn.in = round(config$num_mcmc *0.1), 
             beta0=beta0,
             asmcmc = FALSE)
  y_pred_mean <- predict(mod, newdata = data$Xtest, type = "response", full_posterior = FALSE)
  y_pred_full <- predict(mod, newdata = data$Xtest, type = "response", full_posterior = TRUE)
  predictions <- list(y_pred_mean = y_pred_mean, y_pred_full = y_pred_full)
  saveRDS(data, file = file.path(outdir, "data.Rds"))
  saveRDS(config, file = file.path(outdir, "config_used.Rds"))
  saveRDS(mod, file = file.path(outdir, "model.Rds"))
  saveRDS(predictions, file = file.path(outdir, "predictions.Rds"))
  cat('Saved `BayesVS-GLM` into ', outdir, "\n.")
  
  cat('\nGibbs for Bayes-GLM(Oracle): \n')
  mod_ora <- bayes_glmvs_oracle(x = data$X, y = data$y, z_true=data$z, family = config$family, control = config$hyperparam)
  mod_ora <- fit(mod_ora, num_mcmc = config$num_mcmc, burn.in = round(config$num_mcmc *0.1), 
                 beta0=beta0,
                 asmcmc = FALSE)
  y_pred_mean_ora <- predict(mod_ora, newdata = data$Xtest, type = "response", full_posterior = FALSE)
  y_pred_full_ora <- predict(mod_ora, newdata = data$Xtest, type = "response", full_posterior = TRUE)
  predictions_ora <- list(y_pred_mean = y_pred_mean_ora, y_pred_full = y_pred_full_ora)
  saveRDS(mod_ora, file = file.path(outdir, "model_oracle.Rds"))
  saveRDS(predictions_ora, file = file.path(outdir, "predictions_oracle.Rds"))
  cat('Saved `Bayes-GLM(Oracle)` into ', outdir, "\n")
  
  model_glm <- glm(data$y ~ data$X -1, family = config$family$family)
  saveRDS(model_glm, file = file.path(outdir, "model_glm.Rds"))
  cat('Saved `MLE` into ', outdir, "\n")
  
  
  cat('\nGibbs for Bayes-GLM: \n')
  z1 <- rep(1, config$p)
  mod_z1 <- bayes_glmvs_oracle(x = data$X, y = data$y, z_true=z1, family = config$family, control = config$hyperparam)
  mod_z1 <- fit(mod_z1, num_mcmc = config$num_mcmc, burn.in = round(config$num_mcmc *0.1), 
                beta0=beta0,
                asmcmc = FALSE)
  y_pred_mean_z1 <- predict(mod_z1, newdata = data$Xtest, type = "response", full_posterior = FALSE)
  y_pred_full_z1 <- predict(mod_z1, newdata = data$Xtest, type = "response", full_posterior = TRUE)
  predictions_z1 <- list(y_pred_mean = y_pred_mean_z1, y_pred_full = y_pred_full_z1)
  
  # Save model oracle
  saveRDS(mod_z1, file = file.path(outdir, "model_z1.Rds"))
  saveRDS(predictions_z1, file = file.path(outdir, "predictions_z1.Rds"))
  cat('Saved `Bayes-GLM` into ', outdir, "\n")
}

simula_baselines <- function(config){
  
  set.seed(config$seed)
  
  start_chain <- config$start_chain
  str_n <- paste0("Nmcmc_", config$num_mcmc)
  str_a0 <- paste0("a0_", config$hyperparam$a0)
  str_y0 <- paste0("y0_", config$hyperparam$y0_mean, "_", config$hyperparam$y0_std)
  
  str_param <- paste0("p", config$p)
  if (config$c != 0) str_param <- paste0(str_param, "_c", config$c)
  if (config$d != 0) str_param <- paste0(str_param, "_d", config$d)
  str_param <- paste0(str_param, "_n", config$n)
  
  #str_param <- paste0(str_n, "_", str_a0, "_", str_y0, "_smallbeta", "/", str_param, "/seed_", config$seed)
  str_param <- paste0(str_n, "_", str_a0, "_", str_y0, "/", str_param, "/seed_", config$seed)
  
  # Create folder output simulations_synthetic/.../p..._n.../seed..
  outdir <- file.path(paste0("simulations_with_baselines/", config$hyperparam$method_sir, "_start_", config$start_chain), 
                      config$name, str_param)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cat("Results will be saved into: \n", outdir, "\n")
  
  # Generate data
  data <- gen_data_corr(n=config$n, p=config$p, k=config$k, c=config$c, d=config$d,
                        family=config$family,
                        true_beta=config$true_beta, true_tau=config$true_tau)
  
  # Define starting point of the Gibbs chain
  if(start_chain=='glm'){
    g <- glm(data$y ~ data$X -1 , family=config$family)
    beta0 <- g$coefficients
    
  } else if(start_chain=='zero') {
    beta0 <- rep(0, config$p) 
  }
  
  cat('\nSPIKE AND SLAB model: \n')
  if (config$hyperparam_baselines$priorspike == 'poissonzellner') {
    prior <- PoissonZellnerPrior(predictors = data$X, counts = data$y)
  } else if (config$hyperparam_baselines$priorspike == 'lm') {
    prior <- SpikeSlabPrior(x = data$X, y = data$y)
  } else if (config$hyperparam_baselines$priorspike == 'logitzellner') {
    prior <- LogitZellnerPrior(predictors = data$X, successes = data$y)
  } else {
    prior <- NULL
  }
  y <- data$y
  data_train <- cbind(y, as.data.frame(data$X))
  
  if (config$family$family == "gaussian"){
    mod_spike <- lm.spike(y~.-1, data = data_train, 
                          niter=config$hyperparam_baselines$num_mcmc, 
                          prior = prior, 
                          initial.value=beta0)
    pred_prob_spike <- BoomSpikeSlab::predict.lm.spike(mod_spike, 
                                                       newdata=data$Xtest, 
                                                       type = "response")
    ypred_spike <- apply(pred_prob_spike, 1, mean)
    predictions_spike <- list(y_pred_mean = pred_prob_spike, y_pred_full = ypred_spike)
    
    cat('\nHorseshoe: \n')
    mod_hs <- bayesreg(y~., data=data_train, model ="normal", prior = "horseshoe", 
                       n.samples=config$hyperparam_baselines$num_mcmc_hs, 
                       burnin=config$hyperparam_baselines$burnin_hs)
    ypred_hs <- bayesreg:::predict.bayesreg(mod_hs, 
                                                newdata=as.data.frame(data$Xtest), 
                                                type = "response")
    predictions_hs <- list(y_pred_mean = ypred_hs, y_pred_full = NULL)
    
  } else if (config$family$family == "poisson"){
    mod_spike <- poisson.spike(y~.-1, data = data_train, 
                               niter=config$hyperparam_baselines$num_mcmc, 
                               prior = prior, 
                               initial.value = beta0)
    pred_prob_spike <- BoomSpikeSlab::predict.poisson.spike(mod_spike, 
                                                            newdata=data$Xtest, 
                                                            type = "response")
    pred_prob_spike <- apply(pred_prob_spike, 1, mean)
    ypred_spike <- rpois(length(pred_prob_spike), lambda = pred_prob_spike)
    predictions_spike <- list(y_pred_mean = ypred_spike, lambda = pred_prob_spike)
    
    cat('\nHorseshoe: \n')
    mod_hs <- bayesreg(y~., data=data_train, model ="poisson", prior = "horseshoe", 
                       n.samples=config$hyperparam_baselines$num_mcmc_hs, 
                       burnin=config$hyperparam_baselines$burnin_hs)
    pred_prob_hs <- bayesreg:::predict.bayesreg(mod_hs, 
                                                newdata=as.data.frame(data$Xtest), 
                                                type = "response")
    ypred_hs <- rpois(length(pred_prob_hs), lambda = pred_prob_hs)
    predictions_hs <- list(y_pred_mean = ypred_hs, lambda = pred_prob_hs)
    
  } else if (config$family$family == "binomial"){
    mod_spike <- logit.spike(y~.-1, data = data_train, 
                          niter=config$hyperparam_baselines$num_mcmc, 
                          prior = prior,
                          initial.value = beta0)
    pred_prob_spike <- BoomSpikeSlab::predict.logit.spike(mod_spike, 
                                                          newdata=data$Xtest, 
                                                          type = "response")
    pred_prob_spike <- apply(pred_prob_spike, 1, mean)
    ypred_spike <- rbinom(length(pred_prob_spike), size=1, prob = pred_prob_spike)
    predictions_spike <- list(y_pred_mean = ypred_spike, lambda = pred_prob_spike)
    
    cat('\nHorseshoe: \n')
    y <- as.factor(y)
    data_train <- cbind(y, as.data.frame(data$X))
    mod_hs <- bayesreg(y~., data=data_train, model ="logistic", prior = "horseshoe",
                       n.samples=config$hyperparam_baselines$num_mcmc_hs,
                       burnin=config$hyperparam_baselines$burnin_hs)
    ypred_hs <- bayesreg:::predict.bayesreg(mod_hs,
                                            newdata=as.data.frame(data$Xtest),
                                            type = "class")
    ypred_hs <- as.numeric(as.character(ypred_hs))
    predictions_hs <- list(y_pred_mean = ypred_hs, y_pred_full = NULL)
  
  } else {
    stop("GLM Family for Baselines not implemented yet")
  }
  
  saveRDS(config, file = file.path(outdir, "config_baselines_used.Rds"))
  saveRDS(mod_spike, file = file.path(outdir, "model_spike.Rds"))
  saveRDS(predictions_spike, file = file.path(outdir, "predictions_spike.Rds"))
  cat('Saved `SpikeSlab` into ', outdir, "\n.")
  
  saveRDS(mod_hs, file = file.path(outdir, "model_hs.Rds"))
  saveRDS(predictions_hs, file = file.path(outdir, "predictions_hs.Rds"))
  cat('Saved `Horseshoe` into ', outdir, "\n")
  
}

