gen_data_gamma <- function(n, gamma, family, true_beta, true_tau, 
                          phi_noise=0.001, delta_noise=5){
  
  # Assert to catch errors
  if (!is.numeric(n) || n <= 0) stop("`n` must be a positive integer.")
  if (is.null(true_tau) || true_tau <= 0) stop("`true_tau` must be a positive number.")
  
  # Define full sample size
  n_test <- floor(0.25 * n)
  N <- n + n_test
  
  namefamily <- family$family
  
  cat("Asking to generate data as:  \n", 
      "y ~ f(eta), eta = X1*beta \n",
      "X = (X1, X2) \n",
      "X2 = gamma*X1 + (1-gamma)*W \n",
      "With parameter gamma = ", gamma, "\n\n")
  
  cat("Total", N, "observations: ", n, "train +", n_test, "test.\n\n")
  
  cat("Model family :", namefamily)
  
  p <- 2 
  k <- 1 
  
  X <- matrix(NA, nrow = N, ncol = p)
  e <- matrix(NA, nrow = N, ncol = p)
  z <- rep(0, p)
  beta <- rep(0, p)
  
  L <- matrix(rnorm(N * k), nrow = N)
  
  # Compute eta and response variable 
  eta <- L %*% true_beta
  r <- generate_y_from_eta(eta, family, true_tau)
  
  # Target covariates
  e[, 1] <- matrix(rnorm(N * k, sd = sqrt(phi_noise)), nrow = N)
  e[, 2] <- matrix(rnorm(N * k), nrow = N)
                   
  X[, 1] <- L #+ e[, 1]
  X[, 2] <- gamma*L + (1-gamma)*e[, 2]
  
  z[1:k]  <- rep(1, k)
  
  beta[1:k] <- true_beta

  # Split train/test
  Xtest <- X[(n+1):N, ]
  ytest <- r[(n+1):N]
  
  X <- X[1:n, ]
  y <- r[1:n]
  
  return(list(
    X = X,
    Xtest = Xtest,
    L = L,
    beta = beta,
    z = z,
    y = y,
    ytest = ytest,
    tau = true_tau,  
    k = k,
    gamma = gamma
  ))
}

format_gamma <- function(gamma) {
  if (gamma %% 1 == 0) {
    paste0("gamma", as.integer(gamma))
  } else {
    paste0("gamma", sub("0\\.", "0", as.character(gamma)))
  }
}

simula <- function(config){

  set.seed(config$seed)
  
  start_chain <- config$start_chain
  str_nmcmc <- paste0("Nmcmc_", config$num_mcmc)
  str_alpha <- paste0("alpha_", config$hyperparam$alpha)
  str_a0 <- paste0("a0_", config$hyperparam$a0)
  str_y0 <- paste0("y0_", config$hyperparam$y0_mean, "_", config$hyperparam$y0_std)
  str_n <- paste0("n", config$n)
  str_param <- paste0(str_nmcmc, "_", str_alpha, "_", str_a0, "_", str_y0, "_", str_n)
  
  str_gamma <- format_gamma(config$gamma)

  # Create folder output simulations_correlation/family/params../gamma../seed..
  outdir <- file.path("simulations_correlation_prova", 
                      paste0(config$hyperparam$method_sir, "_", config$hyperparam$proposal,"_start_", config$start_chain), 
                      config$name, 
                      str_param,
                      str_gamma,
                      paste0("seed_", config$seed))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cat("Results will be saved into: \n", outdir, "\n")
  
  # Generate data 
  data <- gen_data_gamma(n=config$n, gamma=config$gamma, family=config$family,
                        true_beta=config$true_beta, true_tau=config$true_tau)
  
  # Define starting point of the Gibbs chain
  if(start_chain=='glm'){
    g <- glm(data$y ~ data$X -1 , family=config$family)
    beta0 <- g$coefficients
  } else if(start_chain=='zero') {
    beta0 <- rep(0, config$p) 
  }

  ### BayesVS-GLM  ============================================================= 
  
  cat('\nGibbs for BayesVS-GLM: \n')
  mod <- bayes_glmvs(x = data$X, y = data$y, family = config$family, control = config$hyperparam)
  mod <- fit(mod, num_mcmc = config$num_mcmc, burn.in = floor(config$num_mcmc/10), 
             beta0=beta0,
             asmcmc = FALSE)
  
  if(is.null(mod@fit$stopped)) {
    
    y_pred_mean <- predict(mod, newdata = data$Xtest, type = "response", full_posterior = FALSE)
    y_pred_full <- predict(mod, newdata = data$Xtest, type = "response", full_posterior = TRUE)
    predictions <- list(y_pred_mean = y_pred_mean, y_pred_full = y_pred_full)
    saveRDS(data, file = file.path(outdir, "data.Rds"))
    saveRDS(config, file = file.path(outdir, "config_used.Rds"))
    saveRDS(mod, file = file.path(outdir, "model.Rds"))
    saveRDS(predictions, file = file.path(outdir, "predictions.Rds"))
    cat('Saved `BayesVS-GLM` into ', outdir, "\n.")
    
    
    ### Bayes-GLM(Oracle)  ======================================================= 
  
    # cat('\nGibbs for Bayes-GLM(Oracle): \n')
    # mod_ora <- bayes_glmvs_oracle(x = data$X, y = data$y, z_true=data$z, family = config$family, control = config$hyperparam)
    # mod_ora <- fit(mod_ora, num_mcmc = config$num_mcmc, burn.in = floor(config$num_mcmc/10), 
    #                beta0=beta0,
    #                asmcmc = FALSE)
    # y_pred_mean_ora <- predict(mod_ora, newdata = data$Xtest, type = "response", full_posterior = FALSE)
    # y_pred_full_ora <- predict(mod_ora, newdata = data$Xtest, type = "response", full_posterior = TRUE)
    # predictions_ora <- list(y_pred_mean = y_pred_mean_ora, y_pred_full = y_pred_full_ora)
    # saveRDS(mod_ora, file = file.path(outdir, "model_oracle.Rds"))
    # saveRDS(predictions_ora, file = file.path(outdir, "predictions_oracle.Rds"))
    # cat('Saved `Bayes-GLM(Oracle)` into ', outdir, "\n")
    
    mod_glm <- glm(data$y ~ data$X -1, family = config$family$family)
    saveRDS(mod_glm, file = file.path(outdir, "model_glm.Rds"))
    cat('Saved `MLE` into ', outdir, "\n")
    
    
    ### Bayes-GLM  ===============================================================
    cat('\nGibbs for Bayes-GLM: \n')
    z1 <- rep(1, config$p)
    mod_z1 <- bayes_glmvs_oracle(x = data$X, y = data$y, z_true=z1, family = config$family, control = config$hyperparam)
    mod_z1 <- fit(mod_z1, num_mcmc = config$num_mcmc, burn.in = floor(config$num_mcmc/10), 
                  beta0=beta0,
                  asmcmc = FALSE)
    y_pred_mean_z1 <- predict(mod_z1, newdata = data$Xtest, type = "response", full_posterior = FALSE)
    y_pred_full_z1 <- predict(mod_z1, newdata = data$Xtest, type = "response", full_posterior = TRUE)
    predictions_z1 <- list(y_pred_mean = y_pred_mean_z1, y_pred_full = y_pred_full_z1)
    
    # Save model oracle
    saveRDS(mod_z1, file = file.path(outdir, "model_z1.Rds"))
    saveRDS(predictions_z1, file = file.path(outdir, "predictions_z1.Rds"))
    cat('Saved `Bayes-GLM` into ', outdir, "\n")
    
    
    #### Table Z =================================================================
    
    Z <- mod@fit$z
    ZZ <- apply(Z, 1, function(x) paste(x, collapse = ""))
    tZ <- table(ZZ)
    stZ <- sort(tZ, decreasing = TRUE)
    saveRDS(stZ, file=file.path(outdir, "table_z.Rds"))
  } 
}

