library(ggplot2)
library(tidyverse)
library(fillpattern)

# palette.colors(palette = "R4")
#"#000000" "#DF536B" "#61D04F" "#2297E6" "#28E2E5" "#CD0BBC" "#F5C710" "#9E9E9E"
model_colors <- c(
  "BayesVS-GLM"       = "#2297E6",
  "Bayes-GLM(Oracle)" = "#28E2E5",
  "Bayes-GLM"         = "#61D04F",
  "MLE"               = "#CD0BBC",  
  "True beta"         = "#DF536B",
  "SpikeSlab"         = "#F28E2B",
  "Horseshoe"         = "#9467BD",
  "Horseshoe+"        = "#E07B39"
)
model_linetypes <- c(
  "BayesVS-GLM"       = "solid",
  "Bayes-GLM(Oracle)" = "dashed",
  "MLE"               = "dotted",
  "SpikeSlab"         = "longdash",
  "Horseshoe"         = "twodash",
  "Horseshoe+"        = "dotdash"
)

model_shapes <- c(
  "BayesVS-GLM"        = 16,  # solid circle
  "Bayes-GLM(Oracle)"  = 17,  # solid triangle
  "MLE"                = 18,  # cross
  "SpikeSlab"          = 15,  # quadrato pieno
  "Horseshoe"          = 2,   # triangolo vuoto
  "Horseshoe+"         = 3    # plus (+)
)

model_size <- c(
  "BayesVS-GLM"        = 4.5,
  "Bayes-GLM(Oracle)"  = 4,
  "MLE"                = 3,    
  "SpikeSlab"          = 3,   
  "Horseshoe"          = 3,
  "Horseshoe+"         = 3
)

load_simulation <- function(folder_path) {
  list(
    path = folder_path,
    data = readRDS(file.path(folder_path, "data.Rds")),
    config = readRDS(file.path(folder_path, "config_used.Rds")),
    model = readRDS(file.path(folder_path, "model.Rds")),
    model_oracle = readRDS(file.path(folder_path, "model_oracle.Rds")),
    predictions = readRDS(file.path(folder_path, "predictions.Rds")),
    predictions_oracle = readRDS(file.path(folder_path, "predictions_oracle.Rds")),
    model_glm = readRDS(file.path(folder_path, "model_glm.Rds")),
    model_z1 = readRDS(file.path(folder_path, "model_z1.Rds")),
    predictions_z1 = readRDS(file.path(folder_path, "predictions_z1.Rds")),
    model_spike = readRDS(file.path(folder_path, "model_spike.Rds")),
    predictions_spike = readRDS(file.path(folder_path, "predictions_spike.Rds")),
    model_hs = readRDS(file.path(folder_path, "model_hs.Rds")),
    predictions_hs = readRDS(file.path(folder_path, "predictions_hs.Rds"))
  )
}


compute_metrics_beta <- function(model, data, which_idx) {
  beta_true <- data$beta

  # if (inherits(model, c("bayes_glmvs", "bayes_glmvs_oracle"))) {
  #   beta_post <- model@fit$beta
  
  if (inherits(model, "bayes_glmvs")) {
    if(inherits(model, "bayes_glmvs_oracle")) {
      z <- t(replicate(nrow(model@fit$beta), model@z_true))
    } else {
      z <- model@fit$z
    }
    beta_post <- model@fit$beta * z
  } else if (inherits(model, "bayesreg")) {
    beta_post <- t(model$beta)
  } else if (inherits(model, "glm.spike")) {
    beta_post <- model$beta
  } else {
    stop(paste0("Unknown type of model: ", paste(class(model), collapse = ", ")))
  }
  
  beta_mean <- apply(beta_post, 2, mean)
  beta_lower <- apply(beta_post, 2, quantile, 0.025)
  beta_upper <- apply(beta_post, 2, quantile, 0.975)
  
  
  if(which_idx == "z1") { 
    # Consider only active (non-zero) coefficients
    active_idx <- which(beta_true != 0)
  
    tibble(
      r_mse_beta1 = mean(((beta_mean[active_idx] - beta_true[active_idx]) / beta_true[active_idx])^2),
      r_mae_beta1 = mean(abs((beta_mean[active_idx] - beta_true[active_idx]) / beta_true[active_idx])),
      r_maxse_beta1 = max(((beta_mean[active_idx] - beta_true[active_idx]) / beta_true[active_idx])^2),
      r_maxae_beta1 = max(abs((beta_mean[active_idx] - beta_true[active_idx]) / beta_true[active_idx])),
      coverage_beta1 = mean(beta_true[active_idx] >= beta_lower[active_idx] & 
                             beta_true[active_idx] <= beta_upper[active_idx])
    )
  } else if(which_idx == "all") {
    # Consider all coefficients
    active_idx <- rep(1, length(beta_true))
    
    tibble(
      mse_beta = mean((beta_mean[active_idx] - beta_true[active_idx] )^2),
      mae_beta = mean(abs(beta_mean[active_idx] - beta_true[active_idx])),
      maxse_beta = max((beta_mean[active_idx] - beta_true[active_idx])^2),
      maxae_beta = max(abs(beta_mean[active_idx] - beta_true[active_idx])),
      coverage_beta = mean(beta_true[active_idx] >= beta_lower[active_idx] & 
                             beta_true[active_idx] <= beta_upper[active_idx])
    )
    
  } 
}



compute_metrics_z <- function(model, data) {
  
  z_true <- data$z
  
  # Return NA if it's not a class of BayesGlm with the z nor a SpikeSLab
  #if (!(inherits(model, c("bayes_glmvs", "bayes_glmvs_oracle"))) || is.null(model@fit$z)) return(tibble(acc_z = NA_real_))
  
  if (inherits(model, c("bayes_glmvs", "bayes_glmvs_oracle"))) {
    if (is.null(model@fit$z)) return(tibble(acc_z = NA_real_))
    z_post <- model@fit$z
  } else if (inherits(model, c("poisson.spike", "lm.spike", "logit.spike"))) {
    
    z_post <- (model$beta != 0) * 1
  } else {
    return(tibble(acc_z = NA_real_))
  }
 
  z_comp <- matrix(z_true, nrow = nrow(z_post), ncol = length(z_true), byrow = TRUE)
  acc_j <- colMeans(z_post == z_comp)
  tibble(acc_z = mean(acc_j))
}


compute_metrics_tau <- function(model, data) {
  
  # Return NA if it's not a class of BayesGlm 
  
  if (!(inherits(model, c("bayes_glmvs", "bayes_glmvs_oracle"))) || is.null(model@fit$tau)) return(tibble())  # o NA
  tau_true <- data$tau
  tau_post <- model@fit$tau
  tibble(
    tau_mean = mean(tau_post),
    tau_sd = sd(tau_post),
    r_mse_tau = ((mean(tau_post) - tau_true)/tau_true)^2,
    r_mae_tau = abs((mean(tau_post) - tau_true)/tau_true)
  )
}


compute_metrics_prediction <- function(predictions, data) {
  
  if (is.null(predictions$y_pred_mean)) return(tibble())
  
  y_test <- data$ytest
  y_pred_mean <- predictions$y_pred_mean
  tibble(
    r_mse_y = mean(((y_pred_mean - y_test)/y_test)^2),
    r_mae_y = mean(abs((y_pred_mean - y_test)/y_test))
  )
}


compute_metrics_single_model <- function(model, predictions, data, config, model_name) {
  
  beta1_metrics <- compute_metrics_beta(model, data, which_idx="z1")
  beta_metrics <- compute_metrics_beta(model, data, which_idx="all")
  z_metrics <- compute_metrics_z(model, data)
  #tau_metrics <- compute_metrics_tau(model, data)
  pred_metrics <- compute_metrics_prediction(predictions, data)
  
  base_info <- tibble(
    model_type = model_name,
    name = config$name,
    n = config$n,
    p = config$p,
    n_over_p = config$n / config$p
  )
  
  #error_metrics <- bind_cols(base_info, beta_metrics, z_metrics, tau_metrics, pred_metrics)
  error_metrics <- bind_cols(base_info, beta1_metrics, beta_metrics, z_metrics, pred_metrics)
  # 
  # if (inherits(model, c("bayes_glmvs", "bayes_glmvs_oracle"))) {
  #   beta_post <- model@fit$beta
  # } else if (inherits(model, "bayesreg")) {
  #   beta_post <- t(model$beta)
  # } else if (inherits(model, "glm.spike")) {
  #   beta_post <- model$beta
  # } else {
  #   stop(paste0("Unknown type of model: ", paste(class(model), collapse = ", ")))
  # }
  # beta_df <- tibble(
  #   sim_name = config$name,
  #   model_type = model_name,
  #   n = config$n,
  #   p = config$p,
  #   beta_mean = list(apply(beta_post, 2, mean)),
  #   beta_median = list(apply(beta_post, 2, median)),
  #   beta_q25 = list(apply(beta_post, 2, quantile, 0.25)),
  #   beta_q75 = list(apply(beta_post, 2, quantile, 0.75)),
  #   beta_sd = list(apply(beta_post, 2, sd)),
  #   beta_true = list(data$beta)
  # )
  # 
  # list(metrics = error_metrics, beta = beta_df)
  list(metrics = error_metrics)
}


compute_metrics_glm <- function(glm_model, glm_predictions=NULL, data, config){
  
  base_info <- tibble(
    model_type = "MLE",
    name = config$name,
    n = config$n,
    p = config$p,
    n_over_p = config$n / config$p
  )
  
  beta_true <- data$beta
  beta_est <- glm_model$coefficients
  vcov_mat <- vcov(glm_model)
  se <- sqrt(diag(vcov_mat)) #standard error
  
  beta_lower <- beta_est - 1.96*se #apply(beta_post, 2, quantile, 0.025)
  beta_upper <- beta_est + 1.96*se #apply(beta_post, 2, quantile, 0.975)
  
  # Consider only active (non-zero) coefficients
  active_idx <- which(beta_true != 0)
  beta1_metrics <- tibble(
    r_mse_beta1 = mean(((beta_est[active_idx] - beta_true[active_idx]) / beta_true[active_idx])^2),
    r_mae_beta1 = mean(abs((beta_est[active_idx] - beta_true[active_idx]) / beta_true[active_idx])),
    r_maxse_beta1 = max(((beta_est[active_idx] - beta_true[active_idx]) / beta_true[active_idx])^2),
    r_maxae_beta1 = max(abs((beta_est[active_idx] - beta_true[active_idx]) / beta_true[active_idx])),
    coverage_beta1 = mean(beta_true[active_idx] >= beta_lower[active_idx] & 
                           beta_true[active_idx] <= beta_upper[active_idx])
    )
  
  # Consider all coefficients
  active_idx <- rep(1, length(beta_true))
  beta_metrics <- tibble(
    mse_beta = mean(((beta_est[active_idx] - beta_true[active_idx]) )^2),
    mae_beta = mean(abs((beta_est[active_idx] - beta_true[active_idx]) )),
    maxse_beta = max(((beta_est[active_idx] - beta_true[active_idx]) )^2),
    maxae_beta = max(abs((beta_est[active_idx] - beta_true[active_idx]) )),
    coverage_beta = mean(beta_true[active_idx] >= beta_lower[active_idx] & 
                           beta_true[active_idx] <= beta_upper[active_idx])
  )
  
  # TODO pred_metrics <- tibble()
  error_metrics <- bind_cols(base_info, beta1_metrics, beta_metrics) #, pred_metrics)
  list(metrics = error_metrics)
}


compute_metrics_from_sim <- function(sim) {
  
  main <- compute_metrics_single_model(sim$model, sim$predictions, 
                                       sim$data, sim$config, "BayesVS-GLM")
  oracle <- compute_metrics_single_model(sim$model_oracle, sim$predictions_oracle, 
                                         sim$data, sim$config, "Bayes-GLM(Oracle)")
  full <- compute_metrics_single_model(sim$model_z1, sim$predictions_z1, 
                                         sim$data, sim$config, "Bayes-GLM")
  
  glm <- compute_metrics_glm(sim$model_glm, NULL,sim$data, sim$config)
  
  spike <- compute_metrics_single_model(sim$model_spike, sim$predictions_spike, 
                                       sim$data, sim$config, "SpikeSlab")
  
  hs <- compute_metrics_single_model(sim$model_hs, sim$predictions_hs, 
                                        sim$data, sim$config, "Horseshoe")
  
  list(
    metrics = bind_rows(main$metrics, oracle$metrics, 
                        glm$metrics, spike$metrics, hs$metrics)
    #beta = bind_rows(main$beta, oracle$beta, spike$beta, hs$beta)
  )
}


plot_beta_boxplot <- function(sim, which_idx="all") {
  
  # Load file
  model <- sim$model
  data <- sim$data
  model2 <- sim$model_z1
  
  # Info 
  beta_true <- data$beta
  beta_post <- model@fit$beta
  beta_post2 <- model2@fit$beta
  
  # To have boxplot only for components with z = 1 ------------------------ 
  # -----------------------------------------------------------------------
  # Determine which indices to include
  idx <- rep(TRUE, length(beta_true))
  if(which_idx=="z1") {
    idx <- beta_true != 0
  }
  
  # Subset beta_true and beta_post and MLE
  beta_true <- beta_true[idx]
  beta_post <- beta_post[, idx, drop = FALSE]
  beta_post2 <- beta_post2[, idx, drop = FALSE]
  beta_glm <- sim$model_glm$coefficients[idx]
  
  # -----------------------------------------------------------------------
  
  # Prepare for boxplot
  p <- length(beta_true)
  beta_levels <- 1:p #paste0("beta", 1:p)
  colnames(beta_post) <- beta_levels
  colnames(beta_post2) <- beta_levels
  beta_labels_expr <- beta_levels #  lapply(which(idx), function(j) bquote(beta[.(j)]))  # for axis
  
  df_long1 <- pivot_longer(as.data.frame(beta_post), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "BayesVS-GLM")
  
  df_long2 <- pivot_longer(as.data.frame(beta_post2), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "Bayes-GLM")
  
  df_long_all <- bind_rows(df_long1, df_long2)
  df_long_all$beta <- factor(df_long_all$beta, levels = beta_levels)
  
  #beta_true_df <- data.frame(beta = 1:p, value = beta_true)
  #beta_glm_df <- data.frame(beta = 1:p, value = beta_glm)
  beta_true_df <- data.frame(
    beta = rep(1:p, 2),
    value = rep(beta_true, 2),
    posterior = c(rep("BayesVS-GLM", p), rep("Bayes-GLM", p))
  )
  beta_glm_df <- data.frame(
    beta = rep(1:p, 2),
    value = rep(beta_glm, 2),
    posterior = c(rep("BayesVS-GLM", p), rep("Bayes-GLM", p))
  )
  beta_mean_df <- data.frame(beta = 1:p, value = apply(beta_post, 2, mean))
  beta_mean_df2 <- data.frame(beta = 1:p, value = apply(beta_post2, 2, mean))
  
  beta_mean_df$posterior <- "BayesVS-GLM"
  beta_mean_df2$posterior <- "Bayes-GLM"
  
  beta_means_all <- bind_rows(beta_mean_df, beta_mean_df2)
  rownames(beta_means_all) <- NULL
  boxplot <- ggplot(df_long_all, aes(x = beta, y = value , fill=posterior)) +
    geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.6), width = 0.4, linewidth = 0.3) +  
    geom_point(data = beta_true_df, aes(x = beta, y = value, group = posterior, color = "True beta"), 
               position = position_dodge(width = 0.6), shape = 18, size=2) +
    geom_point(data = beta_glm_df, aes(x = beta, y = value, group = posterior, color = "MLE"), 
               position = position_dodge(width = 0.6), shape = 8, size=2) +
    theme_minimal(base_size = 14) +
    scale_fill_manual(
      name = "Posterior",
      values = c("BayesVS-GLM" = model_colors[["BayesVS-GLM"]],
                 "Bayes-GLM" = model_colors[["Bayes-GLM"]]),
      guide = guide_legend(override.aes = list(shape = NA))
    ) +
    scale_color_manual(
      name = "Reference",
      values = c("MLE" = model_colors[["MLE"]],
                 "True beta" = model_colors[["True beta"]])
    ) +
    labs(x = "Component j", y = "Value")
  
  
  return(boxplot)
  
}


plot_beta_z_boxplot <- function(sim, which_idx="all", set_ylim=NULL, set_legend=NULL) {
  
  # Load file
  model <- sim$model
  data <- sim$data
  model2 <- sim$model_z1
  
  # Info 
  beta_true <- data$beta
  beta_post <- model@fit$beta
  z_post <- model@fit$z
  beta_post2 <- model2@fit$beta
  beta_per_z <- beta_post*z_post
  
  p <- length(beta_true)
  beta_labels_expr <- 1:p
  
  # To have boxplot only for components with z = 1 ------------------------ 
  # Determine which indices to include
  
  if (which_idx=="z1") {
    idx <- beta_true != 0
  } else if (which_idx=="z0") {
    idx <- beta_true == 0
  } else if (which_idx=="all") {
    idx <- rep(TRUE, length(beta_true))
  }
  
  # Subset beta_true and beta_post and MLE
  beta_true <- beta_true[idx]
  beta_post2 <- beta_post2[, idx, drop = FALSE]
  beta_per_z <- beta_per_z[, idx, drop = FALSE]
  beta_glm <- sim$model_glm$coefficients[idx]
  beta_labels_expr <- beta_labels_expr[idx]
  
  p <- length(beta_true)
  beta_levels <- 1:p 
  
  # -----------------------------------------------------------------------
  
  # Prepare for boxplot
  colnames(beta_per_z) <- beta_levels
  colnames(beta_post2) <- beta_levels
  
  
  df_long1 <- pivot_longer(as.data.frame(beta_per_z), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "BayesVS-GLM")
  
  df_long2 <- pivot_longer(as.data.frame(beta_post2), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "Bayes-GLM")
  
  df_long_all <- bind_rows(df_long1, df_long2)
  df_long_all$beta <- factor(df_long_all$beta, levels = beta_levels)
  
  #beta_true_df <- data.frame(beta = 1:p, value = beta_true)
  #beta_glm_df <- data.frame(beta = 1:p, value = beta_glm)
  beta_true_df <- data.frame(
    beta = rep(beta_levels, 2),
    value = rep(beta_true, 2),
    posterior = c(rep("BayesVS-GLM", p), rep("Bayes-GLM", p))
  )
  beta_glm_df <- data.frame(
    beta = rep(beta_levels, 2),
    value = rep(beta_glm, 2),
    posterior = c(rep("BayesVS-GLM", p), rep("Bayes-GLM", p))
  )
  
  
  beta_mean_df <- data.frame(beta = beta_levels, value = apply(beta_per_z, 2, mean))
  beta_mean_df2 <- data.frame(beta = beta_levels, value = apply(beta_post2, 2, mean))
  beta_mean_df$posterior <- "BayesVS-GLM"
  beta_mean_df2$posterior <- "Bayes-GLM"
  beta_means_all <- bind_rows(beta_mean_df, beta_mean_df2)
  rownames(beta_means_all) <- NULL
  
  
  boxplot <- ggplot(df_long_all, aes(x = beta, y = value , fill=posterior)) +
    geom_boxplot(position = position_dodge(width = 0.8), 
                 width = 0.6, linewidth = 0.3,
                 outlier.colour = "#000000", outlier.size = 1, outlier.alpha = 0.1) +  
    geom_point(data = beta_true_df, aes(x = beta, y = value, group = posterior, color = "True beta"), 
               position = position_dodge(width = 0.8), shape = 8, size=2) +
    geom_point(data = beta_glm_df, aes(x = beta, y = value, group = posterior, color = "MLE"), 
               position = position_dodge(width = 0.8), shape = 18, size=2) +
    scale_fill_manual(
      name = "Posterior",
      values = c("BayesVS-GLM" = model_colors[["BayesVS-GLM"]],
                 "Bayes-GLM" = model_colors[["Bayes-GLM"]]),
      guide = guide_legend(override.aes = list(shape = NA))
    ) +
    scale_color_manual(
      name = "Reference",
      values = c("MLE" = model_colors[["MLE"]],
                 "True beta" = model_colors[["True beta"]])
    ) +
    scale_x_discrete(labels = beta_labels_expr) + 
    labs(x = "Component j", y = "Value") +
    theme_minimal(
      base_size = 14
    )
  
  if (!is.null(set_ylim)) {
    boxplot <- boxplot + coord_cartesian(ylim = set_ylim)
  }
  
  # Optional legend control
  if (!is.null(set_legend)) {
    if (set_legend == "none") {
      boxplot <- boxplot + theme(
        legend.position = "none"        
        )
    } else if (set_legend == "inner") {
      boxplot <- boxplot + theme(
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.box.background = element_blank()
      )
    }
  }
  
  return(boxplot)
  
}


plot_beta_z_boxplot_with_baselines <- function(sim, which_idx="all", set_ylim=NULL, set_legend=NULL) {
  
  # Load file
  model <- sim$model
  data <- sim$data
  model2 <- sim$model_z1
  
  # Info 
  beta_true <- data$beta
  
  # BayesVS-GLM
  beta_post <- model@fit$beta
  z_post <- model@fit$z
  beta_per_z <- beta_post*z_post
  
  # Bayes-GLM
  beta_post2 <- model2@fit$beta
  
  # SpikeSlab & Horseshoe
  beta_post_baseline1 <- sim$model_spike$beta
  beta_post_baseline2 <- t(sim$model_hs$beta)
  
  num_models <- 4
  
  p <- length(beta_true)
  beta_labels_expr <- 1:p
  
  # To have boxplot only for components with z = 1 ------------------------ 
  # Determine which indices to include
  
  if (which_idx=="z1") {
    idx <- beta_true != 0
  } else if (which_idx=="z0") {
    idx <- beta_true == 0
  } else if (which_idx=="all") {
    idx <- rep(TRUE, length(beta_true))
  }
  
  # Subset beta_true and beta_post and MLE
  beta_true <- beta_true[idx]
  beta_post2 <- beta_post2[, idx, drop = FALSE]
  beta_per_z <- beta_per_z[, idx, drop = FALSE]
  beta_post_baseline1 <- beta_post_baseline1[, idx, drop =FALSE] 
  beta_post_baseline2 <- beta_post_baseline2[, idx, drop =FALSE]
  beta_glm <- sim$model_glm$coefficients[idx]
  beta_labels_expr <- beta_labels_expr[idx]
  
  p <- length(beta_true)
  beta_levels <- 1:p 
  
  # -----------------------------------------------------------------------
  
  # Prepare for boxplot
  colnames(beta_per_z) <- beta_levels
  colnames(beta_post2) <- beta_levels
  colnames(beta_post_baseline1) <- beta_levels
  colnames(beta_post_baseline2) <- beta_levels
  
  df_long1 <- pivot_longer(as.data.frame(beta_per_z), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "BayesVS-GLM")
  
  df_long2 <- pivot_longer(as.data.frame(beta_post2), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "Bayes-GLM")
  
  df_long3 <- pivot_longer(as.data.frame(beta_post_baseline1), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "SpikeSlab")
  
  df_long4 <- pivot_longer(as.data.frame(beta_post_baseline2), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "Horseshoe")
  
  
  df_long_all <- bind_rows(df_long1, df_long2, df_long3, df_long4) 
  df_long_all$beta <- factor(df_long_all$beta, levels = beta_levels)
  
  beta_true_df <- data.frame(
    beta = rep(beta_levels, num_models),
    value = rep(beta_true, num_models),
    posterior = c(rep("BayesVS-GLM", p), 
                  rep("Bayes-GLM", p), 
                  rep("SpikeSlab", p),
                  rep("Horseshoe", p))
  )
  beta_glm_df <- data.frame(
    beta = rep(beta_levels, num_models),
    value = rep(beta_glm, num_models),
    posterior = c(rep("BayesVS-GLM", p), 
                  rep("Bayes-GLM", p), 
                  rep("SpikeSlab", p),
                  rep("Horseshoe", p))
  )
  
  posterior_levels <- c("BayesVS-GLM", "Bayes-GLM", "SpikeSlab", "Horseshoe")
  legend_breaks    <- c(posterior_levels, "MLE", "True beta")
  df_long_all$posterior  <- factor(df_long_all$posterior,  levels = posterior_levels)
  beta_true_df$posterior <- factor(beta_true_df$posterior, levels = posterior_levels)
  beta_glm_df$posterior  <- factor(beta_glm_df$posterior,  levels = posterior_levels)
  
  dodge <- 0.9
  
  boxplot <- ggplot(df_long_all, aes(x = beta, y = value , fill=posterior)) +
    geom_boxplot(aes(color=after_scale(fill)),
                 position = position_dodge(width = dodge), alpha=0.7,
                 width = 0.7, linewidth = 0.5,
                 #outlier.colour = "#000000", 
                 outlier.size = 1, outlier.alpha = 0.2) +  
    geom_point(data = beta_true_df, 
               aes(x = beta, y = value, group = posterior, shape="True beta"), 
               position = position_dodge(width = dodge), color = model_colors[["True beta"]], size=2) +
    geom_point(data = beta_glm_df, 
               aes(x = beta, y = value, group = posterior, shape="MLE"), 
               position = position_dodge(width = dodge), color = model_colors[["MLE"]], size=2) +
    theme_minimal(base_size = 14) +
    scale_fill_manual(
      name = "Posterior",
      breaks = posterior_levels,
      values = c(model_colors[posterior_levels])
    ) +
    scale_shape_manual(
      name   = "Reference",
      breaks = c("MLE", "True beta"),
      values = c("MLE" = 18, "True beta" = 8)
    ) +
    guides(
      fill  = guide_legend(order = 1, override.aes = list(shape = NA))#,
      #color = guide_legend(order = 1, override.aes = list(shape = c(NA, NA, NA, NA, 18, 8)))
    ) +
    scale_x_discrete(labels = beta_labels_expr) + 
    labs(x = "Component j", y = "Value")
  
  if (!is.null(set_ylim)) {
    boxplot <- boxplot + coord_cartesian(ylim = set_ylim)
  }
  
  # Optional legend control
  if (!is.null(set_legend)) {
    if (set_legend == "none") {
      boxplot <- boxplot + theme(legend.position = "none")
    } else if (set_legend == "inner") {
      boxplot <- boxplot + theme(
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.box.background = element_blank()
      )
    }
  }
  
  return(boxplot)
  
}


get_accuracy_z_df <- function(sim){
  z_post <- sim$model@fit$z
  z_true <- sim$data$z
  p <- sim$config$p
  z_comp <- matrix(z_true, nrow = nrow(z_post), ncol = length(z_true), byrow = TRUE)
  acc_j <- colMeans(z_post == z_comp)
  tibble(
    index = 1:p,
    z_true = factor(z_true, levels = c(0, 1)),
    acc_j = acc_j,
    seed = sim$config$seed,
    n = sim$config$n,
    p = sim$config$p
  )
}

get_accuracy_zspike_df <- function(sim){
  z_post <- (sim$model_spike$beta !=0)*1
  z_true <- sim$data$z
  p <- sim$config$p
  z_comp <- matrix(z_true, nrow = nrow(z_post), ncol = length(z_true), byrow = TRUE)
  acc_j <- colMeans(z_post == z_comp)
  tibble(
    index = 1:p,
    z_true = factor(z_true, levels = c(0, 1)),
    acc_j = acc_j,
    seed = sim$config$seed,
    n = sim$config$n,
    p = sim$config$p
  )
}

plot_accuracy_z <- function(df) {
  
  p <- max(df$index)
  plotz <- ggplot(df, aes(x = index, y = acc_j, color = z_true)) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("0" = "black", "1" = "red"), name = "True z") +
    scale_x_continuous(breaks = 1:p) +
    labs(
      title = paste0("n = ", df$n[1], ", p = ",  df$p[1]),
      x = "Component j",
      y = "Value"
    ) +
    theme_minimal(base_size = 14)
  return(plotz)
}


plot_boxplot_accuracy_z <- function(df_long, fig_dir=NULL, spike=FALSE) {
  p <- unique(df_long$p)
  n <- unique(df_long$n)
  
 boxplot <- ggplot(df_long, aes(x = factor(index), y = acc_j, color= z_true)) +
    geom_boxplot(fill = "lightblue", alpha = 0.5) +
    scale_color_manual(values = c("0" = "black", "1" = "red"), name = "True z") +
    labs(
      #title = paste0("n = ", n, ", p = ", p),
      x = "Component j",
      y = "Value" # expression("Accuracy of "~z[j])
    ) +
    #coord_cartesian(ylim = c(0.8, 1.01)) +
   theme_minimal(base_size = 14)
  
 if(is.null(fig_dir)) return(boxplot)
 else {
   if(spike){
     namefile <- file.path(fig_dir, paste0("zspike_acc_box_p", p, "_n", n, ".pdf"))
   } else{
     namefile <- file.path(fig_dir, paste0("z_acc_box_p", p, "_n", n, ".pdf"))
   }
   pdf(namefile, width = 7, height = 5)
   print(boxplot)
   dev.off()
 }
}


plot_one_simulation <- function(sim, figure_dir) {
  
  n <- sim$config$n
  p <- sim$config$p
  
  # Plot boxplot beta*z
  plot_boxplot2 <- plot_beta_z_boxplot(sim, set_legend = "inner")
  namefile <- file.path(figure_dir, paste0("beta_z_box_p", p, "_n", n, ".pdf"))
  pdf(namefile, width = 7, height = 5)
  print(plot_boxplot2)
  dev.off()

  # Plot accuracy z
  # df <- get_accuracy_z_df(sim)
  # pz <- plot_accuracy_z(df)
  # namefile <- file.path(figure_dir, paste0("z_acc_p", p, "_n", n, ".pdf"))
  # pdf(namefile, width = 7, height = 5)
  # print(pz)
  # dev.off()
  
  # Plot posterior distributions using plot.bayes_glmvs
  # tv = list(true_beta=sim$data$beta, true_z =sim$data$z, true_tau=sim$data$tau)
  # plot(sim$model, data_x = sim$data$X, true_values=tv, figure_folder=figure_dir)
  
  # Plot posterior distributions using plot.bayes_glmvs as.mcmc
  plot(sim$model, data_x = sim$data$X, asmcmc=TRUE, figure_folder=figure_dir)

}


plot_error_vs_n_over_p <- function(error_metrics, metric_name, y_label = NULL, figure_dir) {
  
  figure_dir <- paste0(figure_dir, "/error_metrics")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  
  if(is.null(y_label)) y_label <- metric_name
  
  plot <- ggplot(error_metrics, aes(x = n_over_p, y = .data[[metric_name]])) +
    geom_point(alpha = 0.6, size = 2) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = error_metrics$n_over_p) +
    labs(x = "n / p", y = y_label) +
    theme_minimal(base_size = 14)
  pdf(file = file.path(figure_dir, paste0(metric_name,".pdf")))
  print(plot)
  dev.off()
}


# To plot metrics + confidence intervals for multiple simulation ≠ seeds
plot_error_vs_n_over_p_seeds <- function(error_metrics, metric_name, 
                                         y_label = NULL, set_legend=NULL,
                                         figure_dir, plottype="bars") {
  
  if(is.null(y_label)) y_label <- metric_name
  
  model_levels <- c("BayesVS-GLM", "Bayes-GLM(Oracle)", "SpikeSlab", "Horseshoe", "MLE")
  error_metrics$model_type <- factor(error_metrics$model_type, levels = model_levels)
  
  plot_data <- error_metrics
  
  # --- Exclude models based on metric_name ---
  if (metric_name %in% c("r_mse_tau", "r_mae_tau")) {
    plot_data <- dplyr::filter(plot_data, model_type != "MLE")
  } else if (metric_name == "acc_z") {
    plot_data <- dplyr::filter(plot_data, !(model_type %in% c("Bayes-GLM(Oracle)", "MLE")))
  }
  
  if (plottype== "ribbon") {
    plot_data <- dplyr::filter(plot_data, model_type %in% c("Bayes-GLM(Oracle)", "BayesVS-GLM", "MLE"))
  }
  
  # --- Metric-specific computations ---
  if (metric_name == "acc_z") {  
    plot_data <- dplyr::filter(plot_data, model_type %in% c("BayesVS-GLM", "SpikeSlab"))
    
    # Calcolo nel dominio logit
    eps <- 1e-9
    plot_data <- plot_data %>%
      mutate(
        clipped_metric = pmin(pmax(.data[[metric_name]], eps), 1 - eps),
        logit_metric   = qlogis(clipped_metric)
      )
    
    summary_df <- plot_data %>%
      group_by(model_type, n_over_p) %>%
      summarise(
        mean_logit = mean(logit_metric, na.rm = TRUE),
        sd_logit   = sd(logit_metric, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        lower = plogis(mean_logit - sd_logit),
        upper = plogis(mean_logit + sd_logit),
        mean_val = plogis(mean_logit)
      )
    
  } else if (metric_name == 'coverage_beta' | metric_name == 'coverage_beta1') {
    
    # Calcolo nel dominio logit
    eps <- 1e-9
    plot_data <- plot_data %>%
      mutate(
        clipped_metric = pmin(pmax(.data[[metric_name]], eps), 1 - eps),
        logit_metric   = qlogis(clipped_metric)
      )
    
    summary_df <- plot_data %>%
      group_by(model_type, n_over_p) %>%
      summarise(
        mean_logit = mean(logit_metric, na.rm = TRUE),
        sd_logit   = sd(logit_metric, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        lower = plogis(mean_logit - sd_logit),
        upper = plogis(mean_logit + sd_logit),
        mean_val = plogis(mean_logit)
      )
    
  } else {
    
    # Calcola media e sd per ciascun modello e n/p
    summary_df <- plot_data %>%
      group_by(model_type, n_over_p) %>%
      summarise(
        mean_val = mean(.data[[metric_name]], na.rm = TRUE), 
        sd_val   = sd(.data[[metric_name]], na.rm = TRUE), 
        .groups = "drop"
      ) %>%
      mutate(
        lower = mean_val - sd_val,
        upper = mean_val + sd_val
      )
  
  }

  # --- Plot ---
  dodge_width <- 8 #ifelse(metric_name=="acc_z", 8, 12)  
  x_breaks <- sort(unique(plot_data$n_over_p))
  
  plot <- ggplot(summary_df, aes(x = n_over_p, y = mean_val, 
                                 color = model_type, fill = model_type, 
                                 linetype = model_type, 
                                 shape = model_type, size = model_type)) 
  if (plottype == "ribbon") {
    plot <- plot + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
      geom_line(linewidth = 1) +
      geom_point() 
      #geom_line(linewidth = 1, position = position_dodge(width = dodge_width)) 
  } else {
    plot <- plot + 
      geom_point(position = position_dodge(width = dodge_width)) +
      geom_errorbar(aes(ymin = lower, ymax = upper), 
                    width = 6, linewidth = 0.6, alpha = 1,
                    position = position_dodge(width = dodge_width)) 
  }
  plot <- plot +
    scale_color_manual(values = model_colors, drop = TRUE) +
    scale_fill_manual(values = model_colors, drop = TRUE) +
    scale_linetype_manual(values = model_linetypes, drop = TRUE) +
    scale_shape_manual(values = model_shapes, drop = TRUE) +
    scale_size_manual(values = model_size, drop = TRUE) +
    scale_x_continuous(breaks = x_breaks) + # , labels = scales::label_number(accuracy = 0.01)) +
    labs(x = "n / p", y = y_label,
         color = "Model", fill = "Model", 
         linetype ="Model",
         shape = "Model", 
         size = "Model") +
    theme_minimal(base_size = 14)
  
  # Add ylim only for acc_z or positive errors
  if (metric_name == "acc_z") {
    plot <- plot + coord_cartesian(ylim = c(NA, 1.01))
  } else if (grepl("coverage", metric_name) ){
    plot <- plot + coord_cartesian(ylim = c(NA, 1.01))
  } else {
    plot <- plot + 
      geom_hline(yintercept = 0, color = "grey50", linetype = "dashed")
  }
  
  # Optional legend control
  if (!is.null(set_legend)) {
    if (set_legend == "none") {
      plot <- plot + theme(legend.position = "none")
    } else if (set_legend == "inner") {
      
      if (metric_name == "acc_z") {
        plot <- plot + theme(
          legend.position = c(0.98, 0.1),
          legend.justification = c("right", "bottom"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.box.background = element_blank(),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.title=element_text(size=14), 
          legend.text=element_text(size=14)
        )
      } else {
        plot <- plot + theme(
          legend.position = c(0.98, 0.98),
          legend.justification = c("right", "top"),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.box.background = element_blank(),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.title=element_text(size=14), 
          legend.text=element_text(size=14)
        )
      }
    }
  } else {
    plot <- plot + theme(
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.title=element_text(size=14), 
      legend.text=element_text(size=14)
    )
  }
  
  if (metric_name=='acc_z'){
    if (plottype == "ribbon") {
      metric_name <- paste0(metric_name, "_ribbon")
    }
    pdf(file = file.path(figure_dir, paste0(metric_name, ".pdf")), 
        width = 9, height = 5)
  } else {
    if (plottype == "ribbon") {
      metric_name <- paste0(metric_name, "_ribbon")
    }
    pdf(file = file.path(figure_dir, paste0(metric_name, ".pdf")))
  }
  
  #pdf(file = file.path(figure_dir, paste0(metric_name, ".pdf")))
  print(plot)
  dev.off()
}




