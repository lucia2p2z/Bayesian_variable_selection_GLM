z_true <- c(1,0)
p <- 2

# palette.colors(palette = "R4")
#"#000000" "#DF536B" "#61D04F" "#2297E6" "#28E2E5" "#CD0BBC" "#F5C710" "#9E9E9E"
model_colors <- c(
  "BayesVS-GLM"       = "#2297E6",
  "Bayes-GLM(Oracle)" = "#28E2E5",
  "Bayes-GLM"         = "#61D04F",
  "MLE"               = "#CD0BBC",  
  "True beta"         = "#DF536B"    
)
model_linetypes <- c(
  "BayesVS-GLM"       = "solid",
  "Bayes-GLM(Oracle)" = "dashed",
  "MLE"               = "dotted"
)

model_shapes <- c(
  "BayesVS-GLM"        = 16,  # solid circle
  "Bayes-GLM(Oracle)"  = 17,  # solid triangle
  "MLE"                = 18   # cross
)

model_size <- c(
  "BayesVS-GLM"        = 1.5,  # solid circle
  "Bayes-GLM(Oracle)"  = 3,  # solid triangle
  "MLE"                = 3    # cross
)


get_accuracy_z_df <- function(sim){
  z_post <- sim$model@fit$z
  z_comp <- matrix(z_true, nrow = nrow(z_post), ncol = length(z_true), byrow = TRUE)
  acc_j <- colMeans(z_post == z_comp)
  tibble(
    index = 1:p,
    z_true = factor(z_true, levels = c(0, 1)),
    acc_j = acc_j,
    seed = sim$config$seed,
    gamma = sim$config$gamma
  )
}


load_simulation <- function(folder_path) {
  list(
    path = folder_path,
    data = readRDS(file.path(folder_path, "data.Rds")),
    config = readRDS(file.path(folder_path, "config_used.Rds")),
    model = readRDS(file.path(folder_path, "model.Rds")),
    #model_oracle = readRDS(file.path(folder_path, "model_oracle.Rds")),
    predictions = readRDS(file.path(folder_path, "predictions.Rds")),
    #predictions_oracle = readRDS(file.path(folder_path, "predictions_oracle.Rds")),
    model_glm = readRDS(file.path(folder_path, "model_glm.Rds")),
    model_z1 = readRDS(file.path(folder_path, "model_z1.Rds")),
    predictions_z1 = readRDS(file.path(folder_path, "predictions_z1.Rds"))
  )
}


# PLOT =====================================================


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
    theme_minimal() +
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
    theme_minimal() +
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

plot_boxplot_accuracy_z <- function(df_long, fig_dir=NULL) {
  p <- unique(df_long$p)
  n <- unique(df_long$n)
  gamma <- unique(df_long$gamma)
  
  boxplot <- ggplot(df_long, aes(x = factor(index), y = acc_j, color= z_true)) +
    geom_boxplot(fill = "lightblue", alpha = 0.5) +
    scale_color_manual(values = c("0" = "black", "1" = "red"), name = "True z") +
    labs(
      title = paste0("gamma = ", gamma),
      x = "Component j",
      y = "Value" # expression("Accuracy of "~z[j])
    ) +
    coord_cartesian(ylim = c(0, 1.01)) +
    theme_minimal() 
  
  if(is.null(fig_dir)) return(boxplot)
  else {
    namefile <- file.path(fig_dir, paste0("z_acc_box_gamma", gamma, ".pdf"))
    pdf(namefile, width = 7, height = 5)
    print(boxplot)
    dev.off()
  }
}

plot_one_simulation <- function(sim, figure_dir) {
  
  gamma <- sim$config$gamma

  # Plot boxplot beta*z
  plot_boxplot2 <- plot_beta_z_boxplot(sim)
  namefile <- file.path(figure_dir, paste0("beta_z_box_gamma", gamma, ".pdf"))
  pdf(namefile, width = 7, height = 5)
  print(plot_boxplot2)
  dev.off()
  
  # Plot accuracy z
  # df <- get_accuracy_z_df(sim)
  # pz <- plot_accuracy_z(df)
  # namefile <- file.path(figure_dir, paste0("z_acc_gamma", gamma, ".pdf"))
  # pdf(namefile, width = 7, height = 5)
  # print(pz)
  # dev.off()
  
  # Plot posterior distributions using plot.bayes_glmvs
  tv = list(true_beta=sim$data$beta, true_z =sim$data$z, true_tau=sim$data$tau)
  plot(sim$model, data_x = sim$data$X, true_values=tv, figure_folder=figure_dir)
  
  # Plot posterior distributions using plot.bayes_glmvs as.mcmc
  plot(sim$model, data_x = sim$data$X, asmcmc=TRUE, figure_folder=figure_dir)
  
}

plot_accuracy_vs_gamma <- function(df_acc, group_by_z, metric) {
  
  # Transform in logit:
  eps <- 1e-9
  df_acc <- df_acc %>%
    mutate(
      clipped_metric = pmin(pmax(acc_j, eps), 1 - eps),
      logit_metric   = qlogis(clipped_metric)
    )
  
  # Group by:
  if(group_by_z) {
    group_data <- df_acc %>% group_by(z_true, gamma)
  } else {
    group_data <- df_acc %>% group_by(gamma)
  }
  
  # Compute summary metrics:
  if(metric=="mean") {
    plot_data <- group_data %>%
      summarise(
        mean_logit = mean(logit_metric, na.rm = TRUE),
        sd_logit   = sd(logit_metric, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        lower = plogis(mean_logit - sd_logit),
        upper = plogis(mean_logit + sd_logit),
        metric_val = plogis(mean_logit)
      )
  } else if(metric == "median") {
    plot_data <- group_data %>% 
      summarise(
        median_logit = median(logit_metric, na.rm = TRUE),
        upper_logit  = quantile(logit_metric, 0.975, na.rm = TRUE),
        lower_logit  = quantile(logit_metric, 0.025, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        lower = plogis(lower_logit),
        upper = plogis(upper_logit),
        metric_val = plogis(median_logit)
      )
  } else { 
    stop('Unknown type of `metric`.')
  }
  
  # Plot: 
  if (group_by_z){
    plt <- ggplot(plot_data, aes(x = factor(gamma), y = metric_val,
                                 color = z_true, fill = z_true,
                                 group = z_true)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0) +
      geom_line(aes(group = z_true), linewidth = 1) 
  } else {
    plt <- ggplot(plot_data, aes(x = factor(gamma), y = metric_val)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, group = 1), alpha = 0.2, linetype = 0) +
      geom_line(aes(group = 1), linewidth = 1) 
  }
  
  plt <- plt +
    geom_point(size = 2) +
    labs(x = "gamma", y = paste0("Accuracy (", metric, ")")) +
    theme_minimal(base_size = 14)
  return(plt)
}

