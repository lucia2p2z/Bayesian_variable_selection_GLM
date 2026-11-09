library(ggplot2)
library(dplyr)
library(pracma)
library(tidyr)
source('utils_beta_sampling.R')


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
  "BayesVS-GLM"        = 1.5,  # solid circle
  "Bayes-GLM(Oracle)"  = 3,    # solid triangle
  "MLE"                = 3,    # cross
  "SpikeSlab"          = 18,   # cross
  "Horseshoe"          = 3,
  "Horseshoe+"         = 3
)

  


plot_beta_z_boxplot_full <- function(model, model2, model_glm, names_variables=NULL) {
  
  # Info 
  beta_post <- model@fit$beta
  beta_post2 <- model2@fit$beta
  z_post <- model@fit$z
  beta_glm <- model_glm$coefficients
  
  beta_per_z <- beta_post*z_post
  
  # Prepare for boxplot
  p <- ncol(beta_post)
  beta_levels <- if (is.null(names_variables)) paste0("beta", 1:p) else tolower(names_variables) #
  colnames(beta_per_z) <- beta_levels
  colnames(beta_post2) <- beta_levels
  beta_labels_expr <- beta_levels 
  
  
  df_long1 <- pivot_longer(as.data.frame(beta_per_z), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "BayesVS-GLM")
  
  df_long2 <- pivot_longer(as.data.frame(beta_post2), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value")%>% mutate(posterior = "Bayes-GLM")
  
  df_long_all <- bind_rows(df_long1, df_long2)
  df_long_all$beta <- factor(df_long_all$beta, levels = beta_levels)
  
  beta_glm_df <- data.frame(
    beta = rep(1:p, 2),
    value = rep(beta_glm, 2),
    posterior = c(rep("BayesVS-GLM", p), rep("Bayes-GLM", p))
  )  

  
  posterior_levels <- c("BayesVS-GLM", "Bayes-GLM")
  df_long_all$posterior  <- factor(df_long_all$posterior,  levels = posterior_levels)
  beta_glm_df$posterior  <- factor(beta_glm_df$posterior,  levels = posterior_levels)
  
  
  boxplot <- ggplot(df_long_all, aes(x = beta, y = value , fill=posterior)) +
    geom_boxplot(aes(color=after_scale(fill)),
                 position = position_dodge(width = 0.8), 
                 width = 0.6, linewidth = 0.3,
                 alpha=0.7,
                 #outlier.colour = "#000000", 
                 outlier.size = 1, outlier.alpha = 0.1) +  
    geom_point(data = beta_glm_df, 
               aes(x = beta, y = value, group = posterior, shape = "MLE"), 
               position = position_dodge(width = 0.8), color = model_colors[["MLE"]], size=2) +
    theme_minimal() +
    scale_fill_manual(
      name = "Posterior",
      breaks = posterior_levels,
      values = c(model_colors[posterior_levels])
    ) +
    scale_shape_manual(
      name   = "Reference",
      breaks = c("MLE"),
      values = c("MLE" = 18)
    ) +
    guides(
      fill  = guide_legend(order = 1, override.aes = list(shape = NA))#,
      #color = guide_legend(order = 1, override.aes = list(shape = c(NA, NA, NA, NA, 18, 8)))
    ) +
    # scale_fill_manual(
    #   name = "Posterior",
    #   values = c("BayesVS-GLM" = model_colors[["BayesVS-GLM"]],
    #              "Bayes-GLM" = model_colors[["Bayes-GLM"]]),
    #   guide = guide_legend(override.aes = list(shape = NA))
    # ) +
    # scale_color_manual(
    #   name = "Reference",
    #   values = c("MLE" = model_colors[["MLE"]])
    # ) +
    scale_x_discrete(labels = beta_labels_expr) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(x=" ",
         y = "Value")
  
  return(boxplot)
  
}

  

plot_beta_z_boxplot <- function(model, model_glm, df) {
  
  # Info 
  beta_post <- model@fit$beta
  z_post <- model@fit$z
  beta_glm <- model_glm$coefficients
  
  beta_per_z <- beta_post*z_post
  
  # Prepare for boxplot
  p <- ncol(df)
  beta_levels <- tolower(names(df)) #paste0("beta", 1:p)
  colnames(beta_per_z) <- beta_levels
  beta_labels_expr <- beta_levels #lapply(1:p, function(j) bquote(beta[.(j)] * z[.(j)]))  # for axis
  
  
  df_long <- pivot_longer(as.data.frame(beta_per_z), cols = everything(), 
                          names_to="beta", values_to="value")
  df_long$beta <- factor(df_long$beta, levels = beta_levels)
  
  beta_mean_df <- data.frame(beta = beta_levels, value = apply(beta_per_z, 2, mean))
  beta_glm_df <- data.frame(beta = beta_levels, value = beta_glm)
  
  boxplot <- ggplot(df_long, aes(x = beta, y = value)) +
    geom_boxplot(fill = model_colors[["BayesVS-GLM"]], alpha = 0.5) +
    geom_point(data = beta_mean_df, aes(x = beta, y = value, color = "Post.mean"), shape = 4, size=2 ) +
    geom_point(data = beta_glm_df, aes(x = beta, y = value, color = "MLE"), shape = 3, size=2 ) +
    theme_minimal() +
    scale_color_manual(
      name = "Legend",
      values = c("Post.beta" = model_colors[["BayesVS-GLM"]], "Post.mean" ="blue", "MLE" =model_colors[["MLE"]])
    ) +
    scale_x_discrete(labels = beta_labels_expr) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    # labs(title = bquote("Posterior distribution of coefficients " *beta * z), 
    #      x=" ",
    #      y = "Value")
    labs(x=" ", y = "Value")
  
  
  return(boxplot)
  
}



plot_beta_boxplot_only_baselines <- function(mod_spike, mod_hs, names_variables=NULL) {

  # Info
  beta_post1 <- mod_spike$beta
  beta_post2 <- t(mod_hs$beta)
  #beta_post3 <- t(mod_hsplus$beta)


  # Prepare for boxplot
  p <- ncol(beta_post1)
  beta_levels <- if (is.null(names_variables)) paste0("beta", 1:p) else tolower(names_variables) #
  colnames(beta_post1) <- beta_levels
  colnames(beta_post2) <- beta_levels
  # colnames(beta_post3) <- beta_levels
  beta_labels_expr <- beta_levels

  df_long1 <- pivot_longer(as.data.frame(beta_post1),
                           cols = everything(),
                           names_to = "beta",
                           values_to = "value") %>% mutate(posterior = "SpikeSlab")

  df_long2 <- pivot_longer(as.data.frame(beta_post2),
                           cols = everything(),
                           names_to = "beta",
                           values_to = "value")%>% mutate(posterior = "Horseshoe")

  # df_long3 <- pivot_longer(as.data.frame(beta_post3),
  #                          cols = everything(),
  #                          names_to = "beta",
  #                          values_to = "value")%>% mutate(posterior = "Horseshoe+")

  df_long_all <- bind_rows(df_long1, df_long2) 
  df_long_all$beta <- factor(df_long_all$beta, levels = beta_levels)

  posterior_levels <- c("SpikeSlab", "Horseshoe")
  df_long_all$posterior  <- factor(df_long_all$posterior,  levels = posterior_levels)
  
  boxplot <- ggplot(df_long_all, aes(x = beta, y = value , fill=posterior)) +
    geom_boxplot(aes(color=after_scale(fill)),
                 position = position_dodge(width = 0.8),
                 width = 0.6, linewidth = 0.3, alpha = 0.7,
                 #outlier.colour = "#000000", 
                 outlier.size = 1, outlier.alpha = 0.1) +
    theme_minimal() +
    scale_fill_manual(
      name = "Posterior",
      values = c("SpikeSlab" = model_colors[["SpikeSlab"]],
                 "Horseshoe" = model_colors[["Horseshoe"]]),
      guide = guide_legend(override.aes = list(shape = NA))
    ) +
    scale_x_discrete(labels = beta_labels_expr) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(x=" ",
         y = "Value")

  return(boxplot)

}



plot_beta_z_boxplot_with_baselines <- function(model, model_bayes, model_glm, 
                                               mod_spike, mod_hs, names_variables) {
  
  # Info 
  beta_post <- model@fit$beta
  beta_post2 <- model_bayes@fit$beta
  z_post <- model@fit$z
  beta_glm <- model_glm$coefficients
  beta_per_z <- beta_post*z_post
  
  beta_post_baseline1 <- mod_spike$beta
  beta_post_baseline2 <- t(mod_hs$beta)
  
  num_models <- 4
  
  # Prepare for boxplot
  p <- ncol(mod@x)  
  beta_levels <- tolower(names_variables)
  colnames(beta_per_z) <- beta_levels
  colnames(beta_post2) <- beta_levels
  colnames(beta_post_baseline1) <- beta_levels
  colnames(beta_post_baseline2) <- beta_levels
  beta_labels_expr <- beta_levels 
  
  df_long1 <- pivot_longer(as.data.frame(beta_per_z), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "BayesVS-GLM")
  
  df_long2 <- pivot_longer(as.data.frame(beta_post2), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value")%>% mutate(posterior = "Bayes-GLM")
  
  
  df_long3 <- pivot_longer(as.data.frame(beta_post_baseline1), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>% mutate(posterior = "SpikeSlab")
  
  df_long4 <- pivot_longer(as.data.frame(beta_post_baseline2), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value")%>% mutate(posterior = "Horseshoe")
  
  # df_long5 <- pivot_longer(as.data.frame(beta_post_baseline3), 
  #                          cols = everything(), 
  #                          names_to = "beta", 
  #                          values_to = "value")%>% mutate(posterior = "Horseshoe+")
   
  df_long_all <- bind_rows(df_long1, df_long2, df_long3, df_long4) 
  df_long_all$beta <- factor(df_long_all$beta, levels = beta_levels)
  
  beta_glm_df <- data.frame(
    beta = rep(1:p, num_models),
    value = rep(beta_glm, num_models),
    posterior = c(rep("BayesVS-GLM", p), 
                  rep("Bayes-GLM", p), 
                  rep("SpikeSlab", p),
                  rep("Horseshoe", p))
  ) 
  
  posterior_levels <- c("BayesVS-GLM", "Bayes-GLM", "SpikeSlab", "Horseshoe")
  df_long_all$posterior  <- factor(df_long_all$posterior,  levels = posterior_levels)
  beta_glm_df$posterior  <- factor(beta_glm_df$posterior,  levels = posterior_levels)
  
  
  boxplot <- ggplot(df_long_all, aes(x = beta, y = value , fill=posterior)) +
    geom_boxplot(aes(color=after_scale(fill)),
                 position = position_dodge(width = 0.8), 
                 width = 0.6, linewidth = 0.3,
                 alpha=0.7,
                 #outlier.colour = "#000000", 
                 outlier.size = 1, outlier.alpha = 0.1) +  
    geom_point(data = beta_glm_df, aes(x = beta, y = value, group = posterior, color = "MLE"), 
               position = position_dodge(width = 0.8), shape = 18, size=2) +
    theme_minimal(base_size = 14) +
    scale_fill_manual(
      name = "Posterior",
      values = c("BayesVS-GLM" = model_colors[["BayesVS-GLM"]],
                 "Bayes-GLM" = model_colors[["Bayes-GLM"]],
                 "SpikeSlab" = model_colors[["SpikeSlab"]],
                 "Horseshoe" = model_colors[["Horseshoe"]]
                 ),
      guide = guide_legend(override.aes = list(shape = NA))
    ) +
    scale_color_manual(
      name = "Reference",
      values = c("MLE" = model_colors[["MLE"]])
    ) +
    scale_x_discrete(labels = beta_labels_expr) +
    #theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    theme(
      #axis.text = element_text(size = 16),      # tick labels (x and y)
      axis.text.x = element_text(angle=55, vjust=1, hjust=1),
      #axis.title = element_text(size = 18),     # axis titles
      #legend.text = element_text(size = 16),    # legend item labels
      #legend.title = element_text(size = 18)    # legend title
    ) + 
    labs(x=" ",
         y = "Value")
  
  return(boxplot)
  
}

