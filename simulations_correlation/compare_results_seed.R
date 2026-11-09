library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

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

base_path <- config$base
subdirs <- list.dirs(base_path, full.names = TRUE, recursive = TRUE)

# Don't read the directory 'figures' if it exists
subdirs <- subdirs[!grepl("figures", subdirs)]

# Keep only the sub-directories that ends with /seed_X)
seed_dirs <- subdirs[grepl("/seed_\\d+$", subdirs)]

# Extract seed (ex. "seed_1") 
extract_seed <- function(path) {
  str_match(path, "seed_(\\d+)")[,1]  
}
#seed_groups <- split(seed_dirs, sapply(seed_dirs, extract_seed))


# Extract d (ex. "d1", "d3") 
extract_gamma <- function(path) {
  str_match(path, "gamma(\\d+)")[,1]  # es. "d3"
}
gamma_groups <- split(seed_dirs, sapply(seed_dirs, extract_gamma))

df_accuracy_tot <- data.frame()

for (g_val in names(gamma_groups)) {
  
  print(g_val)
  
  seed_dirs <- gamma_groups[[g_val]]
  post_prob_z <- data.frame()
  #total_metrics <- data.frame()  # initialize empty
  
  folder_figures <- file.path(base_path, "AAA_figures")
  folder_figures_g <- file.path(folder_figures, g_val)
  dir.create(folder_figures_g, recursive = TRUE, showWarnings = FALSE)
  
  folder_tables <- file.path(base_path, "AAA_tables")
  dir.create(folder_tables, recursive = TRUE, showWarnings = FALSE)
  
  i <- 0 
  all_df <- list()
  accuracy_df <- data.frame()
  for (s in seed_dirs) {
    i <- i + 1
    
    sim <- load_simulation(s)
    
    if(i==1) {
      folder_figures_g_s <- file.path(folder_figures_g, "first_seed")
      dir.create(folder_figures_g_s, recursive = TRUE, showWarnings = FALSE)
      plot_one_simulation(sim, folder_figures_g_s)
    }
    
    plot_one_simulation(sim, s)
    
    mod <- sim$model
    post_prob_z <- rbind(post_prob_z, apply(mod@fit$z, MARGIN=2, FUN=mean))
    
    Z <- mod@fit$z
    ZZ <- apply(Z, 1, function(x) paste(x, collapse = ""))
    tZ <- table(ZZ)
    stZ <- sort(tZ, decreasing = TRUE)
    saveRDS(stZ, file=file.path(s, "table_z.Rds"))
    
    df <- data.frame(
      model = names(tZ),
      freq = as.numeric(tZ)
    )
    
    all_df[[s]] <- df
    
    # Metrics 
    accuracy_df <- bind_rows(accuracy_df, get_accuracy_z_df(sim))

    #metrics <- readRDS(file.path(s, "prediction_metrics.Rds"))
    #metrics_df <- as.data.frame(metrics)
    #seed_num <- as.numeric(str_match(s, "seed_(\\d+)$")[, 2])
    #metrics_df$metric <- rownames(metrics_df)
    #metrics_df$seed <- seed_num
    #total_metrics <- rbind(total_metrics, metrics_df)
  } 
  
  df_accuracy_tot <-  bind_rows(df_accuracy_tot, accuracy_df)
  plot_boxplot_accuracy_z(accuracy_df, folder_figures_g)

  
  combined_df <- do.call(rbind, all_df)
  # sum frequencies of same models across runs
  summary_df <- combined_df %>%
    group_by(model) %>%
    summarise(freq = sum(freq)) %>%
    arrange(desc(freq))
  
  summary_df$prob <- summary_df$freq / sum(summary_df$freq)
  
  saveRDS(summary_df, file = paste0(folder_tables, "/tot_table_z_", g_val, ".Rds"))
  
  
  p <- sim$config$p
  names(post_prob_z) <- 1:p
  prob_to_plot <- post_prob_z %>% pivot_longer(everything(), values_to="value", names_to="variable") 
  # preserve original column order
  prob_to_plot$variable <- factor(prob_to_plot$variable, levels = colnames(post_prob_z))
  
  pdf(file.path(folder_figures_g, paste0('posterior_inclusion_prob_', g_val, '.pdf')))
  plt <- ggplot(prob_to_plot, aes(x=variable, y=value)) + 
    geom_boxplot()+
    labs(x = "Component j", y = "Value", 
         title = paste("Gamma =", sim$config$gamma)) +
    theme_minimal() 
  print(plt)
  dev.off()
  
  
  
  # ============
  
  if(FALSE) {
    rownames(total_metrics) <- NULL
    colnames(post_prob_z) <- names_variables
    
    
    # Reshape to long format with model and value
    long_metrics <- data.frame(
      metric = rep(total_metrics$metric, 3),
      seed   = rep(total_metrics$seed, 3),
      Model  = c(rep("BayesVS-GLM", nrow(total_metrics)),rep("Bayes-GLM", nrow(total_metrics)), rep("MLE", nrow(total_metrics))),
      value  = c(total_metrics$`BayesVS-GLM`,total_metrics$`Bayes-GLM`, total_metrics$MLE)
    )
    
    # Compute mean and sd per metric per model
    agg_metrics <- aggregate(value ~ metric + Model, data = long_metrics,
                             FUN = function(x) c(mean = mean(x), sd = sd(x)))
    agg_metrics$mean <- agg_metrics$value[, "mean"]
    agg_metrics$sd   <- agg_metrics$value[, "sd"]
    agg_metrics$value <- NULL
    
    
    pdf(paste0(base_path,'/error_metrics.pdf'))
    plt <- ggplot(long_metrics[long_metrics$metric != "LL", ], aes(x = metric, y = value, fill = Model)) + 
      geom_boxplot()+ 
      labs(x = "", y = "Value") +
      theme_minimal() +
      scale_fill_manual(values = model_colors)
    print(plt)
    dev.off()
    
    saveRDS(long_metrics, file=paste0(base_path, "/total_metrics.Rds"))
    saveRDS(agg_metrics, file=paste0(base_path, "/aggregated_metrics.Rds"))
    
    #### TABLE FOR LATEX ========================================
    t <- as.data.frame(agg_metrics)
    t$mean <- round(as.numeric(t$mean), 3)
    t$sd <- round(as.numeric(t$sd), 3)
    t$`mean±sd` <- paste0(round(t$mean, 2), " ± ", round(t$sd, 2))
    t$mean <- NULL
    t$sd <- NULL
    
    t_wide <- t %>%
      pivot_wider(
        names_from = Model,
        values_from = `mean±sd`
      )
    
    
  }
  
}


pdf(file.path(folder_figures, "Acc_z_median_by_z.pdf"))
p <- plot_accuracy_vs_gamma(df_accuracy_tot, TRUE, "median")
print(p)
dev.off()

pdf(file.path(folder_figures, "Acc_z_mean_by_z.pdf"))
p <- plot_accuracy_vs_gamma(df_accuracy_tot, TRUE, "mean")
print(p)
dev.off()

pdf(file.path(folder_figures, "Acc_z_median.pdf"))
p <- plot_accuracy_vs_gamma(df_accuracy_tot, FALSE, "median")
print(p)
dev.off()

pdf(file.path(folder_figures, "Acc_z_mean.pdf"))
p <- plot_accuracy_vs_gamma(df_accuracy_tot, FALSE, "mean")
print(p)
dev.off()



