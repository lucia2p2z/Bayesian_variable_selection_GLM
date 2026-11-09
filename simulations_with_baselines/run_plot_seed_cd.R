library(stringr)
library(dplyr)

subdirs <- list.dirs(config$base_path, full.names = TRUE, recursive = TRUE)

# Don't read the directory 'figures' if it exists
subdirs <- subdirs[!grepl("figures", subdirs)]

# Keep only the sub-directories that ends with /seed_X)
seed_dirs <- subdirs[grepl("/seed_\\d+$", subdirs)]
 
# Extract d (ex. "d1", "d3") 
extract_d <- function(path) {
  str_match(path, "d(\\d+)")[,1]  
}
# Extract c (ex. "c1", "c3") 
extract_c <- function(path) {
  str_match(path, "c(\\d+)")[,1]  
}

# Extract c and d (ex. "c1_d10", "c3_d10") 
extract_c_d <- function(path) {
  str_match(path, "c(\\d+)_d(\\d+)")[,1]  
}

# Extract n (ex. "n100", "n500") 
extract_n <- function(path) {
  str_match(path, "n(\\d+)")[,1]  
}

c_d_groups <- split(seed_dirs, sapply(seed_dirs, extract_c_d))

for (c_d_val in names(c_d_groups)) {
  
  subdirs_c_d <- c_d_groups[[c_d_val]]
  
  figure_dir_c_d <- file.path(config$figure_dir, paste0("figures_", c_d_val))
  dir.create(figure_dir_c_d, showWarnings = FALSE, recursive = TRUE)
  
  # Load all simulations
  sims_loaded <- lapply(subdirs_c_d, load_simulation)
  
  # Extract seed from name_directory
  extract_seed <- function(path) {
    as.integer(str_match(path, "seed_(\\d+)$")[,2])
  }
  seed_vals <- sapply(subdirs_c_d, extract_seed)
  
  # Compute metrics and save with the seed value 
  results <- lapply(sims_loaded, compute_metrics_from_sim)
  error_metrics <- bind_rows(
    Map(function(res, seed) {
      mutate(res$metrics, seed = seed)
    }, results, seed_vals)
  )
  
  # beta_df <- bind_rows(
  #   Map(function(res, seed) {
  #     mutate(res$beta, seed = seed)
  #   }, results, seed_vals)
  # )
  # 
  
  # Save everything
  write.csv(error_metrics, file = file.path(config$base_path, paste0("summary_metrics_", c_d_val, ".csv")), row.names = FALSE)
  
  # Plots all simulations with confidence shaded area 
  plot_error_vs_n_over_p_seeds(error_metrics, "r_mse_beta1", y_label = "RelMSE(beta)", set_legend="inner", figure_dir=figure_dir_c_d, plottype="ribbon")
  plot_error_vs_n_over_p_seeds(error_metrics, "r_maxse_beta1", y_label = "RelMaxSE(beta)", set_legend="inner", figure_dir=figure_dir_c_d,  plottype="ribbon")
  plot_error_vs_n_over_p_seeds(error_metrics, "acc_z", y_label = "Accuracy(z)", set_legend=NULL, figure_dir=figure_dir_c_d, plottype="ribbon")
  
  # Plots all simulations with confidence bars 
  plot_error_vs_n_over_p_seeds(error_metrics, "acc_z", y_label = "Accuracy(z)", set_legend=NULL, figure_dir=figure_dir_c_d)
  plot_error_vs_n_over_p_seeds(error_metrics, "r_mse_beta1", y_label = "RelMSE(beta)", set_legend="inner", figure_dir=figure_dir_c_d)
  plot_error_vs_n_over_p_seeds(error_metrics, "r_mae_beta1", y_label = "RelMAE(beta)", set_legend="inner", figure_dir=figure_dir_c_d)
  plot_error_vs_n_over_p_seeds(error_metrics, "r_maxse_beta1", y_label = "RelMaxSE(beta)", set_legend="inner", figure_dir=figure_dir_c_d)
  plot_error_vs_n_over_p_seeds(error_metrics, "r_maxae_beta1", y_label = "RelMaxAE(beta)", set_legend="inner", figure_dir=figure_dir_c_d)
  plot_error_vs_n_over_p_seeds(error_metrics, "coverage_beta1", y_label = "EC(beta)", set_legend=NULL, figure_dir=figure_dir_c_d)
  
  plot_error_vs_n_over_p_seeds(error_metrics, "mse_beta", y_label = "MSE(beta)", set_legend="inner", figure_dir=figure_dir_c_d)
  plot_error_vs_n_over_p_seeds(error_metrics, "maxse_beta", y_label = "MaxSE(beta)", set_legend="inner", figure_dir=figure_dir_c_d)
  plot_error_vs_n_over_p_seeds(error_metrics, "coverage_beta", y_label = "EC(beta)", set_legend=NULL, figure_dir=figure_dir_c_d)
  
  
  # Plots for one simulation: the first one 
  single_folder <- paste0(figure_dir_c_d, "/fig_seed", seed_vals[[1]])
  dir.create(single_folder, showWarnings = FALSE, recursive = TRUE)
  sims_s <- sims_loaded[sapply(sims_loaded, function(s) s$config$seed == seed_vals[[1]])]
  for (s in sims_s){
    plot_one_simulation(s, single_folder)
    
    p <- s$config$p
    n <- s$config$n
    
    # WITH BASELINES
    namefile <- file.path(single_folder, paste0("beta_baselines_z_box_p", p, "_n", n, "_inner.pdf"))
    pdf(namefile, width = 10, height = 6)
    boxz <- plot_beta_z_boxplot_with_baselines(s, set_ylim = NULL, set_legend="inner")
    print(boxz)
    dev.off()
    
    namefile <- file.path(single_folder, paste0("beta_baselines_z0_box_p", p, "_n", n, "_outer.pdf"))
    pdf(namefile, width = 12, height = 5)
    boxz0 <- plot_beta_z_boxplot_with_baselines(s, which_idx="z0", set_ylim = NULL)
    print(boxz0)
    dev.off()
    
    namefile <- file.path(single_folder, paste0("beta_baselines_z0_box_p", p, "_n", n, ".pdf"))
    pdf(namefile, width = 12, height = 5)
    boxz0 <- plot_beta_z_boxplot_with_baselines(s, which_idx="z0", set_ylim = NULL, set_legend="none")
    print(boxz0)
    dev.off()
    
    namefile <- file.path(single_folder, paste0("beta_baselines_z1_box_p", p, "_n", n, ".pdf"))
    pdf(namefile, width = 6, height = 5)
    boxz1 <- plot_beta_z_boxplot_with_baselines(s, which_idx="z1", set_ylim = NULL, set_legend="none")
    print(boxz1)
    dev.off()
    
    namefile <- file.path(single_folder, paste0("beta_baselines_z1_box_p", p, "_n", n, "_outer.pdf"))
    pdf(namefile, width = 8, height = 5)
    boxz1 <- plot_beta_z_boxplot_with_baselines(s, which_idx="z1", set_ylim = NULL)
    print(boxz1)
    dev.off()
  }
  
  # Plot boxplot of accuracy per components of z, for ≠ seed, and same p & n
  n_groups <- split(subdirs_c_d, sapply(subdirs_c_d, extract_n)) 
  for (n in names(n_groups)) {
    # find into sims_loaded the one that corresponds to this n 
    nval <- as.integer(sub("n", "", n))
    sims_n <- sims_loaded[sapply(sims_loaded, function(s) s$config$n == nval)]
    accuracy_df <- bind_rows(lapply(sims_n, get_accuracy_z_df))
    plot_boxplot_accuracy_z(accuracy_df, figure_dir_c_d)
    
    accuracy_spike_df <- bind_rows(lapply(sims_n, get_accuracy_zspike_df))
    plot_boxplot_accuracy_z(accuracy_spike_df, figure_dir_c_d, spike=TRUE)
  }
  
  message("Saved images to:\n  ", figure_dir_c_d)
  
}
