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

extract_seed <- function(path) {
  as.integer(str_match(path, "seed_(\\d+)$")[,2])
}

c_d_groups <- split(seed_dirs, sapply(seed_dirs, extract_c_d))
c_d_val <- names(c_d_groups)[1]

# Extract directory for one value of (c,d) and one seed
subdirs_c_d <- c_d_groups[[c_d_val]]
seed_vals <- sapply(subdirs_c_d, extract_seed)

subdirs_c_d <- subdirs_c_d[sapply(subdirs_c_d, function(s) grepl(seed_vals[[1]], s))]
sims_s <- lapply(subdirs_c_d, load_simulation)


# Plots for one simulation: the first one 
figure_dir_c_d <- file.path(config$figure_dir, paste0("figures_", c_d_val))
dir.create(figure_dir_c_d, showWarnings = FALSE, recursive = TRUE)
single_folder <- paste0(figure_dir_c_d, "/fig_seed", seed_vals[[1]], "/main_figures")
dir.create(single_folder, showWarnings = FALSE, recursive = TRUE)

for (s in sims_s){
  
  p <- s$config$p
  n <- s$config$n
  
  namefile <- file.path(single_folder, paste0("beta_z1_box_p", p, "_n", n, ".pdf"))
  pdf(namefile, width = 5, height = 5)
  boxz1 <- plot_beta_z_boxplot(s, which_idx="z1", set_ylim = c(-0.5, 1.25), set_legend="none")
  print(boxz1)
  dev.off()
  
  namefile <- file.path(single_folder, paste0("beta_z0_box_p", p, "_n", n, "_inner.pdf"))
  pdf(namefile, width = 7, height = 5)
  boxz0 <- plot_beta_z_boxplot(s, which_idx="z0", set_ylim = c(-0.5, 1.25), set_legend="inner")
  print(boxz0)
  dev.off()
  
  namefile <- file.path(single_folder, paste0("beta_z0_box_p", p, "_n", n, ".pdf"))
  pdf(namefile, width = 7, height = 5)
  boxz0 <- plot_beta_z_boxplot(s, which_idx="z0", set_ylim = c(-0.5, 1.25))
  print(boxz0)
  dev.off()
  
  namefile <- file.path(single_folder, paste0("beta_z_box_p", p, "_n", n, "_inner.pdf"))
  pdf(namefile, width = 7, height = 5)
  boxz0 <- plot_beta_z_boxplot(s, set_ylim = c(-0.5, 1.25), set_legend="inner")
  print(boxz0)
  dev.off()
  
  
  # WITH BASELINES
  namefile <- file.path(single_folder, paste0("beta_baselines_z_box_p", p, "_n", n, "_inner.pdf"))
  pdf(namefile, width = 10, height = 6)
  boxz <- plot_beta_z_boxplot_with_baselines(s, set_ylim = c(-0.5, 1.25), set_legend="inner")
  print(boxz)
  dev.off()
  
  namefile <- file.path(single_folder, paste0("beta_baselines_z0_box_p", p, "_n", n, "_outer.pdf"))
  pdf(namefile, width = 12, height = 5)
  boxz0 <- plot_beta_z_boxplot_with_baselines(s, which_idx="z0", set_ylim = c(-0.5, 1.25))
  print(boxz0)
  dev.off()
  
  namefile <- file.path(single_folder, paste0("beta_baselines_z0_box_p", p, "_n", n, ".pdf"))
  pdf(namefile, width = 12, height = 5)
  boxz0 <- plot_beta_z_boxplot_with_baselines(s, which_idx="z0", set_ylim = c(-0.5, 1.25), set_legend="none")
  print(boxz0)
  dev.off()
  
  namefile <- file.path(single_folder, paste0("beta_baselines_z1_box_p", p, "_n", n, ".pdf"))
  pdf(namefile, width = 6, height = 5)
  boxz1 <- plot_beta_z_boxplot_with_baselines(s, which_idx="z1", set_ylim = c(-0.5, 1.25), set_legend="none")
  print(boxz1)
  dev.off()
  
  namefile <- file.path(single_folder, paste0("beta_baselines_z1_box_p", p, "_n", n, "_outer.pdf"))
  pdf(namefile, width = 8, height = 5)
  boxz1 <- plot_beta_z_boxplot_with_baselines(s, which_idx="z1", set_ylim = c(-0.5, 1.25))
  print(boxz1)
  dev.off()
  
  
}



