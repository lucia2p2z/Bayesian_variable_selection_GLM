
base <- "irls2_normal_start_zero/gaussian/Nmcmc_1000_alpha_1_a0_1e-04_y0_1_4_n100/"
rest <- "/seed_1255/model.Rds"
folder <- paste0(base, "/tables_for_claudio_gaussian")
dir.create(folder, recursive = TRUE)

values <- c("gamma05", "gamma09", "gamma095", "gamma099")

for (gammavalue in values){
  path <- paste0(base, gammavalue, rest)
  model <- readRDS(path)
  Z <- model@fit$z
  t <- table(Z[, 1], Z[, 2])
  saveRDS(t, file.path(folder, paste0("table_", gammavalue, ".Rds")))
}



base <- "irls2_normal_start_zero/poisson/Nmcmc_1000_alpha_1_a0_0.001_y0_0_1_n100/"
rest <- "/seed_1255/model.Rds"
folder <- paste0(base, "/tables_for_claudio_poisson")
dir.create(folder, recursive = TRUE)

values <- c("gamma05", "gamma09", "gamma095", "gamma099")

for (gammavalue in values){
  path <- paste0(base, gammavalue, rest)
  model <- readRDS(path)
  Z <- model@fit$z
  t <- table(Z[, 1], Z[, 2])
  saveRDS(t, file.path(folder, paste0("table_", gammavalue, ".Rds")))
}


  
