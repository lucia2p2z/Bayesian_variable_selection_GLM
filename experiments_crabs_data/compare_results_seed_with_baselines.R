library(stringr)
library(dplyr)
library(ggplot2)

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

setwd("~/Documents/PhD-unitn/Claudio-Inigo/Bayesian-model-selection/code-paper/experiments_crabs_data")
base_path <- "resultsGLM_cv_90_10_20260618"

subdirs <- list.dirs(base_path, full.names = TRUE, recursive = TRUE)

# Don't read the directory 'figures' if it exists
subdirs <- subdirs[!grepl("figures", subdirs)]

# Keep only the sub-directories that ends with /seed_X)
seed_dirs <- subdirs[grepl("/seed_\\d+$", subdirs)]

# Extract seed (ex. "seed_1") 
extract_seed <- function(path) {
  str_match(path, "seed_(\\d+)")[,1]  
}
seed_groups <- split(seed_dirs, sapply(seed_dirs, extract_seed))

names_variables <- readRDS(file.path(seed_dirs[[1]], "names_variables.Rds"))
total_metrics <- data.frame()  # initialize empty
post_prob_z <- data.frame()
post_prob_z_spike <- data.frame()

for (s in seed_dirs) {
  
  mod <- readRDS(file.path(s, "mod.Rds"))
  mod_spike <- readRDS(file.path(s, "mod_spike.Rds"))
  
  metrics <- readRDS(file.path(s, "prediction_metrics.Rds"))
  metrics_baselines <- readRDS(file.path(s, "prediction_metrics_baselines.Rds"))
  metrics_df <- as.data.frame(metrics)
  metrics_df <- cbind(metrics_df, as.data.frame(metrics_baselines))
  
  if ("SpikeSlab-GLM" %in% colnames(metrics_df)) {
    metrics_df <- metrics_df |> rename("SpikeSlab" = "SpikeSlab-GLM")
  }
  metrics_df <- metrics_df |> select(-any_of("Horseshoe+"))
  number_of_models <- ncol(metrics_df)
  
  # Add identifiers
  seed_num <- as.numeric(str_match(s, "seed_(\\d+)$")[, 2])
  metrics_df$metric <- rownames(metrics_df)
  metrics_df$seed <- seed_num
  
  total_metrics <- rbind(total_metrics, metrics_df)
  post_prob_z <- rbind(post_prob_z, apply(mod@fit$z, MARGIN=2, FUN=mean))
  
  post_prob_z_spike <- rbind(post_prob_z_spike, colMeans(mod_spike$beta != 0))
} 

rownames(total_metrics) <- NULL
colnames(post_prob_z) <- tolower(names_variables)
colnames(post_prob_z_spike) <- tolower(names_variables)

# ============


# Reshape to long format with model and value

long_metrics <- data.frame(
  metric = rep(total_metrics$metric, number_of_models),
  seed   = rep(total_metrics$seed, number_of_models),
  Model  = c(rep("BayesVS-GLM", nrow(total_metrics)),
             rep("Bayes-GLM", nrow(total_metrics)), 
             rep("MLE", nrow(total_metrics)),
             rep("SpikeSlab", nrow(total_metrics)),
             rep("Horseshoe", nrow(total_metrics))
             ),
  value  = c(total_metrics$`BayesVS-GLM`,
             total_metrics$`Bayes-GLM`, 
             total_metrics$MLE,
             total_metrics$SpikeSlab,
             total_metrics$Horseshoe
             )
)

posterior_levels <- c("BayesVS-GLM", "Bayes-GLM", "SpikeSlab", "Horseshoe", "MLE")
long_metrics$Model  <- factor(long_metrics$Model,  levels = posterior_levels)

# Compute mean and sd per metric per model
agg_metrics <- aggregate(value ~ metric + Model, data = long_metrics,
                         FUN = function(x) c(mean = mean(x), sd = sd(x)))
agg_metrics$mean <- agg_metrics$value[, "mean"]
agg_metrics$sd   <- agg_metrics$value[, "sd"]
agg_metrics$value <- NULL


pdf(paste0(base_path,'/crabs_error_metrics_with_baselines.pdf'), width=9, height = 5)
p <- ggplot(long_metrics[long_metrics$metric != "LL", ], aes(x = metric, y = value, fill = Model)) + 
  geom_boxplot(aes(color=after_scale(fill)),
               alpha=0.7)+ 
  labs(x = "", y = "Value") +
  theme_minimal(base_size=14) +
  scale_fill_manual(values = model_colors) + 
  theme(
    axis.text = element_text(size = 16),      # tick labels (x and y)
    axis.title = element_text(size = 18),     # axis titles
    legend.text = element_text(size = 16),    # legend item labels
    legend.title = element_text(size = 18)    # legend title
  )
print(p)
dev.off()

pdf(paste0(base_path,'/crabs_posterior_inclusion_prob.pdf'))
prob_to_plot <- post_prob_z %>% pivot_longer(everything(), values_to="value", names_to="variable") 
# preserve original column order
prob_to_plot$variable <- factor(prob_to_plot$variable, levels = colnames(post_prob_z))

p <- ggplot(prob_to_plot, aes(x=variable, y=value)) + 
  geom_boxplot()+
  labs(x = "", y = "Value") +
  theme_minimal(base_size = 14) + 
  #theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  theme(
    axis.text = element_text(size = 16),      # tick labels (x and y)
    axis.text.x = element_text(angle=25, vjust=1, hjust=1),
    axis.title = element_text(size = 18),     # axis titles
    legend.text = element_text(size = 16),    # legend item labels
    legend.title = element_text(size = 18)    # legend title
  )
  #labs(title= "Posterior inclusion probabilities", x = "", y = "Value") +
print(p)
dev.off()

pdf(paste0(base_path,'/posterior_inclusion_prob_spike.pdf'))
prob_to_plot <- post_prob_z_spike %>% pivot_longer(everything(), values_to="value", names_to="variable") 
# preserve original column order
prob_to_plot$variable <- factor(prob_to_plot$variable, levels = colnames(post_prob_z))

p <- ggplot(prob_to_plot, aes(x=variable, y=value)) + 
  geom_boxplot()+
  labs(x = "", y = "Value") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 16),      # tick labels (x and y)
    axis.text.x = element_text(angle=25, vjust=1, hjust=1),
    axis.title = element_text(size = 18),     # axis titles
    legend.text = element_text(size = 16),    # legend item labels
    legend.title = element_text(size = 18)    # legend title
  )
#labs(title= "Posterior inclusion probabilities", x = "", y = "Value") +
print(p)
dev.off()

saveRDS(long_metrics, file=paste0(base_path, "/total_metrics_with_baselines.Rds"))
saveRDS(agg_metrics, file=paste0(base_path, "/aggregated_metrics_with_baselines.Rds"))


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
