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
  "True beta"         = "#DF536B"    
)
model_linetypes <- c(
  "BayesVS-GLM"       = "solid",
  "Bayes-GLM(Oracle)" = "dashed",
  "MLE"               = "dotted"
)

base_path <- "resultsGLM_cv_90_10_20251104"
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

names_variables <- readRDS(file.path(base_path, "names_variables.Rds"))
total_metrics <- data.frame()  # initialize empty
post_prob_z <- data.frame()

for (s in seed_dirs) {
  
  metrics <- readRDS(file.path(s, "prediction_metrics.Rds"))
  mod <- readRDS(file.path(s, "mod.Rds"))
  metrics_df <- as.data.frame(metrics)
  
  # Add identifiers
  seed_num <- as.numeric(str_match(s, "seed_(\\d+)$")[, 2])
  metrics_df$metric <- rownames(metrics_df)
  metrics_df$seed <- seed_num
  
  total_metrics <- rbind(total_metrics, metrics_df)
  post_prob_z <- rbind(post_prob_z, apply(mod@fit$z, MARGIN=2, FUN=mean))
} 

rownames(total_metrics) <- NULL
colnames(post_prob_z) <- tolower(names_variables)

# ============


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
p <- ggplot(long_metrics[long_metrics$metric != "LL", ], aes(x = metric, y = value, fill = Model)) + 
  geom_boxplot()+ 
  labs(x = "", y = "Value") +
  theme_minimal() +
  scale_fill_manual(values = model_colors)
print(p)
dev.off()

pdf(paste0(base_path,'/posterior_inclusion_prob.pdf'))
prob_to_plot <- post_prob_z %>% pivot_longer(everything(), values_to="value", names_to="variable") 
# preserve original column order
prob_to_plot$variable <- factor(prob_to_plot$variable, levels = colnames(post_prob_z))

p <- ggplot(prob_to_plot, aes(x=variable, y=value)) + 
  geom_boxplot()+
  labs(x = "", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) 
  #labs(title= "Posterior inclusion probabilities", x = "", y = "Value") +
print(p)
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
