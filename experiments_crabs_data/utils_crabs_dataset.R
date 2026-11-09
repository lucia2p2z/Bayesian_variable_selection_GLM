## Auxiliary functions for 
# the preprocessing, plot and performing experiments 
# on the Heart Disease dataset. 

# Columns description: 
# https://www.kaggle.com/datasets/redwankarimsony/heart-disease-data?phase=FinishSSORegistration&returnUrl=/datasets/redwankarimsony/heart-disease-data/versions/6?resource=download&SSORegistrationToken=CfDJ8NfOyeaVXHJCjrjiquH0Of5aiG9nCJRwFB21zFCvHKojQ1sSJZVjt7ZKT00BskOSvXxv1vnoDvCTyBT_W7_3RUwqi0Af0crXa48cPYvQZyyFRGNZHylQCS8fCGxaX2p3pasNae5OI8BRLsQjdCwZaULNHMNbW1gIWgqQK92tR7ppapmZ-nIxLg55GHsag4a-il_DaPLQwHnT1d_-fzfd2uoNEQwhwXR9tta6IUqN_Wkb8FIF8hdagdP__QZVdHDkXqfZQU_qFVjm8H6wHySTEoHet34ac2oqUBwe0emTB8eFdFSgkoNYHyx2J37r8BIrgh8JT9BKz8B56qaPNEdHrF7m78vsRutu&DisplayName=Lucia%20Filippozzi 
library(corrplot)
library(tidyr)


plot_data_distribution <- function(df, folder) {
  
  ### Plot data distribution
  
  pdf(paste0(folder,'/satellites.pdf'))
  barplot(prop.table(table(df$Satellites)))
  dev.off()
  
  pdf(paste0(folder,'/dark.pdf'))
  barplot(prop.table(table(df$Dark)))
  dev.off()
  
  pdf(paste0(folder,'/goodspine.pdf'))
  barplot(prop.table(table(df$GoodSpine)))
  dev.off()
  
  pdf(paste0(folder,'/rep1.pdf'))
  hist(df$Rep1, freq = FALSE, main="", xlab = "Rep 1")
  dev.off()
  
  pdf(paste0(folder,'/rep2.pdf'))
  hist(df$Rep2, freq = FALSE, main="", xlab = "Rep 2")
  dev.off()
  
  
  # Correlation Heatmap
  pdf(paste0(folder,'/correlation_numeric_features.pdf'))
  numeric_cols <- df[, sapply(df, is.numeric)]
  cor_matrix <- cor(numeric_cols, use = "complete.obs")
  corrplot(cor_matrix, method = "circle", type = "upper")
  dev.off()
  
  # Numeric Covariates Distribution Plot
  pdf(paste0(folder,'/distribution_numeric_features.pdf'))
  df_long <- numeric_cols[,!colnames(numeric_cols) %in% "num"] %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")
  
  # Create faceted boxplots
  p <- ggplot(df_long, aes(x = "", y = Value)) +
    geom_boxplot(fill = "lightblue", color = "black", outlier.alpha = 0.4) +
    facet_wrap(~ Variable, scales = "free_y") +  # separate scale for each feature
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1.2, "lines")  
    ) +
    labs(title = "Boxplots of Numeric Features",
         x = NULL,
         y = "Value")
  print(p)
  dev.off()
}


preprocess <- function(dataset, folder=NULL){
  
  df <- dataset 
  
  names(df) <- gsub(" ", "_", names(df))
  names(df) <- gsub("-", "_", names(df))
  
  # 4 - Encode the target variable
  df$target <- df$Satellites
  
  if(!is.null(folder)){
    plot_data_distribution(df, folder)
  }
  
  # 5 -  remove satellites variable
  df$Satellites <- NULL

  return(df)
}


prepare_dataset_for_BayesGlm <- function(covariates) {
  x_numeric <- model.matrix(~ ., data = covariates)  # -1 removes intercept so all columns are explicit
  
  # Convert to data frame (optional)
  x_numeric <- as.data.frame(x_numeric)
  x_numeric$`(Intercept)` <- NULL
  
  x_numeric$Darkyes <- x_numeric$Darkyes+1
  x_numeric$GoodSpineyes <- x_numeric$GoodSpineyes+1
  return(x_numeric)
}



