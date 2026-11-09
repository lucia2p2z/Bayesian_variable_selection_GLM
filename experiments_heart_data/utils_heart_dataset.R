## Auxiliary functions for 
# the preprocessing, plot and performing experiments 
# on the Heart Disease dataset. 

# Columns description: 
# https://www.kaggle.com/datasets/redwankarimsony/heart-disease-data?phase=FinishSSORegistration&returnUrl=/datasets/redwankarimsony/heart-disease-data/versions/6?resource=download&SSORegistrationToken=CfDJ8NfOyeaVXHJCjrjiquH0Of5aiG9nCJRwFB21zFCvHKojQ1sSJZVjt7ZKT00BskOSvXxv1vnoDvCTyBT_W7_3RUwqi0Af0crXa48cPYvQZyyFRGNZHylQCS8fCGxaX2p3pasNae5OI8BRLsQjdCwZaULNHMNbW1gIWgqQK92tR7ppapmZ-nIxLg55GHsag4a-il_DaPLQwHnT1d_-fzfd2uoNEQwhwXR9tta6IUqN_Wkb8FIF8hdagdP__QZVdHDkXqfZQU_qFVjm8H6wHySTEoHet34ac2oqUBwe0emTB8eFdFSgkoNYHyx2J37r8BIrgh8JT9BKz8B56qaPNEdHrF7m78vsRutu&DisplayName=Lucia%20Filippozzi 
library(corrplot)
library(tidyr)


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



plot_data_distribution <- function(df, folder) {
  
  ### Plot data distribution
  pdf(paste0(folder,'/sex.pdf'))
  barplot(prop.table(table(df$sex)))
  dev.off()
  
  pdf(paste0(folder,'/cp_chestpain.pdf'))
  barplot(prop.table(table(df$cp)))
  dev.off()
  
  pdf(paste0(folder,'/fbs_fastingbloodsugar.pdf'))
  barplot(prop.table(table(df$fbs)))
  dev.off()
  
  pdf(paste0(folder,'/rest_ecg.pdf'))
  barplot(prop.table(table(df$restecg)))
  dev.off()
  
  pdf(paste0(folder,'/thalch_maxheartrate.pdf'))
  hist(df$thalch, freq = FALSE, main="", xlab = "Maximum heart rate achieved")
  dev.off()
  
  pdf(paste0(folder,'/age.pdf'))
  hist(df$age, freq = FALSE, main="", xlab = "Age")
  dev.off()
  
  pdf(paste0(folder,'/cholesterol.pdf'))
  hist(df$chol, freq = FALSE, main="", xlab = "Cholesterol")
  dev.off()
  
  pdf(paste0(folder,'/trestbps_restingbloodpressure.pdf'))
  hist(df$trestbps, freq = FALSE, main="", xlab = "Resting blood pressure")
  dev.off()
  
  pdf(paste0(folder,'/trestbps_vs_chol.pdf'))
  ggplot(df, aes(x=chol, y=trestbps)) +  geom_point(aes(col=target))
  dev.off() 
  
  pdf(paste0(folder,'/target_before_preprocess.pdf'))
  barplot(prop.table(table(df$num)), xlab = 'Heart diseases')
  dev.off()
  
  pdf(paste0(folder,'/target_after_preprocess.pdf'))
  barplot(prop.table(table(df$target)), xlab = 'Presence of Heart disease')
  dev.off()
  
  # Age Distribution Plot
  pdf(paste0(folder,'/age_distribution_by_disease.pdf'))
  p <- ggplot(df, aes(x = factor(num), y = age, fill = factor(num))) +
    geom_boxplot() +
    labs(title = "Age Distribution by Heart Disease Status", 
         x = "Heart Disease (0: No, >= 1: Yes)", 
         y = "Age") +
    theme_minimal()
  print(p)
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
  # We preprocess the dataset by 
  # 1. removing "id" feature and the features where too many entries are missing. 
  #    (slope of the peak ST segment, number of major vessels, and thalassemia background)
  # 2. removing records where at least one feature is missing. 
  # 3.Transform as.factor some variables 
  # 4. Transform the target variable 'num' into binary variable:
  #    0: No heart disease (absence of disease).
  #    1-4: Presence of heart disease, with varying severity 
  # 5. perform one-hot encoding of the two categorical (and non binary) features
  #    (chest pain type and resting electrochardiographic results) 
  # 6. Normalize features per center
  
  
  df <- dataset 
  
  names(df) <- gsub(" ", "_", names(df))
  names(df) <- gsub("-", "_", names(df))
  names(df) <- tolower(names(df))
  
  # 1 - Remove useless columns and columns with NaNs > 30%
  df[df == ""] <- NA
  na_prop <- colMeans(is.na(df))  # proportion of NA per column
  cols_to_remove <- names(na_prop[na_prop > 0.3])
  cat("Columns with % of NaNs > 30%: ", cols_to_remove) # slope, ca,  thal
  df <- df[, !(names(df) %in% cols_to_remove)]
  df$id <- NULL 
  
  # 2 - Drop rows with NaNs  
  df <- na.omit(df)
  # Drop row with anomaly n resting blood pressure = 0
  df <- df[df$trestbps>0,]
  
  # 3 - Transform as factors 
  df$sex <- as.factor(df$sex)
  df$dataset <- as.factor(df$dataset)
  df$cp  <- as.factor(df$cp)
  df$fbs <- as.factor(df$fbs)
  df$restecg <- as.factor(df$restecg)
  df$exang <- as.factor(df$exang)

  # 4 - Encode the target variable
  df$target <- factor(ifelse(df$num >= 1, "Disease", "NoDisease"))
  
  if(!is.null(folder)){
    plot_data_distribution(df, folder)
  }
  
  # 5 -  remove num variable and dataset variable
  df$num <- NULL
  df$dataset <- NULL
  
  df[] <- lapply(df, function(x) {
    if (is.factor(x)) {
      levels(x) <- gsub("[- ]", "_", levels(x))
      x
    } else {
      x
    }
    })
  #df <- df[-nrow(df), ]
  return(df)
}


prepare_dataset_for_BayesGlm <- function(covariates) {
  x_numeric <- model.matrix(~ ., data = covariates)  # -1 removes intercept so all columns are explicit
  
  # Convert to data frame (optional)
  x_numeric <- as.data.frame(x_numeric)
  x_numeric$`(Intercept)` <- NULL
  
  cols <- c("sexMale","cpatypical_angina", "cpnon_anginal", "cptypical_angina",
            "fbsTRUE", "restecgnormal",  "restecgst_t_abnormality", "exangTRUE")
  for(c in cols){
    x_numeric[[c]] <- x_numeric[[c]] + 1
  }
  
  return(x_numeric)
}


# The following function may be removed if useless: 
plot_beta_boxplot_full <- function(model, model2, model_glm, df) {
  
  # Info 
  beta_post <- model@fit$beta
  beta_post2 <- model2@fit$beta
  z_post <- model@fit$z
  beta_glm <- model_glm$coefficients
  

  # Prepare for boxplot
  p <- ncol(df)
  beta_levels <- names(df) #paste0("beta", 1:p)
  colnames(beta_post) <- beta_levels
  colnames(beta_post2) <- beta_levels
  beta_labels_expr <- beta_levels 
  
  
  df_long1 <- pivot_longer(as.data.frame(beta_post), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>%
    mutate(posterior = "BayesGLM")
  
  df_long2 <- pivot_longer(as.data.frame(beta_post2), 
                           cols = everything(), 
                           names_to = "beta", 
                           values_to = "value") %>%
    mutate(posterior = "Bayes Full Model")
  
  df_long_all <- bind_rows(df_long1, df_long2)
  df_long_all$beta <- factor(df_long_all$beta, levels = beta_levels)
  
  
  beta_mean_df <- data.frame(beta = beta_levels, value = apply(beta_post, 2, mean))
  beta_mean_df2 <- data.frame(beta = beta_levels, value = apply(beta_post2, 2, mean))
  beta_glm_df <- data.frame(beta = beta_levels, value = beta_glm)
  
  beta_mean_df$posterior <- "BayesGLM"
  beta_mean_df2$posterior <- "Bayes Full Model"
  
  beta_means_all <- bind_rows(beta_mean_df, beta_mean_df2)
  rownames(beta_means_all) <- NULL
  boxplot <- ggplot(df_long_all, aes(x = beta, y = value , fill=factor(posterior))) +
    geom_boxplot(alpha = 0.5) + #fill = model_colors[["BayesGLM"]], alpha = 0.5) +
    #geom_point(data = beta_means_all, aes(x = beta, y = value, color = posterior), shape = 4, size=2 ) +
    #geom_point(data = beta_glm_df, aes(x = beta, y = value, color = "MLE"), shape = 3, size=2 ) +
    theme_minimal() +
    scale_fill_manual(
      name = "Posterior",
      values = c("BayesGLM" = model_colors[["BayesGLM"]],
                 "Bayes Full Model" = model_colors[["BayesGLM_oracle"]])
    ) +
    # scale_color_manual(
    #   name = "Legend",
    #   values = c("beta*z" = "blue", 
    #              "beta_full" = "darkorange", 
    #              "MLE" = model_colors[["MLE"]])
    # ) +
    scale_x_discrete(labels = beta_labels_expr) +
    # theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(title = bquote("Posterior distribution of coefficients " *beta ), 
         x=" ",
         y = "Value")

  return(boxplot)
  
}
