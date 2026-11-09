library(corrplot)

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


plot_data_distribution <- function(dataset, folder){
  
  pdf(paste0(folder,'/distribution_data.pdf'))
  p <- dataset %>%  
    pivot_longer(everything(), values_to="Value", names_to="Variable") %>% 
    ggplot() +
    geom_boxplot(aes(x=Variable, y=Value)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  
  print(p)
  dev.off()
  
  
  cor_matrix <- cor(dataset, use = "complete.obs")
  pdf(paste0(folder,'/corr.pdf'))
  corrplot(cor_matrix, method = "square", type = "upper", diag = FALSE,
           addCoef.col = 'black', number.cex = 0.5)
  dev.off()
}


preprocess <- function(dataset, folder=NULL){
  
  df <- dataset 
  
  names(df) <- gsub(" ", "_", names(df))
  names(df) <- gsub("-", "_", names(df))
  
  # 1 - Encode the target variable
  df$target <- df$MORT
  
  # 2 - remove the old variable
  df$MORT <- NULL
  
  # 3 - standardize  
  df <- as.data.frame(scale(df))
  
  # 4- remove clear outliers of HC and NOX
  df <- df[!df$NOX > 6,] 
  df <- df[!df$HC > 6,] 
  
  if(!is.null(folder)){
    plot_data_distribution(df, folder)
  }
  
  # Remove extremeluy correlateed variables
  df$HC <- NULL
  
  return(df)
}



prepare_dataset_for_BayesGlm <- function(covariates) {
  x_numeric <- model.matrix(~ ., data = covariates)  # -1 removes intercept so all columns are explicit
  
  # Convert to data frame (optional)
  x_numeric <- as.data.frame(x_numeric)
  x_numeric$`(Intercept)` <- NULL
  
  return(x_numeric)
}


