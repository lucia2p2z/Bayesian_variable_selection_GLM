source("gibbs_glm_general.R")
source("utils_plot.R")

# Heritage from "bayes_glmvs" 
#
setClass("bayes_glmvs_oracle",
         contains="bayes_glmvs",
         slots = c(z_true = "numeric"))


bayes_glmvs_oracle <- function(x, y, family, control, z_true, ...) {
  new("bayes_glmvs_oracle",
      x = x, y = y, family = family,
      control = control, 
      z_true = z_true,
      fit = NULL)
}


setMethod("fit", "bayes_glmvs_oracle", 
          function(object, num_mcmc = 1000, burn.in = 0, beta0=NULL, asmcmc=FALSE, ...) {
            
            if (is.null(beta0)) {
              beta0 <- rep(0, ncol(object@x))
            }
            
            fit_res <- glm_gibbs_sir_oracle(
              x = object@x,
              y = object@y,
              alpha = object@control$alpha,
              beta0 = beta0,
              y0 = object@control$y0,
              a0 = object@control$a0,
              z_true = object@z_true,
              hyperparam_prior_tau = object@control$paramtau,
              num_mcmc = num_mcmc,
              burn.in = burn.in,
              family = object@family,
              method_sir = object@control$method_sir,
              proposal = object@control$proposal
            )
            
            object@fit <- fit_res
            if (asmcmc){
              fit_res_mcmc <- as.mcmc(object, ...)
              object@fit <- fit_res_mcmc
            } 
            validObject(object)
            return(object)
          })


setMethod("summary", "bayes_glmvs_oracle", 
          function(object) {

            cat("Bayesian GLM with Oracle z\n")
            cat("-----------------------------------\n")
            cat("Family:       ", object@family$family, "\n")
            cat("Observations: ", nrow(object@x), "\n")
            cat("Covariates:   ", ncol(object@x), "\n")
            cat("Prior alpha:  ", object@control$alpha, "\n")
            cat("True z:  ", object@z_true, "\n")
            
            
            if (is.null(object@fit)) {
              cat("\nModel not yet fitted.\n")
              return(invisible(NULL))
            } else if(is.mcmc(object@fit)) {
              cat("\nSummary of the fit results - as.mcmc:\n")
              summary(object@fit)
            } else { 
            
              
              beta <- object@fit$beta
              beta_means <- colMeans(beta)
              beta_sds <- apply(beta, 2, sd)
              
              cat("\nMCMC:", nrow(beta), "iterations (after burn-in)\n")
              
              cat("\nPosterior Means and SDs for Beta:\n")
              beta_table <- data.frame(
                Mean = round(beta_means, 3),
                SD = round(beta_sds, 3)
              )
              colnames(beta_table) <- colnames(object@x)
              print(beta_table)
              
            }
          }
          )


setMethod("plot", "bayes_glmvs_oracle", 
          function(x, data_x, asmcmc=FALSE, true_values=NULL, figure_folder=NULL) {
            if (is.null(x@fit)) {
              stop("Model has not been fitted.")
            }
            
            if(is.mcmc(x@fit)) {   # if the fit is saved as.mcmc 
              print('Plotting as.mcmc because it is ')
              plot(x@fit, ask=F)
            } else if(asmcmc) {    # if the fit is NOT as.mcmc but I want it
              print('Plotting as.mcmc because it is not but I want it ')
              message("Fit is not an mcmc object. Converting on the fly...")
              plot(as.mcmc(x), ask=F)
            } else {               # otherwise do my plot
              print('Plotting my own way')
              
              if(x@family$family == 'gaussian'){
                tau <- median(x@fit$tau)
              }else if(x@family$family == 'poisson' || x@family$family == 'binomial'){
                tau <- 1
              }
              
              n <- nrow(data_x)
              
              plot_posterior_betas(x@fit, x@control$y0, x@control$a0, tau, x@family, data_x, true_values[['true_beta']],
                                   figure_folder)
              
              if(x@family$family == 'gaussian'){
                s <- x@control$paramtau$shape
                r <- x@control$paramtau$rate
                plot_posterior_taus(x@fit, s, r, true_values[['true_tau']], figure_folder, n)
              } 
            }
          }
          )


setMethod("predict", "bayes_glmvs_oracle", 
          function(object, newdata, type = c("response", "link"), full_posterior=F){
            
            if (is.null(object@fit)) stop("Model has not been fitted.")
            if (missing(newdata)) stop("You must provide 'newdata' for prediction.")
            type <- match.arg(type)
            
            newdata <- as.matrix(newdata)
            family <- object@family$family
            
            if(full_posterior) {
              z <- t(replicate(nrow(object@fit$beta), object@z_true))
              eta <- newdata %*% t(object@fit$beta * z)
            } else {
              beta_mean <- colMeans(object@fit$beta) * object@z_true
              eta <- newdata %*% beta_mean
            }
            
            if(type=='link') {
              return(eta)
            } else if(type =='response') {
              return(object@family$linkinv(eta))
            } 
            
          }
)


