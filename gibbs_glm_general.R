#library(invgamma)
source('utils_beta_sampling.R')

# Gibbs sampler for a GLM with GENERAL LIKELIHOOD 
# without any assumptions X and integrating out the variable c: 
# z|c \sim \prod_j Bernoulli(c_j)
# c|alpha \sim Dirichlet(alpha)
# Args:
#   x: (matrix dim (n x k)) Covariates.
#   y: (vector of length n) A response variable.
#   alpha: (scalar) controlling the prior for the mixture component weights.
#   beta0: (vector of length k) initial values for the regression coefficients.
#   phi0:  (scalar) Initial value for phi.
#   y0: (vector of length n) Parameter of the prior on Beta.
#   a0: (scalar) Parameter of the prior on Beta.
#   shape_phi0: (scalar) Shape parameter for the prior on phi.
#   num_mcmc:  (scalar) Number of MCMC iterations.
#   burn.in:  (scalar) Number of initial MCMC iterations to discard (default is 0).
#   family: family {stats} provide a convenient way to specify the details of 
#             the models used by functions such as glm. 
#             Examples: gaussian(), poisson(), Gamma(),
# Return:
#   results: list of posterior sample of interests. 
glm_gibbs_sir_oracle <- function(x, y, alpha, beta0, y0, a0, z_true,
                                 hyperparam_prior_tau, 
                                 num_mcmc, burn.in=0, family, 
                                 method_sir='irls2', proposal='tstudent'){
  
  
  k <- ncol(x)
  n <- nrow(x)
  beta <- matrix(NA, num_mcmc, k)
  a0_beta_post <- rep(NA, num_mcmc)
  y0_beta_post <- matrix(NA, num_mcmc, n)
  tau <- rep(NA, num_mcmc)
  
  # Initial starting points of the chain
  tau[1] <- 1
  beta[1,] <- beta0 
  a0_beta_post[1] <- a0
  y0_beta_post[1,] <- y0
  
  control <- c()
  if(family$family == "gaussian"){
    shape_phi0 <- hyperparam_prior_tau$shape
    rate_phi0 <- hyperparam_prior_tau$rate
  } else if(family$family == "Gamma"){
    stop("Gamma - Not implemented yet")
  } else if(family$family == "poisson" || family$family == "binomial"){
    control <- 1 
  }
  
  pb <- txtProgressBar(min = 2, max = num_mcmc, style = 3)
  
  for (i in (2:num_mcmc)) {
    
    setTxtProgressBar(pb, i)
    
    if(i==2) print(paste0("\nProposal used = ", proposal))
    
    coeff <- family$a(tau[i-1])*a0*tau[i-1]
    y0_beta_post[i, ] <- (y +  coeff*y0)/(1+coeff)
    a0_beta_post[i] <- (1 +  coeff)/family$a(tau[i-1])
    
    x0 <- as.matrix(x[,z_true==0])
    x1 <- as.matrix(x[,z_true==1])
    
    #browser()
    if(ncol(x0)!=0){
      betastart <-  beta[i-1, z_true==0] # rep(0, ncol(x0)) # beta[i-1, z_true==0]
      sirbeta0 <- my_sir_laplace(num_samples=20, num_resamples=1, x=x0, 
                                 y0=y0, a0=a0, 
                                 tau=tau[i-1], family=family,
                                 beta0=betastart,
                                 method_beta = method_sir, proposal=proposal, deg_of_freedom=10) 
      beta[i, z_true==0] <- sirbeta0$beta_samples
    }
    if(ncol(x1)!=0){
      betastart <- beta[i-1, z_true==1]  # rep(0, ncol(x1))  beta[i-1, z_true==1]
      #browser()
      sirbeta1 <- my_sir_laplace(num_samples=20, num_resamples=1, x=x1, 
                                 y0=y0_beta_post[i, ], a0=a0_beta_post[i], 
                                 tau=tau[i-1], family=family,
                                 beta0=betastart,
                                 method_beta = method_sir, proposal=proposal, deg_of_freedom=10) 
      beta[i, z_true==1] <- sirbeta1$beta_samples
    } 

    # tau
    if(family$family == "gaussian"){ 
      shape_post <- shape_phi0 + n/2
      control$shape <- shape_post
      
      rate_post <- 0.5*sum( drop(y - x%*%(beta[i,]*z_true))^2) + rate_phi0
      control$rate <- rate_post
    }
    tau[i] <- family$sampling.tau(control)
  }
  
  close(pb)
  
  if (burn.in!=0){
    beta <- beta[-(1:burn.in),]
    tau <- tau[-(1:burn.in)]
    a0_beta_post <- a0_beta_post[-(1:burn.in)]
    y0_beta_post <- y0_beta_post[-(1:burn.in),]
  }
  
  result <- list(
    beta = beta,
    tau = tau,
    a0_beta_post = a0_beta_post,
    y0_beta_post = y0_beta_post
  )
}





