library(parallel)
library(BayesGlm)
source("utils_data.R")
source("utils_plot_experiments.R")
source("run_simulation_correlation.R")
source("bayes_glmvs_oracle_class_s4.R")

NUMSIM <- 20
vect_seed <- 1230 + 1:NUMSIM

n <- 500
p <- 2
k <- 1
gamma_vect <- c(0)
  #c(0, 0.5, 0.7, 0.9, 0.95, 0.99) #, 0.92, 0.93, 0.04, 0.95, 0.96, 0.97, 0.98, 0.99)

  
for(gamma in gamma_vect) {
  
  set.seed(1)
  true_tau <- rinvgamma(1, shape=3, rate=1) 
  true_beta <- runif(n = k, min = 0.1, max = 2) 
  
  cat("Beta true =", true_beta,"\n\n")
  cat("Simulation with n =",n, " - p =",p,"\n")
  
  f <- function(seed) {
    
    y0_mean <- 0.5
    y0_std <- 0.5
    config <- list(
      name = "binomial",
      family = binomial.bayes(),
      n = n,
      p = p,
      gamma=gamma,
      seed = seed,
      true_beta = true_beta, 
      true_tau = true_tau,
      hyperparam = list(
        alpha = 0.2,
        paramtau = 1,
        y0_mean = y0_mean,
        y0_std = y0_std,
        y0 = rbinom(n, 1, 0.5), 
        a0 = 0.05,
        method_sir='irls2',
        proposal = 'normal'
      ),
      num_mcmc = 1000,
      start_chain = 'zero'
    )
    
    simula(config)
    
  }
  f(1239) 
  #mclapply(X=vect_seed, FUN=f, mc.cores = 20)
}



