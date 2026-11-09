library(parallel)
library(BayesGlm)
library(BoomSpikeSlab)
library(bayesreg)
source("utils_data.R")
source("run_simulation.R")
source("bayes_glmvs_oracle_class_s4.R")

NUMSIM <- 10
vect_seed <- 1234 + 1:NUMSIM

n_vect <- c(100, 500, 1000)
p <- 20 
d_vect <- c(p%/%2)
c_vect <- c(2)

for(n in n_vect) {
  
  for(c in c_vect) {
    
    for(d in d_vect) {
    
      k <- p - 2*c - d
      num_mcmc <- 5000
      set.seed(1)
      true_tau <- rinvgamma(1, shape=3, rate=1) 
      true_beta <- runif(n = k, min = 1, max = 6) 
      
      cat("Beta true =", true_beta,"\n\n")
      cat("Simulation with n =",n, "in [", n_vect, "] - c =",c,"\n")
      
      f <- function(seed) {
        
        y0_mean <- 1
        y0_std <- 4
        config <- list(
          name = "gaussian_cd",
          family = gaussian.bayes(),
          n = n,
          p = p,
          k = k,
          c = c,
          d = d,
          seed = seed,
          true_beta = true_beta, 
          true_tau = true_tau,
          hyperparam = list(
            alpha = 1,
            paramtau = list(shape=3, rate =1),
            y0_mean = y0_mean,
            y0_std = y0_std,
            y0 = rnorm(n, y0_mean, y0_std),
            a0 = 1e-4,
            method_sir = 'irls2',
            proposal = 'normal'
          ),
          num_mcmc = num_mcmc,
          start_chain = 'glm',
          hyperparam_baselines = list(
            priorspike ='lm',
            num_mcmc = num_mcmc,
            num_mcmc_hs = num_mcmc,
            burnin_hs=1e4
          )
        )
        
        simula(config)  # our model
        simula_baselines(config). # additional baselines
        
      }
      mclapply(X=vect_seed, FUN=f, mc.cores = 10)
    }
  }
}
