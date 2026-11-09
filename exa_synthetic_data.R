source("utils_data.R")
source("bayes_glmvs_oracle_class_s4.R")

#====================================================================
# EXAMPLE 1: Synthetic dataset - Gaussian model
#====================================================================


# GENERATE DATA.  --------------------------------
fam <- gaussian.bayes()
n <- 100
p <- 6           # number of total covariates 
c <- 1            # 2c = number of correlated variables
d <- 1            # number of noisy covariates
k <- p - 2*c - d  # number of important covariates

set.seed(1234)
true_tau <- rinvgamma(1, shape=3, rate=1) 
true_beta <- runif(n = k, min = 1, max = 6) 
data <-  gen_data_corr(n, p, k, c, d, fam, true_beta, true_tau, 
                       phi_noise=0.001, delta_noise=5)

x <- data$X
y <- data$y
true_z <- data$z

data$L <- NULL
data$Xtest <- NULL
data$ytest <- NULL
data$d  <- NULL
data$c  <- NULL
data$k  <- NULL

GaussianGlmvsData <- data
save(GaussianGlmvsData, file = 'experiments/data/GaussianGlmvsData.rda')

# MCMC --------------------------------------

start_time = Sys.time()
set.seed(123)
num_mcmc <- 500
hyperparam <- list(
  alpha = 1,
  paramtau = list(shape=3, rate =1),
  y0 = rnorm(nrow(x), 10, 5),
  a0 = 1e-5,
  method_sir = 'irls2'
)
mod <- bayes_glmvs(x = data$X, y = data$y, family = fam, args = hyperparam)
mod <- fit(mod, num_mcmc = num_mcmc, burn.in = floor(num_mcmc/10), asmcmc =F)

summary(mod)

plot(mod, data_x=x, true_values =list(true_beta=true_beta, true_z=true_z, true_tau=true_tau))

end_time = Sys.time()

cat("ELAPSED TIME = ", end_time - start_time)



#====================================================================
# EXAMPLE 2: Synthetic dataset - Poisson model
#====================================================================


# GENERATE DATA.  --------------------------------

fam <- poisson.bayes()
n <- 1000
p <- 10           # number of total covariates 
c <- 1            # 2c = number of correlated variables
d <- 1            # number of noisy covariates
k <- p - 2*c - d  # number of important covariates

set.seed(1234)
true_tau <- 1
true_beta <- runif(n = k, min = 0.5, max = 1)
data <-  gen_data_corr(n, p, k, c, d, fam, true_beta, true_tau, phi_noise=0.001, delta_noise=5)

x <- data$X
y <- data$y
true_z <- data$z

data$L <- NULL
data$Xtest <- NULL
data$ytest <- NULL
data$d  <- NULL
data$c  <- NULL
data$k  <- NULL
data$tau <- NULL

PoissonGlmvsData <- data
save(PoissonGlmvsData, file = 'experiments/data/PoissonGlmvsData.rda')


# MCMC --------------------------------------

set.seed(123)
num_mcmc <- 100
hyperparam <- list(
  alpha = 1,
  paramtau = 1,
  y0 = rnorm(n, 5, 1),
  a0 = 1e-6,
  method_sir = 'irls2'
)

mod <- bayes_glmvs(x = x, y = y, family = fam, args = hyperparam)
mod <- fit(mod, num_mcmc = num_mcmc, burn.in = floor(num_mcmc/10), asmcmc =F)
summary(mod)

plot(mod, data_x=x, true_values =list(true_beta=true_beta, true_z=true_z, true_tau=true_tau))


mod2 <- bayes_glmvs_oracle(x = x, y = y, family = fam, args = hyperparam, z_true=true_z)
mod2 <- fit(mod2, num_mcmc = num_mcmc, burn.in = floor(num_mcmc/10), asmcmc =F)
summary(mod2)
plot(mod2, data_x = x, true_values =list(true_beta=true_beta, true_tau=true_tau))



