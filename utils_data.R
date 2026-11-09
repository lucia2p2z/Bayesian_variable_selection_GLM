skewness <- function(data){
  n <- nrow(data)
  mean_data <- mean(data) 
  sd_data <- sqrt(var(data))
  
  skew <- apply(data, 2, function(z) (n * sum((z - mean_data)^3)) / ((n - 1) * (n - 2) * sd_data^3))
  
  return(skew)
}


generate_data <- function(n, k, model, xmin=NULL, xmax=NULL, 
                          tau=NULL, z=NULL, beta=NULL) {
  
  if (!model %in% c("poisson", "gaussian")) {
    stop(paste0("Model '", model, "' not implemented. Choose 'poisson' or 'gaussian'."))
  }
  if (!is.numeric(n) || n <= 0 || !is.numeric(k) || k <= 0) {
    stop("Both 'n' and 'k' must be positive integers.")
  }
  
  
  if (is.null(z)) z <- rbinom(k, 1, 0.5)
  
  
  if (model == 'poisson'){
    if (is.null(xmin)) xmin <- -3
    if (is.null(xmax)) xmax <- 3
    if (is.null(beta)) beta <- runif(n = k, min = 0.5, max = 1)
  } else if (model=='gaussian'){
    if (is.null(xmin)) xmin <- -10
    if (is.null(xmax)) xmax <- 10
    if (is.null(beta)) beta <- runif(n = k, min = 1, max = 6)
  } 
  
  beta <- beta*z 
  
  x <- matrix(runif(n*k, min = xmin, max = xmax), nrow=n)
  eta <- x %*%beta
  
  if (model == 'poisson') {
    mu <- exp(eta)
    y <- rpois(n, mu)
  } else if (model == 'gaussian') {
    if (is.null(tau)) tau <- rinvgamma(1, shape=3, rate=1)
    mu <- eta
    y <- mu + rnorm(n, 0, sqrt(tau))
  } 
  
  return(list(
    x = x,
    y = y,
    beta = beta,
    z = z,
    tau = if (model == "gaussian") tau else 1
  ))
  
}


orth_norm_x <- function(x){
  x <- apply(x, 2, function(z) z-mean(z))
  svd(x)$u
}


generate_y_from_eta <- function(eta, family, tau) {
  
  n <- length(eta)
  namefamily <- family$family
  
  # To prevent numerical overflow:
  #eta <- pmin(pmax(eta, -15), 15)
  mu <- family$linkinv(eta) 
  
  if (namefamily == "gaussian") {
    y <- rnorm(n, mean = mu, sd = sqrt(tau))
  } else if (namefamily == "poisson") {
    #mu <- exp(eta)
    y <- rpois(n, lambda = mu)
  } else if (namefamily == "binomial") {
    # mu <- plogis(eta)
    y <- rbinom(n, size = 1, prob = mu)
  } else {
    stop(paste("Model ", namefamily, "not implemented yet... Family must be one of 'gaussian', 'poisson', or 'binomial'"))
  }
  
  return(y)
}



gen_data_corr <- function(n, p, k, c, d, family, true_beta, true_tau, 
                          phi_noise=0.001, delta_noise=5){
  
  # Assert to catch errors
  if (!is.numeric(n) || n <= 0) stop("`n` must be a positive integer.")
  if (!is.numeric(k) || k <= 0) stop("`k` must be a positive integer.")
  if (!is.numeric(c) || c < 0) stop("`c` must be a non-negative integer.")
  if (!is.numeric(d) || d < 0) stop("`d` must be a non-negative integer.")
  if (p != k + 2*c + d) stop("`p` must be  = k + 2c + d")
  if (is.null(true_tau) || true_tau <= 0) stop("`true_tau` must be a positive number.")
  if (is.null(true_beta) || length(true_beta) != k) {
    stop("Length of `true_beta` must be equal to `k` (number of active variables).")
  } 
  
  # Define full sample size
  n_test <- floor(0.25 * n)
  N <- n + n_test
  
  namefamily <- family$family
  
  cat("Asking to generate", p, "covariates:  \n", 
      k, "independent, \n", 
      2*c, "correlated with the first", c, "\n",
      d, " noisy independent.\n\n")
  
  cat("Total", N, "observations: ", n, "train +", n_test, "test.\n\n")
  
  cat("Model family :", namefamily)
  
  X <- matrix(NA, nrow = N, ncol = p)
  e <- matrix(NA, nrow = N, ncol = p)
  z <- rep(0, p)
  beta <- rep(0, p)
  
  L <- matrix(rnorm(N * k), nrow = N)
  
  # Compute eta and response variable 
  eta <- L %*% true_beta
  r <- generate_y_from_eta(eta, family, true_tau)
  
  # Target covariates
  e[, 1:k] <- matrix(rnorm(N * k, sd = sqrt(phi_noise)), nrow = N)
  X[, 1:k] <- L + e[, 1:k]
  z[1:k]  <- rep(1, k)
  beta[1:k] <- true_beta
  
  # Correlated covariates
  if (c != 0) {
    e[, (k+1):(k+2*c)] <- matrix(rnorm(N * (2*c), sd = sqrt(delta_noise)), nrow = N)
    X[, (k+1):(k+2*c)] <- cbind(L[,1:c],L[,1:c]) + e[, (k+1):(k+2*c)]
  }
  
  # Noisy independent covariates
  if (d != 0) {
    e[, (k+2*c+1):(k+2*c+d)] <- matrix(rnorm(N * d), nrow = N)
    X[, (k+2*c+1):(k+2*c+d)] <- e[, (k+2*c+1):(k+2*c+d)]
  }
  
  # Split train/test
  Xtest <- X[(n+1):N, ]
  ytest <- r[(n+1):N]
  
  X <- X[1:n, ]
  y <- r[1:n]
  
  return(list(
    X = X,
    Xtest = Xtest,
    L = L,
    beta = beta,
    z = z,
    y = y,
    ytest = ytest,
    tau = true_tau,  
    k = k,
    c = c,
    d = d
  ))
}




