### FUNCTIONS TO SAMPLE FROM TARGET DISTIRBUTION D(beta; y0, a0)
#
library(mvtnorm)

# KERNEL of the DISTRIBUTION D(beta: y0, a0)
# D(beta; y0, a0) \propto exp{ a_0 [ y_0 X beta - J^\top exp(X beta)] }
D <- function(y0, a0, tau, x, beta, family, logscale=TRUE){
  
  beta <- drop(beta) 
  
  d <- a0*tau*( t(y0)%*%(x%*%beta) - sum(family$b(x%*%beta)) )   
  if(logscale){
    return(d)
  } else {
    return(exp(d))
  }
}


# Compute the hessian of the log-likelihood of beta 
# and evaluate it in beta0.
compute_hessian <- function(beta0, a0, tau, x, family){
  
  a0*tau*t(x) %*% diag(drop(family$bd2(x%*%beta0))) %*% x   
}

# SAMPLING IMPORATNCE RESAMPLING 
# proposal distribution used: Normal((1,1), diag(4,4))
my_sir <- function(num_samples, x, y0, a0, tau, family){
  
  # sampling:
  # TODO ERRORE QUA IN SAMPLES 
  p <- ncol(x)
  mu <- rep(1,p)
  sigma <- rep(4,p)
  s <- rnorm(num_samples, mean=mu, sd = sigma)
  samples <- matrix(s, nrow=num_samples, ncol=p, byrow=TRUE)
  
  # importance:
  logweights <- sapply(1:nrow(samples), function(i) {
    beta <- samples[i, ]  # Extract each row of 'samples' as a 2-dimensional vector
    D(y0=y0, a0=a0, tau=tau, x=x, beta=beta, family=family, logscale=TRUE) - dmvnorm(beta, mean=mu, sigma = diag(sigma, p, p), log=TRUE)
  })
  
  logsumexp <- max(logweights) + log(sum(exp(logweights - max(logweights))))
  weights <- exp(logweights - logsumexp)

  # Resampling: 
  indices <- sample(nrow(samples), size=num_samples, replace=TRUE, prob=weights)
  newsamples <- samples[indices, ]
  
  sir <- list(
    samples=samples,
    weights=weights,
    beta_samples=newsamples
  )
  return(sir)
}


# METROPOLIS-HASTINGS METHOD
# Proposal Normal(betaold, 0.1)
my_mh <- function(num_samples, x, y0, a0, binit, tau, family){
  
  betaold <- binit
  
  i <- 1 
  j <- 0 
  p <- ncol(x)
  beta_samples <- matrix(NA, nrow=num_samples, ncol=p)
  
  num_acc_rej <- list()
  beta_samples[1,] <- binit
  
  sigma <- diag(0.01, p, p)
  
  tot_iter <- 0 
  
  while (i < num_samples) {
    
    tot_iter <- tot_iter + 1 
    betanew <- rmvnorm(1, mean=betaold, sigma=sigma) 
    
    logp <- D(y0=y0, a0=a0, tau=tau, x=x, beta=betanew, family = family, logscale=TRUE) - D(y0=y0, a0=a0, tau = tau, x=x, beta=betaold, family = family, logscale=TRUE)
    logu <- log(runif(1))
    
    if(logu < logp){
      i <- i +1 
      accept <- TRUE
      beta_samples[i,] <- betanew
      betaold <- betanew 
    } 
  }
  
  acc_rate <- num_samples/tot_iter
  mh <-list(
    beta_samples=beta_samples,
    acc_rate=acc_rate
  )
  return(mh)
}


# LAPLACE APPROXIMATION WITHIN SIR
# Args: 
#   proposal ='tstudent' OR 'normal' Proposal distribution to be chosen between: 
#       Both with  Mean = beta_hat = MLE Cov= -Hessian(beta_hat)^-1
#   method = 'optim' OR 'IRLS' to decide which method we want to use to find beta_hat
my_sir_laplace <- function(num_samples, num_resamples, x, y0, a0, tau, family,
                           method_beta, beta0=NULL, proposal='tstudent', 
                           deg_of_freedom=3,
                           tol_laplace=10e-6, maxiter_laplace=1000){
  
  n <- nrow(x)
  p <- ncol(x)
  J <- matrix(1, n, 1)
  
  if(method_beta =='optim') {
    opt <- optim(par=rep(0,p),\(beta) -D(y0=y0, a0=a0, tau=tau, x=x, beta=beta, family=family, logscale=TRUE), hessian = TRUE)
    beta_hat <- opt$par
  } else if(method_beta =='glm'){
    # negative values of the response variable y0 are not allowed when fitting a Poisson GLM:
    # y0[y0 < 0] <- 0
    mod <- glm(y0~x -1, family = family)
    beta_hat <- mod$coefficients
  }else if(method_beta =='irls') {
    if(is.null(beta0)) beta0 <- 0
    res <- my_irls(beta0=beta0, x, y0, tol=1e-6, maxiter = 1000)
    beta_hat <- res$beta_hat
  }else if(method_beta =='irls2'){
    if(is.null(beta0)) beta0 <- rep(0, p)
    res <- my_irls2(beta0, x, y0, family=family, tol=1e-6, maxiter = 1000)
    beta_hat <- res$beta_hat
  } else{
    stop(paste("Unknown method: '", method_beta, "'. Valid methods are: 'optim', 'glm', 'IRLS'.", sep=""))
  }
  
  
  hess <- compute_hessian(beta_hat, a0 = a0, tau=tau, x = x, family = family)
  I <- diag(1, p)
  cov <- solve(hess + 1e-4*I)  # diag(nrow=p) 
  
  # sampling: 
  beta_hat <- drop(beta_hat)
  
  if(proposal=='tstudent'){
    
    if (is.na(deg_of_freedom)) {
      deg_of_freedom <- 3
      warning(paste("deg_of_freedom is NA but required for 'tstudent' proposal. It has been set to:", deg_of_freedom))
    } else if(deg_of_freedom < 3) {
      stop(paste("deg_of_freedom must be greater than 2. You provided: ", deg_of_freedom))
    } 
      
    sigma <-  cov*(deg_of_freedom-2)/deg_of_freedom
    # alternative using Multivariate T, knowing that cov = sigma*v/(v-2)
    samples <- rmvt(num_samples, 
                    df = deg_of_freedom,
                    sigma = sigma, 
                    type = 'shifted', 
                    delta = beta_hat)
    
    # importance:
    logweights <- sapply(1:nrow(samples), function(i) {
      beta <- samples[i, ]  # Extract each row of 'samples' as a 2-dimensional vector
      D(y0=y0, a0=a0, tau=tau, x=x, beta=beta,family=family, logscale=TRUE) - dmvt(beta, df = deg_of_freedom, sigma = sigma, type = 'shifted', delta = beta_hat, log = TRUE)
    })
  } else if(proposal=='normal') {
    samples <- rmvnorm(num_samples, 
                       mean = beta_hat,
                       sigma = cov)
    # importance:
    logweights <- sapply(1:nrow(samples), function(i) {
      beta <- samples[i, ]  # Extract each row of 'samples' as a 2-dimensional vector
      D(y0=y0, a0=a0, tau = tau, x=x, beta, family= family, logscale=TRUE) - dmvnorm(beta, mean = beta_hat, sigma = cov, log = TRUE)
    })
    
  } else{
    stop(paste("Unknown proposal: '", proposal, "'. Valid proposal distributions are: 'tstudent', 'normal'.", sep=""))
  }
  
  logsumexp <- max(logweights) + log(sum(exp(logweights - max(logweights))))
  weights <- exp(logweights - logsumexp)
  

  # Resampling: 
  indices <- sample(nrow(samples), size=num_resamples, replace=TRUE, prob=weights)
  newsamples <- samples[indices, ]
  
  sir <- list(
    samples=samples,
    weights=weights,
    beta_samples=newsamples
  )
  return(sir)
}


# ITERATIVE REWEIGHTED LEAST SQUARE only for 1dim data
my_irls <- function(beta0, x, y0, tol=1e-6, maxiter = 1000) {
  
  diff = Inf
  betaold = drop(beta0)
  i <- 0 
  betas <- list()
  diffs <- list()
  
  while (diff> tol & i < maxiter) {
    i <- i+1
    c <- drop(exp(x %*% betaold))*(x)  # appears several times in the computations 
    
    W <- c*x 
    H <- y0*x + W %*% betaold - c
    
    betanew <- drop(solve(t(W) %*% W)%*% t(W) %*% H)
    
    diff <- abs(betanew - betaold)
    betas[[i]] <- betanew
    diffs[[i]] <- diff
    
    betaold <- betanew 
  }
  res <- list(
    beta_hat = betanew,
    betas = betas,
    diffs = diffs
  )
  return(res)
}

# ITERATIVE REWEIGHTED LEAST SQUARE for multidimensional data (or even scalar)
my_irls2 <- function(beta0, x, y0, family, tol=1e-6, maxiter = 1000) {
  
  p <- ncol(x)
  diff <- Inf
  betaold <- drop(beta0)
  i <- 0 
  betas <- list()
  diffs <- list()
  
  while (diff > tol & i < maxiter) {
    #cat('diff > tol', diff > tol, "\n")
    #cat('i < maxiter', i < maxiter, "\n\n")
    i <- i+1
    #cat('irls i =', i, "\n ")
    
    Deriv1 <- family$bd1(x%*%betaold)
    Deriv2 <- diag(drop(family$bd2(x%*%betaold)))
    res <- y0 - Deriv1
    
    #cat("Det(Deriv2)=", round(det(Deriv2)), "\n")
    
    I <- diag(1, p)
    #if(det(solve(t(x) %*% Deriv2 %*% x + 1e-2*I ) ) < 1e-15) {
      #cat('Determinant = 0 : ', round(solve(det(t(x) %*% Deriv2 %*% x ))))
      #browser()
    #}
    betanew <- drop(solve(t(x) %*% Deriv2 %*% x + 1e-4*I) %*% t(x) %*% res + betaold)

    diff <- ifelse(p==1, abs(betanew - betaold), Norm(betanew - betaold) )
    betas[[i]] <- betanew
    diffs[[i]] <- diff
    
    betaold <- betanew 
  }
  res <- list(
    beta_hat = betanew,
    betas = betas,
    diffs = diffs
  )
  return(res)
}



if(FALSE){
  NUMSAMPLES <- 10000
  k <- 4 
  # sir <- my_sir_laplace(NUMSAMPLES, x, y0, a0, beta0_laplace = c(1,1), method_beta = 'optim', proposal = 'tstudent')
}
