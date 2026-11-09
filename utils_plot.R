library(ggplot2)
library(dplyr)
library(pracma)
library(tidyr)
source('utils_beta_sampling.R')



plot_1d_beta_samples_from_object <- function(sampler_output, x, y0, a0, tau, family,
                                             method_name, D_fun = D,
                                             logscale = FALSE, grid_points = 1000) {
  library(ggplot2)
  
  # Extract samples
  beta_samples <- sampler_output$beta_samples
  if (!is.numeric(beta_samples)) {
    stop("Expected numeric vector for 1D beta_samples.")
  }
  
  # Create grid
  beta_grid <- seq(min(beta_samples), max(beta_samples), length.out = grid_points)
  
  # Evaluate D over grid
  Dvalues <- sapply(beta_grid, function(beta) {
    D_fun(y0 = y0, a0 = a0, tau = tau, x = x, beta = beta, family = family, logscale = logscale)
  })
  
  # Normalize
  dx <- diff(beta_grid)[1]
  Z <- sum(Dvalues) * dx
  Dvalues <- Dvalues / Z
  
  # Create data frames
  d_kernel <- data.frame(beta_grid = beta_grid, Dvalues = Dvalues)
  dfsir <- as.data.frame(sampler_output)
  
  # Plot
  p <- ggplot(dfsir, aes(x = beta_samples)) + 
    ggtitle(paste("Samples from", method_name)) +
    geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "white", bins = 30) +
    geom_density(alpha = .2, fill = "#FF6666") +
    geom_line(data = d_kernel, aes(x = beta_grid, y = Dvalues), colour = "blue") +
    theme_minimal()
  
  return(p)
}




plot_2d_beta_samples_from_object <- function(sampler_output, x, y0, a0, tau, family,
                                             method_name = "SIR", D_fun = D,
                                             logscale = FALSE, bins = 10,
                                             show_density = 0, grid_points = 100) {
  # Extract beta_samples and convert to data.frame
  beta_samples <- as.data.frame(sampler_output$beta_samples)
  if (ncol(beta_samples) != 2) {
    stop("This function expects 2D beta samples (i.e., two columns).")
  }
  colnames(beta_samples) <- c("beta1", "beta2")
  
  # Grid setup
  b1_seq <- seq(min(sampler_output$beta_samples[,1]), max(sampler_output$beta_samples[,1]), length.out = grid_points)
  b2_seq <- seq(min(sampler_output$beta_samples[,2]), max(sampler_output$beta_samples[,2]), length.out = grid_points)
  grid <- expand.grid(b1 = b1_seq, b2 = b2_seq)
  
  # Compute D on grid
  grid$Dval <- apply(grid, 1, function(b) {
    D_fun(y0 = y0, a0 = a0, tau = tau, x = x, beta = b, family = family, logscale = logscale)
  })
  
  # Normalize Dval
  dx <- diff(b1_seq)[1]
  dy <- diff(b2_seq)[1]
  area_element <- dx * dy
  Z <- sum(grid$Dval) * area_element
  grid$Dval <- grid$Dval / Z
  
  # Prepare sample data
  df <- as.data.frame(sampler_output$beta_samples)
  colnames(df) <- c("beta1", "beta2")
  
  # Plot
  p <- ggplot() +
    geom_point(data = df, aes(x = beta1, y = beta2), alpha = 0.3, size = 0.7, color = "black") +
    geom_contour(data = grid, aes(x = b1, y = b2, z = Dval), colour = "blue", bins = bins) +
    labs(title = paste("Samples from", method_name),
         x = expression(beta[1]), y = expression(beta[2])) +
    theme_minimal()
  
  if (show_density) {
    p <- p + stat_density_2d(data = beta_samples,
                             aes(x = beta1, y = beta2, fill = after_stat(level)),
                             geom = "polygon", alpha = 0.3)
  }
  
  return(p)
}

plot_posterior_betas <- function(results, y0, a0, tau, family, x, 
                                 true_values=NULL, figure_folder) {
  
  p <- ncol(results$beta)
  n <- nrow(x)
  
  beta_fixed <- apply(results$beta, 2, median)
  
  if (!is.null(figure_folder)) {
    pdf(file = file.path(figure_folder, paste0("posterior_betas_p", p, "_n", n, ".pdf")), width = 8, height = 6)
    par(mfrow = c(1,1))
  } else {
    par(mfrow = c(ceiling(p / 2), 2))
  }
  
  for (j in 1:p) {

    hist(results$beta[,j], main = ' ', breaks = 40, prob = TRUE, xlab = bquote(beta[.(j)]))
    lines(density(results$beta[,j]), col = "blue", lwd = 2)

    # marginal density:
    xgrid <- seq(min(results$beta[,j]), max(results$beta[,j]), length = 40)
    
    D_slice <- function(bj) {
      beta_tmp <- beta_fixed
      beta_tmp[j] <- bj
      as.numeric(D(y0 = y0, a0 = a0, tau = tau, x = x, beta = beta_tmp, family = family, logscale = FALSE))
    }
    
    mean_post <- mean(results$beta[, j])
    sd_post <- sd(results$beta[, j])
    xgrid_norm <- seq(mean_post - 5 * sd_post, mean_post + 5 * sd_post, length.out = 500)
    #xgrid_norm <- seq(-10, 10, length.out = 500)
    
    prior_vals_raw <- vapply(xgrid, D_slice, numeric(1)) # for plot 
    prior_vals_normgrid <- vapply(xgrid_norm, D_slice, numeric(1))  # for normalization
    norm_const <- trapz(xgrid_norm, prior_vals_normgrid)
    
    if (!is.na(norm_const) && norm_const > 0) {
      prior_vals <- prior_vals_raw / norm_const
      lines(xgrid, prior_vals, col = "red", lwd = 2)
    } else {
      print(paste('Not showing prior because normalization constant = ', norm_const))
    }

    if (!is.null(true_values)) {
      abline(v = true_values[j], col = "green", lwd = 2)
      legend('topright', legend = c("posterior", "prior", paste("true = ", round(true_values[j],2))),
           col = c("blue", "red", "green"), lty = 1, cex = 0.8)
    } else{
      legend('topright', legend = c("posterior", "prior"),
             col = c("blue", "red"), lty = 1, cex = 0.8)
    }
  }
  
  if (is.null(figure_folder)) {
    mtext(expression(paste('Posterior distribution of ', beta)),  side = 3, line=-2, outer = TRUE)
    par(mfrow = c(1, 1))
  } else {
    dev.off()
  }

}



plot_posterior_z <- function(results, true_values=NULL, figure_folder, n=NULL) {
  
  p <- ncol(results$beta)
  
  if (!is.null(figure_folder)) {
    pdf(file = file.path(figure_folder, paste0("posterior_z_p", p, "_n", n, ".pdf")), width = 8, height = 6)
    par(mfrow = c(1,1))
  } else {
    par(mfrow = c(ceiling(p / 2), 2))
  }
  
  for (j in 1:p) {
    
    barplot(prop.table(table(factor(results$z[,j], levels=c(0,1)))),  xlab = bquote(z[.(j)]))
    if (!is.null(true_values)) {
      abline(v = true_values[j] + 1, col = "green", lwd = 2)
      legend("topright", legend = "True z", col = "green", lty = 1, lwd = 2, cex = 0.9)
    }
  }
  
  if (is.null(figure_folder)) {
    mtext(expression(paste('Posterior distribution of ', z)),  side = 3, line=-2, outer = TRUE)
    par(mfrow = c(1, 1))
    z_mean <- apply(results$z, 2, mean)
    z_median <- apply(results$z, 2, median)
    
    if (!is.null(true_values)) {
      z_table <- data.frame(
        True  = true_values,
        Post.Mean = z_mean, 
        Post.Median = z_median
      )
    } else {
      z_table <- data.frame(
        Post.Mean = z_mean, 
        Post.Median = z_median
      )
    }
    print(z_table)
  } else {
    dev.off()
  }
  
}


plot_posterior_taus <- function(results, shape, rate, true_value=NULL, figure_folder, n=NULL) {
  
  p <- ncol(results$beta)
  
  if (!is.null(figure_folder)) {
    pdf(file = file.path(figure_folder, paste0("posterior_tau_p", p, "_n", n, ".pdf")), width = 8, height = 6)
  }
  
  hist(results$tau, main = 'posterior tau', prob = TRUE, xlab = expression(paste(tau)))
  lines(density(results$tau), col = "blue", lwd = 2)
  
  xgrid <- seq(0, max(results$tau), length = 40)
  fun <- invgamma::dinvgamma(xgrid, shape=shape, rate=rate)
  lines(xgrid, fun, col = "red", lwd = 2)
  
  if (!is.null(true_value)) {
    abline(v=true_value, col="green", lwd = 2)
    legend('topright', legend=c("posterior", "prior", paste("true = ", round(true_value,2)) ), col=c("blue","red", "green"), lty=1, cex=0.9)
  } else {
    legend('topright', legend=c("posterior", "prior" ), col=c("blue","red"), lty=1, cex=0.9)
  }
  
  if (!is.null(figure_folder)) {
    dev.off()
  }
  
}
  
  

