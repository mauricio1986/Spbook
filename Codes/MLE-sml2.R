sarML <- function(formula, data, listw)
{
  require("spdep")
  # Model Frame: This part is standard in R to obtain
  # the variables using formula and data argument.
  callT <- match.call(expand.dots = TRUE)
  mf <- callT
  m  <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame()) # final model frame
  
  # Get variables and Globals
  y  <- model.response(mf)        # Get dependent variable from mf
  X  <- model.matrix(formula, mf) # Get X from mf
  n  <- nrow(X)                   # Number of spatial units
  k  <- ncol(X)                   # Number of regressors
  Wy <- lag.listw(listw, y)       # Spatial lag
  W  <- listw2mat(listw)          # listw to matrix
  
  # Generate auxiliary regressions 
  # See Algorithm 3.1
  ols_0 <- lm(y ~ X - 1)
  ols_L <- lm(Wy ~ X - 1)
  e_0   <- residuals(ols_0)
  e_L   <- residuals(ols_L)
  
  # Get eigenvalues to constraint the optimization
  omega <- eigenw(listw)
  
  # Maximize concentrated log-likelihood
  rho_space <- if (is.complex(omega)) 1 / range(Re(eig)) else 1 / range(omega)
  opt_lc <- optimize(f = logLik_sar,   # This function is below
                     lower = rho_space[1] + .Machine$double.eps,
                     upper = rho_space[2] - .Machine$double.eps,
                     maximum = TRUE, 
                     e_0 = e_0, e_L = e_L, omega = omega, n = n)
  # Obtain rho_hat from concentrated log-likelihood
  rho_hat <- opt_lc$maximum
  
  # Generate estimates
  A          <- (diag(n) - rho_hat * W)
  Ay         <- crossprod(t(A), y)
  beta_hat   <- solve(crossprod(X)) %*% crossprod(X, Ay) # See Equation (3.25)
  error      <- Ay - crossprod(t(X), beta_hat)
  sigma2_hat <- crossprod(error) / n                     # See Equation (3.26)
  
  # Hessian
  C       <- crossprod(t(W), solve(A)) # C = WA^{-1}
  alpha   <-  sum(omega ^ 2 / ((1 - rho_hat * omega) ^ 2))
  if (is.complex(alpha)) alpha <- Re(alpha)
  b_b     <- drop(1 / sigma2_hat) * crossprod(X) # k * k
  b_rho   <- drop(1 / sigma2_hat) * (t(X) %*% C %*% X %*% beta_hat) # k * 1
  sig_sig <- n / (2 * sigma2_hat ^ 2) # 1 * 1
  sig_rho <- drop(1 / sigma2_hat) * sum(diag(C)) # 1 * 1
  rho_rho <- sum(diag(crossprod(C))) +  alpha +
    drop(1 / sigma2_hat) * crossprod(C %*% X %*% beta_hat) # 1*1
  row_1   <- cbind(b_b, rep(0, k), b_rho)
  row_2   <- cbind(t(rep(0, k)), sig_sig, sig_rho)
  row_3   <- cbind(t(b_rho), sig_rho, rho_rho)
  Hessian <- rbind(row_1, row_2, row_3)
  std.err <- sqrt(diag(solve(Hessian)))
  
  # Table of coefficients
  all_names          <- c(colnames(X), "sigma2", "rho")
  all_coef           <- c(beta_hat, sigma2_hat, rho_hat)
  z                  <- all_coef / std.err
  p                  <- pnorm(abs(z), lower.tail = FALSE) * 2 
  sar_table          <- cbind(all_coef, std.err, z, p)
  cat(paste("\nEstimates from SAR Model \n\n"))
  colnames(sar_table) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  rownames(sar_table) <- all_names
  printCoefmat(sar_table)
}

logLik_sar <- function(rho, e_0, e_L, omega, n)
{
  # This function returns the concentrated log L for maximization
  
  #Generate determinant using Ord's approximation
  det    <- if (is.complex(omega)) Re(prod(1 - rho * omega)) else prod(1 - rho * omega)
  e_diff <- e_0 - rho * e_L
  sigma2 <- crossprod(e_diff) / n
  
  #Log-Likelihood function
  l_c    <- - (n / 2) - (n / 2) * log(2 * pi) - (n / 2) * log(sigma2) + log(det)
  return(l_c)
}