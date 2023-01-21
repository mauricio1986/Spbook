######## SLM by MLE #####

# Create log-likelihood function for SLM ----
sml_ll <- function(theta, y, X, W, gradient = TRUE, hessian = TRUE){
  # Global
  K <- ncol(X)
  N <- nrow(X)
  
  # Extract parameters
  betas  <- theta[1:K]
  rho    <- theta[K + 1]
  sig.sq <- theta[K + 2]
  
  # Make residuals
  A   <- diag(N) -  rho * W
  Ay  <- A %*% y
  Xb  <- X %*% betas
  res <- Ay - Xb
  
  # Make log-likelihood
  detA <- det(A)
  ll   <- -0.5 * N * log(2 * pi * sig.sq) - 0.5 * crossprod(res) / sig.sq + log(detA)
  
  # Gradient
  if (gradient){
    C           <-  W %*% solve(A)
    grad.betas  <- (1 / sig.sq) * t(X) %*% res
    grad.rho    <- - sum(diag(C)) + (1 / sig.sq) * t(res) %*% W %*% y
    grad.sig.sq <- (1 / (2 * sig.sq ^2)) * (t(res) %*% res - N * sig.sq)
    attr(ll, 'gradient') <- c(grad.betas, grad.rho, grad.sig.sq)
  }
  # Hessian
  if (hessian){
    H    <- matrix(NA, nrow = (K + 2), ncol = (K + 2))
    h_bb <- - (1 / sig.sq) * t(X) %*% X
    h_bs <- - (1 / sig.sq ^ 2) * t(X) %*% res
    h_br <- - (1 / sig.sq) * t(X) %*% W %*% y 
    h_ss <- (N / (2 * sig.sq ^ 2)) - (1 / sig.sq ^ 3) * t(res) %*% res
    h_sr <-  - t(res) %*% W %*% y / sig.sq ^ 2
    h_rr <- - sum(diag(C %*% C)) - (1 / sig.sq) * (t(y) %*% t(W) %*% W %*% y)
    H[1:K, 1:K]     <- h_bb
    H[1:K, K + 1]   <- h_bs
    H[1:K, K + 2]   <- h_br
    H[K + 1, 1:K]   <- t(h_bs)
    H[K + 1, K + 1] <- h_ss
    H[K + 1, K + 2] <- h_sr
    H[K + 2, 1:K]   <- t(h_br)
    H[K + 2, K + 1] <- h_sr
    H[K + 2, K + 2] <- h_rr
    attr(ll, 'hessian') <- H
  }
  return(ll)
}

slm_ml <- function(formula, data, listw, 
                   gradient = TRUE, 
                   hessian  = TRUE, ...){
  require("maxLik")
  require("spdep")
  # Model Frame: This part is standard in R to obtain
  #              the variables using formula and data argument.
  callT <- match.call(expand.dots = TRUE)
  mf <- callT
  m  <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame()) # final model frame
  nframe     <- length(sys.calls())
  
  # Get variables and globals
  y  <- model.response(mf)        # Get dependent variable from mf
  X  <- model.matrix(formula, mf) # Get X from mf
  W  <- listw2mat(listw)          # listw to matrix
  K  <- ncol(X)
  
  # Starting values
  b_hat <- coef(lm(y ~ X - 1))
  start <- c(b_hat, 0, 1)
  names(start) <- c(colnames(X), "rho", "sig.sq")
    
  # Restricted optimization: A %*% theta + B >= 0: Constraint rho and sigma2
  sym          <- all(W == t(W))
  omega        <- eigen(W, only.values = TRUE, symmetric = sym)
  lambda_space <- if (is.complex(omega$values)) 1 / range(Re(omega$values)) else 1 / range(omega$values)
  
  A <- rbind(c(rep(0, K), 1, 0),
             c(rep(0, K), -1, 0), 
             c(rep(0, K), 0, 1))
  B <- c(-1L * (lambda_space[1] + sqrt(.Machine$double.eps)), 
                lambda_space[2] - sqrt(.Machine$double.eps), 
         -1L* sqrt(.Machine$double.eps))
  callT$constraints <- list(ineqA = A, ineqB = B)
  
  # Optimization default controls if not added by user
  if (is.null(callT$method))  callT$method  <- 'bfgs'
  if (is.null(callT$iterlim)) callT$iterlim <- 100000
  opt <- callT
  m <- match(c('method', 'print.level', 'iterlim',
               'tol', 'ftol', 'steptol', 'fixed', 'constraints', 
               'control', 'finalHessian', 'reltol', 'rho', 'outer.iterations', 'outer.eps'),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  opt$start     <- start
  opt[[1]]      <- as.name('maxLik')
  opt$logLik    <- as.name('sml_ll')
  opt$gradient  <- gradient
  opt$hessian   <- hessian
  opt[c('y', 'W', 'X')] <- list(as.name('y'), 
                                as.name('W'), 
                                as.name('X'))
  out <- eval(opt, sys.frame(which = nframe))
  return(out)
}

# Load data ----
data(oldcol, package="spdep")
listw <- spdep::nb2listw(COL.nb, style = "W")

# Use our function ----
test1 <- slm_ml(CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw)
summary(test1)

# Use lagsarlm from spatialreg
library("spatialreg")
sreg <- lagsarlm(CRIME ~ INC + HOVAL, data = COL.OLD, listw = listw)
summary(sreg)

