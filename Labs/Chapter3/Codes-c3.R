## ----sim_set_up1----------------------------------------------------------------------------------------------------------------------------------------
# Global parameters
library("spdep")
library("spatialreg")
set.seed(123)                                   # Set seed
S       <- 100                                  # Number of simulations
n       <- 225                                  # Spatial units
rho     <- 0.7                                  # True rho
w       <- cell2nb(sqrt(n), sqrt(n))            # Create artificial W matrix
iw      <- invIrM(w, rho)                       # Compute inverse of (I - rho*W)
rho_hat <- vector(mode = "numeric", length = S) # Vector to save results.


## ----loop_sim_1-----------------------------------------------------------------------------------------------------------------------------------------
# Loop for simulation
for (s in 1:S) {
  e <- rnorm(n, mean = 0 , sd = 1) # Create error term
  y <- iw %*% e                    # True DGP
  Wy <- lag.listw(nb2listw(w), y)  # Create spatial lag
  out <- lm(y ~ Wy)                # Estimate OLS
  rho_hat[s] <- coef(out)["Wy"]    # Save results
}


## ----sum-loop-1-----------------------------------------------------------------------------------------------------------------------------------------
# Summary of rho_hat
summary(rho_hat)


## ----ols-rho-sim-F, eval = FALSE------------------------------------------------------------------------------------------------------------------------
# # Plot density of estimated rho_hat.
# plot(density(rho_hat),
#      xlab = expression(hat(rho)),
#      main = "")
# abline(v = rho, col = "red")


## ----ols-rho-sim, echo = FALSE, message = FALSE, fig.align='center', out.width = '8cm', out.height = '8cm'----------------------------------------------
plot(density(rho_hat),
     xlab = expression(hat(rho)), 
     main = "",
     xlim = c(0.6, 1.3)
     )
abline(v = rho, col = "red")


## ----message = FALSE, warning=FALSE---------------------------------------------------------------------------------------------------------------------
# Load packages
library("spdep")
library("tmap")
library("spatialreg")
library("memisc")            # Package for tables
library("RColorBrewer") 
library("classInt")
source("getSummary.sarlm.R") # Function for spdep models


## ----load-columbus, warning = FALSE---------------------------------------------------------------------------------------------------------------------
# Load data
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], 
                           quiet = TRUE)
col.gal.nb <- poly2nb(columbus)


## ----spatial-crime-false, eval = FALSE------------------------------------------------------------------------------------------------------------------
# # Spatial distribution of crime
# tm_shape(columbus) +
#   tm_polygons("CRIME",
#               title = "Crime Status",
#               parlette = "RdYBu"
#               ) +
#   tm_layout(frame = FALSE)


## ----spatial-crime, echo = FALSE, message = FALSE, fig.align='center', out.width = '9cm', out.height = '9cm'--------------------------------------------
par(mar = c(0.1, 0.1, 0.1, 0.1))
tm_shape(columbus) +
  tm_polygons("CRIME",
              title = "Crime Status",
              parlette = "RdYBu"
              ) +
  tm_layout(frame = FALSE) 


## ----moran-oldcol---------------------------------------------------------------------------------------------------------------------------------------
# Moran's I test
set.seed(1234)
listw <- nb2listw(col.gal.nb, style = "W")
moran.mc(columbus$CRIME, listw = listw, 
           nsim = 99, alternative = 'greater')


## ----models-crime---------------------------------------------------------------------------------------------------------------------------------------
# Models
columbus$lag.INC   <- lag.listw(listw, 
                         columbus$INC)   # Create spatial lag of INC
columbus$lag.HOVAL <- lag.listw(listw, 
                         columbus$HOVAL) # Create spatial lag of HOVAL
ols <- lm(CRIME ~ INC + HOVAL, 
          data =  columbus)        
slx <- lm(CRIME ~ INC + HOVAL + lag.INC + lag.HOVAL, 
          data =  columbus)
slm <- lagsarlm(CRIME ~ INC + HOVAL, 
                data = columbus,
                listw, 
                method = "eigen")
sdm <- lagsarlm(CRIME ~ INC + HOVAL, 
                data = columbus,
                listw, 
                method = "eigen",
                type = "mixed")
sem <- errorsarlm(CRIME ~ INC + HOVAL, 
                data = columbus,
                listw,
                method = "eigen")
sac <- sacsarlm(CRIME ~ INC + HOVAL, 
                data = columbus,
                listw,
                method = "eigen")


## ----slx-spdep, eval = FALSE----------------------------------------------------------------------------------------------------------------------------
# slx2 <- lmSLX(CRIME ~ INC + HOVAL,
#               data = columbus,
#               listw)
# summary(slx2)


## ----echo = FALSE, results = 'asis', warning=FALSE------------------------------------------------------------------------------------------------------
table_1 <- mtable("OLS" = ols,
                  "SLX" = slx,
                  "SLM" = slm,
                  "SDM" = sdm,
                  "SEM" = sem, 
                  "SAC" = sac, 
       summary.stats = c("AIC", "N"),
       coef.style = "default")
table_1 <- relabel(table_1,
                   "(Intercept)" = "\\emph{Constant}",
                   "rho" = "$\\rho$",
                   "lambda" = "$\\lambda$",
                   "lag.INC" = "$W.INC$",
                   "lag.HOVAL" = "$W.HOVAL$") 
toLatex(table_1, compact = TRUE, useBooktabs =  TRUE)


## ----y-hat-pre------------------------------------------------------------------------------------------------------------------------------------------
# The predicted values
rho       <- slm$rho                                # Estimated rho from SLM model
beta_hat  <- coef(slm)[-1]                          # Estimated parameters
A         <- invIrW(listw, rho = rho)               # (I - rho*W)^{-1}
X         <- cbind(1, columbus$INC, columbus$HOVAL) # Matrix of observed variables
y_hat_pre <- A %*% crossprod(t(X), beta_hat)        # y hat


## ----y-hat-post-----------------------------------------------------------------------------------------------------------------------------------------
# The post-predicted values
col_new <- columbus # copy the data frame

# Change the income value
col_new[30, "INC"] <- 14.906

# The predicted values
X_d        <- cbind(1, col_new$INC, col_new$HOVAL)
y_hat_post <- A %*% crossprod(t(X_d), beta_hat)


## ----diff-predicts--------------------------------------------------------------------------------------------------------------------------------------
# The difference
delta_y         <- y_hat_post - y_hat_pre
col_new$delta_y <- delta_y

# Show the effects
summary(delta_y)
sum(delta_y)


## ----predicted-effect-evalF, eval = FALSE---------------------------------------------------------------------------------------------------------------
# tm_shape(col_new) +
#   tm_polygons(col = "delta_y",
#               title = "Impact",
#               breaks = c(-Inf,-0.05, Inf),
#               palette = c("red", "blue")
#               ) +
#   tm_layout(frame = FALSE)


## ----predicted-effect, echo = FALSE, message = FALSE, fig.align='center', out.width = '10cm', out.height = '10cm'---------------------------------------
# Breaks
tm_shape(col_new) +
  tm_polygons(col = "delta_y",
              title = "Impact",
              breaks = c(-Inf,-0.05, Inf),
              palette = c("red", "blue")
              ) +
  tm_layout(frame = FALSE) 


## ----predicted-effect2-evalF, eval = FALSE--------------------------------------------------------------------------------------------------------------
# # Plot the magnitude of the ME
# pal5    <- brewer.pal(6, "Spectral")
# cats5   <- classIntervals(col_new$delta_y, n = 5, style = "jenks")
# colors5 <- findColours(cats5, pal5)
# plot(col_new, col = colors5)
# legend("topleft", legend = round(cats5$brks, 2), fill = pal5, bty = "n")


## ----predicted-effect2, echo = FALSE, message = FALSE, fig.align='center', out.width = '10cm', out.height = '10cm'--------------------------------------
tm_shape(col_new) +
  tm_polygons(col = "delta_y",
              title = "Impacts",
              palette = "OrRd", 
              style = "quantile",
              ) +
  tm_layout(frame = FALSE) 


## ----using-impacts--------------------------------------------------------------------------------------------------------------------------------------
spatialreg:::impacts.Sarlm(slm, listw = listw)


## ----impacts-by-hand------------------------------------------------------------------------------------------------------------------------------------
## Construct S_r(W) = A(W)^-1 (I * beta_r + W * theta_r)
Ibeta <- diag(length(listw$neighbours)) *  coef(slm)["INC"] 
S <- A %*% Ibeta

ADI <- sum(diag(S)) / nrow(A)
ADI

n     <- length(listw$neighbours)
Total <- crossprod(rep(1, n), S) %*% rep(1, n) / n
Total

Indirect <- Total - ADI
Indirect


## ----impacts-with-se------------------------------------------------------------------------------------------------------------------------------------
# Compute standard errors of impacts
im_obj <- spatialreg:::impacts.Sarlm(slm, listw = listw, R = 200)
summary(im_obj, zstats = TRUE, short = TRUE)


## ----compute-me1----------------------------------------------------------------------------------------------------------------------------------------
# Impacts using traces. 
W <- as(nb2listw(col.gal.nb, style = "W"), "CsparseMatrix")
trMC <- trW(W, type = "MC")
im <- spatialreg:::impacts.Sarlm(slm, tr = trMC, R = 100)
summary(im, zstats =  TRUE, short = TRUE)


## ----compute-cumme1-------------------------------------------------------------------------------------------------------------------------------------
# Cumulative impacts
im2   <- spatialreg:::impacts.Sarlm(slm, tr = trMC, R = 100, Q = 5)
sums2 <- summary(im2, zstats = TRUE, reportQ = TRUE, short =  TRUE)
sums2


## ----log-like-func--------------------------------------------------------------------------------------------------------------------------------------
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
    grad.sig.sq <- (1 / (2 * sig.sq ^2 )) * (t(res) %*% res - N * sig.sq)
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
    H[1:K, K + 1]   <- h_br
    H[1:K, K + 2]   <- h_bs
    H[K + 1, 1:K]   <- t(h_br)
    H[K + 1, K + 1] <- h_rr
    H[K + 1, K + 2] <- h_sr
    H[K + 2, 1:K]   <- t(h_bs)
    H[K + 2, K + 1] <- h_sr
    H[K + 2, K + 2] <- h_ss
    attr(ll, 'hessian') <- H
  }
  return(ll)
}


## ----slm-ml-function, message=FALSE---------------------------------------------------------------------------------------------------------------------
library("maxLik")
slm.ml <- function(formula, data, W, 
                   gradient = TRUE, 
                   hessian  = TRUE, ...){
  require("maxLik")
  # Model Frame: This part is standard in R to obtain
  #              the variables using formula and data argument.
  callT    <- match.call(expand.dots = TRUE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf       <- eval(mf, parent.frame()) # final model frame
  nframe   <- length(sys.calls())
  
  # Get variables and globals
  y  <- model.response(mf)        # Get dependent variable from mf
  X  <- model.matrix(formula, mf) # Get X from mf
  K  <- ncol(X)
  
  # Starting values
  ols.init    <- lm(y ~ X - 1)
  b.init      <- coef(ols.init)
  sigma2.init <- sum(residuals(ols.init)^2) / ols.init$df.residual
  rho.init    <- cor(W %*% y, y)
  start       <- c(b.init, rho.init, sigma2.init)
  names(start) <- c(colnames(X), "rho", "sig.sq")
  
  # Optimization default controls if not added by user
  if (is.null(callT$method))  callT$method  <- 'bfgs'
  if (is.null(callT$iterlim)) callT$iterlim <- 100000
    
  # Restricted optimization if BFGS: A %*% theta + B >= 0: Constraint rho and sigma2
  if (callT$method == "bfgs"){
    sym          <- all(W == t(W))
    omega        <- eigen(W, only.values = TRUE, symmetric = sym)
    lambda_space <- if (is.complex(omega$values)) 1 / range(Re(omega$values)) else 1 / range(omega$values)
    A <- rbind(c(rep(0, K), 1, 0),
               c(rep(0, K), -1, 0), 
               c(rep(0, K), 0, 1))
    B <- c(-1L * (lambda_space[1] + sqrt(.Machine$double.eps)), 
                  lambda_space[2] - sqrt(.Machine$double.eps), 
           -1L * sqrt(.Machine$double.eps))
   callT$constraints <- list(ineqA = A, ineqB = B)
  }

  # Optimization
  opt <- callT
  m <- match(c('method', 'print.level', 'iterlim',
               'tol', 'ftol', 'steptol', 'fixed', 'constraints', 
               'control', 'finalHessian', 'reltol', 'rho', 
               'outer.iterations', 'outer.eps'),
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


## ----generate-DGP-ml------------------------------------------------------------------------------------------------------------------------------------
# Generate DGP
set.seed(1)
n      <- 529
rho    <- 0.6
W.nb2  <- cell2nb(sqrt(n), sqrt(n))
W      <- nb2mat(W.nb2)

# Exogenous variables
x1     <- rnorm(n)
x2     <- rnorm(n)
x3     <- rnorm(n)

# DGP parameters
b0 <- 0 ; b1 <- -1; b2 <- 0; b3 <- 1
sigma2 <- 2
epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma2))

# Simulate the dependent variable
y <- solve(diag(n) -  rho * W) %*% (b0 + b1*x1 + b2*x2 + b3*x3 + epsilon)

# Data as data.frame
data <- as.data.frame(cbind(y, x1, x2, x3))
names(data) <- c("y", "x1", "x2", "x3")


## ----check-func1, message = FALSE, cache = TRUE---------------------------------------------------------------------------------------------------------
# Use our function 
start <- Sys.time()
sml.mle <- slm.ml(y ~ x1 + x2 +  x3, data = data, W = W)
summary(sml.mle)
print(Sys.time()- start)

start <- Sys.time()
sml.mle.nr <- slm.ml(y ~ x1 + x2 +  x3, data = data, W = W, method = "nr")
summary(sml.mle.nr)
print(Sys.time()- start)


## ----conc-ml--------------------------------------------------------------------------------------------------------------------------------------------
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


## ----secondfunc-----------------------------------------------------------------------------------------------------------------------------------------
sar.mle.con <- function(formula, data, W)
{
  # Model Frame: This part is standard in R to obtain
  # the variables using formula and data argument.
  callT <- match.call(expand.dots = TRUE)
  mf <- callT
  m  <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame()) # final model frame
  
  # Get variables and globals
  y  <- model.response(mf)        # Get dependent variable from mf
  X  <- model.matrix(formula, mf) # Get X from mf
  n  <- nrow(X)                   # Number of spatial units
  k  <- ncol(X)                   # Number of regressors
  Wy <- W %*% y                   # Spatial lag
  
  # Generate auxiliary regressions 
  # See Algorithm 3.1
  ols_0 <- lm(y ~ X - 1)
  ols_L <- lm(Wy ~ X - 1)
  e_0   <- residuals(ols_0)
  e_L   <- residuals(ols_L)
  
  # Get eigenvalues to constraint the optimization
  sym          <- all(W == t(W))
  omega        <- eigen(W, only.values = TRUE, symmetric = sym)
  
  # Maximize concentrated log-likelihood
  rho_space <- if (is.complex(omega$values)) 1 / range(Re(omega$values)) else 1 / range(omega$values)
  opt_lc <- optimize(f = logLik_sar,   # This function is below
                     lower = rho_space[1] + .Machine$double.eps,
                     upper = rho_space[2] - .Machine$double.eps,
                     maximum = TRUE, 
                     e_0 = e_0, e_L = e_L, omega = omega$values, n = n, 
                     tol = .Machine$double.eps)
  # Obtain rho_hat from concentrated log-likelihood
  rho_hat <- opt_lc$maximum
  
  # Generate estimates
  A          <- (diag(n) - rho_hat * W)
  Ay         <- crossprod(t(A), y)
  beta_hat   <- solve(crossprod(X)) %*% crossprod(X, Ay) 
  error      <- Ay - crossprod(t(X), beta_hat)
  sigma2_hat <- crossprod(error) / n                    
  
  # Save results
  out <- structure(
    list(
      callT = callT,
      rho_hat = rho_hat, 
      beta_hat = beta_hat,
      sigma2_hat = sigma2_hat, 
      A = A, 
      W = W,
      X = X, 
      omega = omega 
    ),
    class = "slmc.mle"
  )
  
return(out)
}


## ----vcov.slmc.mle--------------------------------------------------------------------------------------------------------------------------------------
vcov.slmc.mle <- function(object, ...){
  rho_hat    <- object$rho_hat
  beta_hat   <- object$beta_hat
  sigma2_hat <- object$sigma2_hat
  A          <- object$A
  W          <- object$W
  X          <- object$X
  omega      <- object$omega$values
  k          <- ncol(X)
  
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

  return(solve(Hessian))
}


## ----summary-smlc.mle-----------------------------------------------------------------------------------------------------------------------------------
# S3 methods for summary
summary.slmc.mle <- function(object, 
                              table = TRUE, 
                              digits = max(3, .Options$digits - 3), 
                              ...){
    X       <- object$X
    n       <- nrow(X)
    k       <- ncol(X)
    df      <- n - (k + 1)
    b       <- c(object$beta_hat, object$sigma2_hat, object$rho_hat)
    names(b) <- c(colnames(X), "sigma2", "Wy")
    std.err <- sqrt(diag(vcov(object)))
    z       <- b / std.err
    p       <- 2 * pt(-abs(z), df = df)
    CoefTable <- cbind(b, std.err, z, p)
    colnames(CoefTable) <- c("Estimate", "Std.Error", "t-value", "Pr(>|t|)")
    result <- structure(
      list(
        CoefTable = CoefTable, 
        digits    = digits, 
        call      = object$call),
      class = 'summary.slmc.mle'
    )
    return(result)
}

print.summary.slmc.mle <- function(x, 
                                   digits = x$digits, 
                                   na.print = "", 
                                   symbolic.cor = p > 4, 
                                   signif.stars = getOption("show.signif.stars"), 
                                   ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
      
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digit = digits, P.value = TRUE, has.Pvalue = TRUE)
  invisible(NULL)
}


## -------------------------------------------------------------------------------------------------------------------------------------------------------
slm2 <- sar.mle.con(y ~ x1 + x2 + x3, data = data, W = W)
summary(slm2)


## ----echo = FALSE, results = 'asis', warning=FALSE------------------------------------------------------------------------------------------------------
sarlm <- lagsarlm(y~ x1 + x2 + x3, 
                   data = data, 
                   listw = mat2listw(W, style = "W"), 
                   method = "eigen")
coefs <- cbind(sml.mle$estimate, 
            sml.mle.nr$estimate, 
           c(slm2$beta_hat, slm2$rho_hat, slm2$sigma2_hat),
           c(sarlm$coefficients, sarlm$rho, sarlm$s2)
           )
ses <- cbind(sqrt(diag(vcov(sml.mle))), 
             sqrt(diag(vcov(sml.mle.nr))),
             sqrt(diag(vcov(slm2)))[c(1, 2, 3, 4, 6, 5)], 
             c(sqrt(diag(vcov(sarlm)))[c(2, 3, 4, 5, 1)], 0)
             )
Tab <- cbind(coefs, ses)
rownames(Tab) <- c("b0", "b1", "b2", "b3", "rho", "sigma2")
colnames(Tab) <- c("F-BFGS", "F-NR", "CMLE", "R", 
                   "F-BFGS", "F-NR", "CMLE", "R")
toLatex(Tab, compact = TRUE, useBooktabs =  TRUE, digits = 5)

