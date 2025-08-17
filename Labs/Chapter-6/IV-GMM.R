## ----make.H------------------------------------------------------------------------------------------------------------------------------
# Function that creates WXs
make.H <- function(W, X, l = 3){
  # This function creates the instruments (WX, ...,W^lX)
  # Drop constant (if any)
  names.x <- colnames(X)
  if (names.x[1] == "(Intercept)") X <- matrix(X[ ,-1], 
                                               dim(X)[1], 
                                               dim(X)[2] - 1) #Drop first column
  names.x <- names.x[which(names.x != "(Intercept)")]
  # Create lagged X variables
  sq1 <- seq(1, ncol(X) * l, ncol(X))
  sq2 <- seq(ncol(X), ncol(X) * l, ncol(X))
  Hmat <- matrix(NA, nrow = nrow(X), ncol = ncol(X) * l)
  names.ins <- c()
  for (i in 1:l) {
    Hmat[, sq1[i]:sq2[i]] <- as.matrix(W %*% X)
    X <- Hmat[, sq1[i]:sq2[i]]
    names.ins <- c(names.ins, 
                   paste(paste(replicate(i, "W"), collapse = ""), 
                         names.x, sep = "*"))
  }
  colnames(Hmat) <- names.ins
  return(Hmat)
}


## ----generate-DGP-2sls-------------------------------------------------------------------------------------------------------------------
# Generate DGP
library("spatialreg")
library("spdep")
set.seed(1986)
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


## ----test.make.H-------------------------------------------------------------------------------------------------------------------------
X <- cbind(1, x1, x2, x3)
colnames(X) <- c("(Intercept)", "x1", "x2", "x3")
H <- make.H(W = W, X = X, l = 3) 
head(H)


## ----slm.2sls----------------------------------------------------------------------------------------------------------------------------
# Main function to estimate S2SLS estimator
slm.2sls <- function(formula, data, W, instruments = 2){
  # Model Frame
  callT    <- match.call(expand.dots = TRUE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf       <- eval(mf, parent.frame())
  
  # Get variables and globals
  y  <- model.response(mf)
  X  <- model.matrix(formula, mf)
  n  <- nrow(X)
  Wy <- W %*% y
  sn <- nrow(W)
  if (n != sn) stop("number of spatial units in W is different to the number of data")
  
  # Generate matrix of instruments H = [X, WX, ... ]
  # and select LI vars
  H <- cbind(X, make.H(W = W, X = X, l = instruments))
  H <- H[, qr(H)$pivot[seq_len(qr(H)$rank)]]
  
  # Get S2SLS estimates
  Z           <- cbind(X, Wy)
  colnames(Z) <- c(colnames(X), "Wy")

  # Create projection matrix
  HH    <- crossprod(H)
  PH    <- H %*% solve(HH) %*% t(H)
  
  # Compute S2SLS coefficients
  Z_hat  <- PH %*% Z
  b_2sls <- solve(crossprod(Z_hat)) %*% crossprod(Z_hat, y)
  y_hat  <- Z %*% b_2sls
  e_hat  <- y - y_hat
  
  # Save results
  results <- structure(
    list(
      coefficients = b_2sls, 
      call         = callT, 
      X            = X, 
      H            = H, 
      Z            = Z, 
      y            = y, 
      PH           = PH, 
      e_hat        = e_hat
    ), 
    class = 'mys2sls'
  )
}


## ----vcov.mys2sls------------------------------------------------------------------------------------------------------------------------
# S3 Method vcov
vcov.mys2sls <- function(object, tse = c("homo", "rob"), ...){
    tse    <- match.arg(tse)
    n      <- nrow(object$Z)
    df     <- n - ncol(object$Z)
    Q.HZ   <- (t(object$H) %*% object$Z) / n
    Q.HH.i <- solve(crossprod(object$H) / n)
    if (tse == "homo"){
      s2  <- crossprod(object$e_hat) / df
      var <- drop(s2) * solve(t(Q.HZ) %*% Q.HH.i %*% Q.HZ) / n
    } else {
      Delta.hat <- 0
      for (i in 1:nrow(object$Z)){
        Delta.hat <- Delta.hat + drop(object$e_hat[i] ^2) * tcrossprod(object$H[i, ])
      }
      bread  <- solve(t(Q.HZ) %*% Q.HH.i %*% Q.HZ)
      cheese <- t(Q.HZ) %*% Q.HH.i %*% (Delta.hat / n)  %*% Q.HH.i %*% Q.HZ
      var    <- (bread %*% cheese %*% bread) / n
    }
    return(var)
}


## ----summary.mys2sls---------------------------------------------------------------------------------------------------------------------
# S3 method summary
summary.mys2sls <- function(object, 
                              tse = c("homo", "rob"),
                              table = TRUE, 
                              digits = max(3, .Options$digits - 3), 
                              ...){
    tse       <- match.arg(tse)
    n         <- nrow(object$Z)
    df        <- n - ncol(object$Z)
    b         <- object$coefficients
    std.err   <- sqrt(diag(vcov(object, tse = tse)))
    z         <- b / std.err
    p         <- 2 * pt(-abs(z), df = df)
    CoefTable <- cbind(b, std.err, z, p)
    colnames(CoefTable) <- c("Estimate", "Std.Error", "t-value", "Pr(>|t|)")
    result <- structure(
      list(
        CoefTable = CoefTable, 
        digits    = digits, 
        call      = object$call),
      class = 'summary.mys2sls'
    )
    return(result)
}

# S3 method print.summary 
print.summary.mys2sls <- function(x, 
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


## ----estimate-2sls-code------------------------------------------------------------------------------------------------------------------
# Estimate S2SLS model
s2sls.e <- slm.2sls(y ~ x1 + x2 + x3, data = data, W = W, instruments = 2)
summary(s2sls.e)
summary(s2sls.e, tse = "rob")


## ----check-2sls-spatialreg---------------------------------------------------------------------------------------------------------------
spreg1 <- stsls(y ~ x1 + x2 + x3, data = data, listw = mat2listw(W, style = "W"))
spreg2 <- stsls(y ~ x1 + x2 + x3, data = data, listw = mat2listw(W, style = "W"),
                robust = TRUE, HC = "HC0")
cbind(as.numeric(coef(spreg1)[c(2, 3, 4, 5, 1)]), as.numeric(s2sls.e$coefficients))
cbind(sqrt(diag(spreg1$var)[c(2, 3, 4, 5, 1)]),
      sqrt(diag(vcov(s2sls.e))),
      sqrt(diag(spreg2$var)[c(2, 3, 4, 5, 1)]),
      sqrt(diag(vcov(s2sls.e, tse = "rob"))))


## ----slm.b2sls---------------------------------------------------------------------------------------------------------------------------
# Main function to estimate S2SLS estimator
slm.b2sls <- function(formula, data, W, instruments = 2){
  # Model frame setup
  callT    <- match.call(expand.dots = TRUE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf       <- eval(mf, parent.frame())
  
  # Get variables and check dimensions
  y  <- model.response(mf)
  X  <- model.matrix(formula, mf)
  n  <- nrow(X)
  Wy <- W %*% y
  sn <- nrow(W)
  if (n != sn) stop("number of spatial units in W is different to the number of data")
  
  # Obtain first (and consistent)  estimates
  H <- cbind(X, make.H(W = W, X = X, l = instruments))
  
  # Compute S2SLSE
  Z           <- cbind(X, Wy)
  colnames(Z) <- c(colnames(X), "Wy")
  HH          <- crossprod(H, H)
  PH          <- H %*% solve(HH) %*% t(H)
  Z_hat       <- PH %*% Z
  b_2sls      <- solve(crossprod(Z_hat)) %*% crossprod(Z_hat, y)
  
  # Second step: BS2SLS estimation
  beta.hat <- b_2sls[1:ncol(X)]
  rho.hat  <- drop(tail(b_2sls, n = 1L))
  H.lee    <- W %*% solve(diag(n) - rho.hat * W) %*% X %*% beta.hat
  H.star   <- cbind(X, H.lee)
  b_2sls   <- solve(crossprod(H.star, Z)) %*% crossprod(H.star, y)
  
  # Compute residuals
  y_hat  <- Z %*% b_2sls
  e_hat  <- y - y_hat
  
  # Save results
  results <- structure(
    list(
      coefficients = b_2sls, 
      call         = callT, 
      X            = X, 
      H            = H.star, 
      Z            = Z, 
      y            = y, 
      PH           = PH, 
      e_hat        = e_hat,
      W            = W
    ), 
    class = 'mybs2sls'
  )
}


## ----vcov.mybs2sls-----------------------------------------------------------------------------------------------------------------------
# S3 Method vcov
vcov.mybs2sls <- function(object, estimate = c("initial", "final"), ...){
  estimate <- match.arg(estimate)
  X        <- object$X
  n        <- nrow(X)
  k        <- ncol(X)
  df       <- n - (k + 1)
  s2       <- crossprod(object$e_hat) / df
  if (estimate == "final"){
    b       <- object$coefficients
    b.hat   <- b[1:k]
    rho.hat <- drop(tail(b, n = 1L))
    W       <- object$W
    H.lee   <- W %*% solve(diag(n) - rho.hat * W) %*% X %*% b.hat
    H.star  <- cbind(X, H.lee)
  } else {
    H.star  <- object$H
  }
  var    <- drop(s2) * solve(crossprod(H.star) / n) / n
  return(var)
}


## ----summary.mybs2sls--------------------------------------------------------------------------------------------------------------------
# S3 methods for summary
summary.mybs2sls <- function(object,
                             estimate = c("initial", "final"), 
                             table = TRUE, 
                             digits = max(3, .Options$digits - 3), 
                             ...){
  estimate  <- match.arg(estimate)
  n         <- nrow(object$Z)
  df        <- n - ncol(object$Z)
  b         <- object$coefficients
  std.err   <- sqrt(diag(vcov(object, estimate = estimate)))
  z         <- b / std.err
  p         <- 2 * pt(-abs(z), df = df)
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std.Error", "t-value", "Pr(>|t|)")
  result <- structure(
    list(
      CoefTable = CoefTable, 
      digits    = digits, 
      call      = object$call),
    class = 'summary.mybs2sls'
  )
  return(result)
}

print.summary.mybs2sls <- function(x, 
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


## ----b2sls-test1-------------------------------------------------------------------------------------------------------------------------
# Using the BS2SLS estimator
b2sls <- slm.b2sls(y ~ x1 + x2 + x3, data = data, W = W, instruments = 2)
summary(b2sls, estimate = "initial")
summary(b2sls, estimate = "final")


## ----moments.lee2007---------------------------------------------------------------------------------------------------------------------
# Trace function
tr <- function(A) return(sum(diag(A)))

# Create moments as in Lee 2007
moments.lee2007 <- function(theta, y, X, H, W){
  # y: vector of dependent variables
  # X: matrix of exogenous variables
  # H: matrix of linear instrument
  # W: Spatial weight matrix
  k        <- ncol(X)
  n        <- nrow(X)
  rho      <- theta[1L]
  beta     <- theta[(2L:(k + 1))]
  I        <- diag(n)
  S        <- I -  rho * W
  epsi     <- crossprod(t(S), y) - crossprod(t(X), beta)
  P1       <- W
  #W2       <- crossprod(t(W), W)
  W2       <- tcrossprod(W, W)
  P2       <- W2 - (tr(W2) / n) * I
  g.lin    <- crossprod(H, epsi)            # k_x * 1
  g.q1     <- crossprod(epsi, P1) %*% epsi  # 1*1
  g.q2     <- crossprod(epsi, P2) %*% epsi  # 1*1
  g        <- rbind(g.q1, g.q2, g.lin)
  
  # Gradient 
  P1s <- P1 + t(P1)
  P2s <- P2 + t(P2)
  # (2 * k_k)* (k + 1)
  D   <-  -1 * rbind(crossprod(epsi, P1s), 
                     crossprod(epsi, P2s), 
                     t(H)) %*% cbind(W %*% y, X)
  
  # Return results (note that they are divided by n)
  out <- list(g = g /n , D = D / n)
  return(out)
}


## ----Jmin--------------------------------------------------------------------------------------------------------------------------------
# Objective function to minimize
Qmin <- function(start, y, X, H, W, Psi, gradient){
  # Thus function returns the negative of:
  #  - objective function g'Psi g
  #  - gradient of g'Upsi g
  g.hat <- moments.lee2007(theta = start, y = y, X = X, H = H, W = W)
  Q     <- -1 * crossprod(g.hat$g, Psi) %*% g.hat$g
  if (gradient){
    D  <- g.hat$D # D.hat is (2 * k_k)* (k + 1)
    Gr <- -2 * crossprod(D, Psi) %*% g.hat$g
    attr(Q, 'gradient') <- as.vector(Gr)
  }
  return(Q)
}


## ----make.vmom---------------------------------------------------------------------------------------------------------------------------
# Create var-cov of moments: Omega
make.vmom <- function(b.hat, y, X, H, W){
  k        <- ncol(X)
  n        <- nrow(X)
  k_x      <- ncol(H)
  rho      <- b.hat[1L]
  beta     <- b.hat[(2L:(k + 1))]
  I        <- diag(n)
  S        <- I -  rho * W
  epsi     <- crossprod(t(S), y) - crossprod(t(X), beta)
  sigma2   <- as.numeric(crossprod(epsi) / n)
  P1       <- W
  #W2       <- crossprod(t(W), W)
  W2       <- tcrossprod(t(W), W)
  P2       <- W2 - (tr(W2) / n) * I
  P1s      <- P1 + t(P1)
  P2s      <- P2 + t(P2)
  
  # Construct V: (2 + k.x) * (2 + k.x)
  V22 <- (1 / sigma2) * crossprod(H)         # k_x + k_x
  Delta <- matrix(0, nrow = 2, ncol = 2)
  Delta[1, 1] <- tr(P1 %*% P1s)
  Delta[1, 2] <- tr(P1 %*% P2s)
  Delta[2, 1] <- tr(P2 %*% P1s)
  Delta[2, 2] <- tr(P2 %*% P2s)
  V <- matrix(0, nrow = (k_x + 2), ncol = (k_x + 2))
  V[1:2, 1:2] <- Delta
  V[3:(k_x + 2), 3:(k_x + 2)] <- V22
  V <- sigma2^2 * V
  # Construct first part of Omega
  omega   <- cbind(diag(P1), diag(P2)) # n * 2
  mu4.hat <- sum(epsi^4) / n
  mu3.hat <- sum(epsi^3) / n
  Vp1 <- matrix(0, nrow = (k_x + 2), ncol = (k_x + 2))
  Vp1[1:2, 1:2]         <- (mu4.hat - 3 * sigma2^2) * crossprod(omega)   # 2 * 2 
  Vp1[1:2, 3:(k_x + 2)] <- mu3.hat * crossprod(omega, H)
  Vp1[3:(k_x + 2), 1:2] <- mu3.hat * crossprod(H, omega)
  Omega                 <- Vp1 + V
  Omega                 <- Omega / n
}


## ----slm.gmm-----------------------------------------------------------------------------------------------------------------------------
# Function to estimate the SLM using GMME or OGMME
slm.gmm <- function(formula, data, W, instruments = 2, 
                    estimator = c("gmm", "ogmm"), 
                    gradient = TRUE){
  # Model Frame
  callT    <- match.call(expand.dots = TRUE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf       <- eval(mf, parent.frame())
  
  # Estimator
  estimator <- match.arg(estimator)
  
  # Get variables and globals
  y  <- model.response(mf)
  X  <- model.matrix(formula, mf)
  n  <- nrow(X)
  Wy <- W %*% y
  sn <- nrow(W)
  if (n != sn) stop("number of spatial units in W is different to the number of data")
  
  # Linear instruments
  H <- cbind(X, make.H(W = W, X = X, l = instruments))
  
  # Starting values for optimization
  start <- c(cor(Wy, y), coef(lm(y ~ X - 1)))
  names(start) <- c("Wy", colnames(X))
  
  # GMM estimator with weighting matrix using a identity matrix  
  k_x <- ncol(H)
  Psi <- diag(2 + k_x)
  require("maxLik")
  opt <- maxLik(logLik = Qmin, 
                start = start, 
                method = "bfgs", 
                y = y, 
                X = X, 
                H = H, 
                W = W, 
                Psi = Psi, 
                gradient = gradient,
                print.level = 3, 
                finalHessian = FALSE)
  
  # OGMM: GMM estimator with weighting matrix using the inverse of the var-cov of moment functions
  if (estimator == "ogmm"){
    # Compute Omega.hat/n
    Omega.hat <- make.vmom(coef(opt), y = y, X = X, H = H, W = W)
    Psi  <- solve(Omega.hat)
    opt <- maxLik(logLik = Qmin, 
                  start = coef(opt), 
                  method = "bfgs", 
                  y = y, 
                  X = X, 
                  H = H, 
                  W = W, 
                  Psi = Psi, 
                  gradient = gradient,
                  print.level = 3, 
                  finalHessian = FALSE)
  }
  
  
  results <- structure(
    list(
      coefficients = coef(opt),
      call = callT, 
      X = X, 
      H = H, 
      y = y, 
      Psi = Psi, 
      W = W, 
      estimator = estimator
    ),
    class = "gmm.slm"
  )
  return(results)
}


## ----make.D------------------------------------------------------------------------------------------------------------------------------
# Create D matrix for asymptotic distribution
make.D <- function(rho, beta, y, X, H, W){
  n      <- nrow(X)
  k      <- ncol(X)
  k_x    <- ncol(H)
  I      <- diag(n)
  S      <- I -  rho * W
  epsi   <- crossprod(t(S), y) - crossprod(t(X), beta)
  sigma2 <- as.numeric(crossprod(epsi) / n)
  P1     <- W
  #W2     <- crossprod(t(W), W)
  W2     <- tcrossprod(W, W)
  P2     <- W2 - (tr(W2) / n) * diag(n)
  P1s    <- P1 + t(P1)
  P2s    <- P2 + t(P2)
  G      <- W %*% solve(S)
  
  # Gen D
  D <- matrix(0, nrow = (k_x + 2) , ncol = k + 1)
  rownames(D) <- c("q1", "q2", colnames(H))
  colnames(D) <- c("Wy", colnames(X))
  D[1, 1] <- sigma2 * tr(P1s %*% G)
  D[2, 1] <- sigma2 * tr(P2s %*% G)
  D[3:(k_x + 2), 1] <- t(H) %*% (G %*% X %*% beta) 
  D[3:(k_x + 2), 2:(k + 1)]  <- t(H) %*% X
  return(D)
}


## ----vcov.gmm.slm------------------------------------------------------------------------------------------------------------------------
# Variance-covariance matrix
vcov.gmm.slm <- function(object, D = c("population", "gradient"), ...){
  estimator <- object$estimator
  D.type <- match.arg(D)
  X      <- object$X
  H      <- object$H
  y      <- object$y
  W      <- object$W
  k      <- ncol(X)
  n      <- nrow(X)
  b.hat  <- object$coefficients
  rho    <- b.hat[1L]
  beta   <- b.hat[(2L:(k + 1))]
  # Matrix D is (k_x + 2) * (k + 1)
  if (estimator == "gmm"){
    if (D.type == "population"){
      D     <- make.D(rho = rho, beta = beta, y = y, X = X, H = H, W = W) / n
    }
    if (D.type == "gradient"){
      D     <- moments.lee2007(b.hat, y = y, X = X, H = H, W = W)$D 
    }
    Omega <- make.vmom(b.hat = b.hat, y = y, X = X, H = H, W = W)
    var   <- solve(crossprod(D)) %*% t(D) %*% Omega %*% D %*% solve(crossprod(D)) / n
  }
  if (estimator == "ogmm"){
    if (D.type == "population"){
      D     <- make.D(rho = rho, beta = beta, y = y, X = X, H = H, W = W) /n
    }
    if (D.type == "gradient"){
      D     <- moments.lee2007(b.hat, y = y, X = X, H = H, W = W)$D
    }
    Psi <- object$Psi
    var <- solve(t(D) %*% Psi %*% D) / n
  }
  return(var)
}


## ----summary.gmm.slm---------------------------------------------------------------------------------------------------------------------
summary.gmm.slm <- function(object, 
                            D = c("population", "gradient"), 
                             table = TRUE, 
                             digits = max(3, .Options$digits - 3),
                             ...){
  D.type <- match.arg(D)
  n       <- nrow(object$X)
  df      <- n - length(object$coefficients)
  b       <- object$coefficients
  std.err <- sqrt(diag(vcov(object, D = D.type)))
  z       <- b / std.err
  p       <- 2 * pt(-abs(z), df = df)
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std.Error", "t-value", "Pr(>|t|)")
  result <- structure(
    list(
      CoefTable = CoefTable, 
      digits    = digits, 
      call      = object$call),
    class = 'summary.gmm.slm'
  )
  return(result)
}

print.summary.gmm.slm <- function(x, 
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


## ----test-gmm-slm, cache = TRUE----------------------------------------------------------------------------------------------------------
# Test 
gmm <- slm.gmm(y ~ x1 + x2 + x3, data = data, instruments = 2, W = W)

summary(gmm, D = "population")
summary(gmm, D = "gradient")

ogmm <- slm.gmm(y ~ x1 + x2 + x3, data = data, instruments = 2, W = W, 
                estimator = "ogmm")
summary(ogmm, D = "population")
summary(ogmm, D = "gradient")


## ----echo = FALSE, results = 'asis', warning=FALSE---------------------------------------------------------------------------------------
mle <- sar.mle.con(y ~ x1 + x2 + x3, data = data, W = W)
coefs <- cbind(c(mle$beta_hat, mle$rho_hat),
               s2sls.e$coefficients, 
               b2sls$coefficients,
               gmm$coefficients[c(2, 3, 4, 5, 1)], 
               ogmm$coefficients[c(2, 3, 4, 5, 1)]
           )
Tab <- cbind(coefs)
rownames(Tab) <- c("b0", "b1", "b2", "b3", "rho")
colnames(Tab) <- c("mle", "s2sls", "bs2sls", "gmm", "ogmm")
toLatex(Tab, compact = TRUE, useBooktabs =  TRUE, digits = 4)


## ----echo = FALSE, results = 'asis', warning=FALSE---------------------------------------------------------------------------------------
ses <- cbind(sqrt(diag(vcov(mle)))[c(1, 2, 3, 4, 6)],
             sqrt(diag(vcov(s2sls.e))), 
             sqrt(diag(vcov(s2sls.e, tse = "rob"))),
             sqrt(diag(vcov(b2sls, estimate = "initial"))), 
             sqrt(diag(vcov(b2sls, estimate = "final"))), 
             sqrt(diag(vcov(gmm)))[c(2, 3, 4, 5, 1)], 
             sqrt(diag(vcov(gmm, D = "gradient")))[c(2, 3, 4, 5, 1)], 
              sqrt(diag(vcov(ogmm)))[c(2, 3, 4, 5, 1)], 
             sqrt(diag(vcov(ogmm, D = "gradient")))[c(2, 3, 4, 5, 1)]
             )
Tab <- cbind(ses)
rownames(Tab) <- c("b0", "b1", "b2", "b3", "rho")
colnames(Tab) <- c("mle", "s2sls-ho", "s2sls-he", "bs2sls-i", "bs2sls-f", "gmm-p", "gmm-g", "ogmm-p", "ogmm-g")
toLatex(Tab, compact = TRUE, useBooktabs =  TRUE, digits = 4)


## ----mom.sem-----------------------------------------------------------------------------------------------------------------------------
# Function that generates g and G
mom.sem <- function(u, M){
  # This function generates the moment conditions
  # inputs: Consistent residuals and W matrix
  n       <- length(u)
  u_l     <- W %*% u
  u_ll    <- W %*% u_l
  trMM    <- sum(diag(crossprod(M)))
  uu      <- crossprod(u)
  uul     <- crossprod(u, u_l)
  uull    <- crossprod(u, u_ll)
  ullul   <- crossprod(u_ll, u_l)
  ullull  <- crossprod(u_ll, u_ll)
  ulul    <- crossprod(u_l, u_l)
  ulull   <- crossprod(u_l, u_ll)
  G       <- matrix(0, 3, 3)
  G[1, 1] <- 2 * uul
  G[2, 1] <- 2 * ullul
  G[3, 1] <- uull + ulul
  G[1, 2] <- -ulul 
  G[2, 2] <- -ullull
  G[3, 2] <- -ulull 
  G[1, 3] <- 1 * n
  G[2, 3] <- trMM 
  G <- G / n
  g <- c(uu, ulul, uul) / n
  list(G = G, g = g)
}  


## ----Qn----------------------------------------------------------------------------------------------------------------------------------
# Function to be optimized
Qn <- function(par, mom, verbose = verbose){
    # par has lambda and sigma
    upsilon <- mom$g - mom$G %*% c(par[1], par[1]^2, par[2])
    upup <- crossprod(upsilon)
    if (verbose) 
      cat("function:", upup, "lambda:", par[1], "sig2:", 
          par[2], "\n")
   return(upup)
}


## ----sem.sfgls---------------------------------------------------------------------------------------------------------------------------
# Main function for the FSGLSE
sem.sfgls <- function(formula, data, M, verbose = TRUE){
  # Model Frame
  callT    <- match.call(expand.dots = TRUE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  # Get variables and globals
  y  <- model.response(mf)
  X  <- model.matrix(formula, mf)
  n  <- nrow(X)
  k  <- ncol(X)
  sn <- nrow(M)
  if (n != sn) stop("number of spatial units in W is different to the number of data") 
  
  # First Step: Obtain consistent residuals from OLS
  ols   <- lm(y ~ X - 1)
  u.hat <- residuals(ols)
  
  # Generate Moments
  mom.hat <- mom.sem(u = u.hat, M = M)
  
  # Initial values for lambda and sigma
  scorr <- crossprod(W %*% u.hat, u.hat) / crossprod(u.hat)
  par  <- c(scorr, var(u.hat))
  
  #Optimization
  opt <- nlminb(start =  par, Qn, mom = mom.hat, verbose = verbose)
  
  # B FSGLS
  lambda.hat <- opt$par[1L]
  ys      <- y - drop(lambda.hat) * W %*% y
  Xs      <- X - drop(lambda.hat) * W %*% X
  b.hat   <- solve(crossprod(Xs)) %*% crossprod(Xs, ys)
  
  # Residuals
  e.hat <- ys - Xs %*% b.hat
  
  # Save results
  results <- structure(
    list(
      coefficients = c(b.hat, lambda.hat), 
      call         = callT, 
      X            = X, 
      y            = y, 
      Xs           = Xs, 
      e.hat        = e.hat
    ), 
    class = 'myfs2sls'
  )
}


## ----generate-DGP-sem1-------------------------------------------------------------------------------------------------------------------
# Generate DGP
library("spatialreg")
library("spdep")
set.seed(1)
n      <- 529
lambda <- 0.6
M.nb2  <- cell2nb(sqrt(n), sqrt(n))
M      <- nb2mat(M.nb2)

# Exogenous variables
x1     <- rnorm(n)
x2     <- rnorm(n)
x3     <- rnorm(n)

# DGP parameters
b0 <- 0 ; b1 <- -1; b2 <- 0; b3 <- 1
sigma2 <- 2
epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma2))

# Simulate the dependent variable
y <- b0 + b1*x1 + b2*x2 + b3*x3 + solve(diag(n) -  lambda * M) %*% epsilon

data <- as.data.frame(cbind(y, x1, x2, x3))
names(data) <- c("y", "x1", "x2", "x3")


## ----------------------------------------------------------------------------------------------------------------------------------------
vcov.myfs2sls <- function(object, ...){
  sigma2 <- crossprod(object$e.hat) / n
  var    <- drop(sigma2) * solve(crossprod(object$Xs))
  return(var)
}
summary.myfs2sls <- function(object, 
                              table = TRUE, 
                              digits = max(3, .Options$digits - 3), 
                              ...){
    n       <- nrow(object$X)
    K       <- ncol(object$X)
    df      <- n - K
    b       <- object$coefficients[1:K]
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
      class = 'summary.myfs2sls'
    )
    return(result)
}

print.summary.myfs2sls <- function(x, 
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


## ----b.fs2sls----------------------------------------------------------------------------------------------------------------------------
b.fs2sls <- sem.sfgls(y ~ x1 + x2 + x3, data = data, M = M, verbose = FALSE)
summary(b.fs2sls)


## ----------------------------------------------------------------------------------------------------------------------------------------
sem_mm    <- GMerrorsar(y ~ x1 + x2 + x3, 
                        data = data,
                        listw = mat2listw(M, style = "W"),
                        verbose = FALSE,
                        legacy =  TRUE)
summary(sem_mm)


## ----moments.lee2010---------------------------------------------------------------------------------------------------------------------
moments.lee2010 <- function(theta, y, X, M){
  # Theta: beta and lambda
  k        <- ncol(X)
  n        <- nrow(X)
  beta     <- theta[1:k]
  lambda   <- tail(theta, n = 1L)
  I        <- diag(n)
  R        <- I -  lambda * M
  u.n      <- y - crossprod(t(X), beta)
  epsi     <- R %*% u.n
  P1       <- M
  M2       <- crossprod(t(M), M)
  P2       <- M2 - (tr(M2) / n) * diag(n)
  g.lin    <- crossprod(X, epsi)            # k * 1
  g.q1     <- crossprod(epsi, P1) %*% epsi  # 1*1
  g.q2     <- crossprod(epsi, P2) %*% epsi  # 1*1
  g        <- rbind(g.lin, g.q1, g.q2)
  
  # Gradient 
  P1s <- P1 + t(P1)
  P2s <- P2 + t(P2)
  D   <-  -1 * rbind(t(X), 
                     crossprod(epsi, P1s), 
                     crossprod(epsi, P2s) 
                     ) %*% cbind(R %*% X, M %*% u.n)
  
  # Return results (note that they are divided by n)
  out <- list(g = g /n , D = D / n)
  return(out)
}


## ----Q.sem-------------------------------------------------------------------------------------------------------------------------------
Q.sem <- function(start, y, X, M, Psi, gradient){
  g.hat <- moments.lee2010(theta = start, y = y, X = X, M = M)
  Q <- -1 * crossprod(g.hat$g, Psi) %*% g.hat$g
  if (gradient){
    D <- g.hat$D
    Gr <- -2 * crossprod(D, Psi) %*% g.hat$g
    attr(Q, "gradient") <- as.vector(Gr)
  }
  return(Q)
}


## ----make.vmom.sem-----------------------------------------------------------------------------------------------------------------------
make.vmom.sem <- function(b.hat, y, X, M){
  k        <- ncol(X)
  n        <- nrow(X)
  beta     <- b.hat[1:k]
  lambda   <- tail(b.hat, n = 1L)
  I        <- diag(n)
  R        <- I -  lambda * M
  u.n      <- y - crossprod(t(X), beta)
  epsi     <- R %*% u.n
  sigma2   <- as.numeric(crossprod(epsi) / n)
  P1       <- M
  M2       <- crossprod(t(M), M)
  P2       <- M2 - (tr(M2) / n) * diag(n)
  P1s      <- P1 + t(P1)
  P2s      <- P2 + t(P2)
  
  # Construct V: (2 + k) * (2 + k)
  V11 <- (1 / sigma2) * crossprod(X)         # k + k
  Delta <- matrix(0, nrow = 2, ncol = 2)
  Delta[1, 1] <- tr(P1s %*% P1)
  Delta[1, 2] <- tr(P1s %*% P2)
  Delta[2, 1] <- tr(P2s %*% P1)
  Delta[2, 2] <- tr(P2s %*% P2)
  V <- matrix(0, nrow = (k + 2), ncol = (k + 2))
  V[1:k, 1:k] <- V11
  V[(k + 1):(k + 2), (k + 1):(k + 2)] <- Delta
  V <- sigma2^2 * V
  # Construct first part of Omega
  omega   <- cbind(diag(P1), diag(P2)) # n * 2
  mu4.hat <- sum(epsi^4) / n
  mu3.hat <- sum(epsi^3) / n
  Vp1 <- matrix(0, nrow = (k + 2), ncol = (k + 2))
  Vp1[(k + 1):(k + 2), (k + 1):(k + 2)] <- (mu4.hat - 3 * sigma2^2) * crossprod(omega)   # 2 * 2 
  Vp1[(k + 1):(k + 2), 1:k] <- mu3.hat * crossprod(omega, X)
  Vp1[1:k, (k + 1):(k + 2)] <- mu3.hat * crossprod(X, omega)
  Omega                 <- Vp1 + V
  Omega                 <- Omega / n
}


## ----sem.gmm-----------------------------------------------------------------------------------------------------------------------------
sem.gmm <- function(formula, data, M, 
                    estimator = c("gmm", "ogmm"), 
                    gradient = TRUE){
  # Model Frame
  callT    <- match.call(expand.dots = TRUE)
  mf       <- callT
  m        <- match(c("formula", "data"), names(mf), 0L)
  mf       <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf       <- eval(mf, parent.frame())
  
  # Estimator
  estimator <- match.arg(estimator)
  
  # Get variables and globals
  y  <- model.response(mf)
  X  <- model.matrix(formula, mf)
  n  <- nrow(X)
  sn <- nrow(M)
  if (n != sn) stop("number of spatial units in W is different to the number of data")
  
  # Starting values for optimization
  ols.e <- lm(y ~ X - 1)
  ols.r <- residuals(ols.e)
  start <- c(coef(ols.e), cor(M %*% ols.r, ols.r))
  names(start) <- c(colnames(X), "Mu")
  
  # GMM estimator with weighting matrix using a identity matrix  
  k   <- ncol(X)
  Psi <- diag(k + 2)
  require("maxLik")
  opt <- maxLik(logLik = Q.sem, 
                start  = start, 
                method = "bfgs", 
                y      = y, 
                X      = X, 
                M      = M, 
                Psi    = Psi, 
                gradient = gradient,
                print.level = 3, 
                finalHessian = FALSE)
  
  # OGMM: GMM estimator with weighting matrix using the inverse of the var-cov of moment functions
  if (estimator == "ogmm"){
    Omega.hat <- make.vmom.sem(coef(opt), y = y, X = X, M = M)
    Psi       <- solve(Omega.hat)
    opt       <- maxLik(logLik = Q.sem, 
                        start = coef(opt), 
                        method = "bfgs", 
                        y = y, 
                        X = X, 
                        M = M, 
                        Psi = Psi, 
                        gradient = gradient,
                        print.level = 3, 
                        finalHessian = FALSE)
  }
  
  results <- structure(
    list(
      coefficients = coef(opt),
      call         = callT, 
      X            = X, 
      y            = y, 
      Psi          = Psi, 
      M            = M, 
      estimator    = estimator
    ),
    class = "gmm.sem"
  )
  return(results)
}


## ----make.D.sem--------------------------------------------------------------------------------------------------------------------------
make.D.sem <- function(b.hat, y, X, M){
  k        <- ncol(X)
  n        <- nrow(X)
  beta     <- b.hat[1:k]
  lambda   <- tail(b.hat, n = 1L)
  I        <- diag(n)
  R        <- I -  lambda * M
  u.n      <- y - crossprod(t(X), beta)
  epsi     <- R %*% u.n
  sigma2   <- as.numeric(crossprod(epsi) / n)
  P1       <- M
  M2       <- crossprod(t(M), M)
  P2       <- M2 - (tr(M2) / n) * diag(n)
  P1s      <- P1 + t(P1)
  P2s      <- P2 + t(P2)
  Q        <- M %*% solve(R)
  
  # Generate D
  D <- matrix(0, nrow = (k + 2) , ncol = k + 1)
  rownames(D) <- c(colnames(X), "q1", "q2")
  colnames(D) <- c(colnames(X), "Mu")
  D[k + 1, k + 1] <- sigma2 * tr(P1s %*% Q)
  D[k + 2, k + 1] <- sigma2 * tr(P2s %*% Q)
  D[1:k, 1:k]     <- t(X) %*% R %*% X
  return(D)
}


## ----vcov.gmm.sem------------------------------------------------------------------------------------------------------------------------
vcov.gmm.sem <- function(object, D = c("population", "gradient"), ...){
  estimator <- object$estimator
  D.type <- match.arg(D)
  X      <- object$X
  y      <- object$y
  M      <- object$M
  k      <- ncol(X)
  n      <- nrow(X)
  b.hat  <- object$coefficients
  
  if (estimator == "gmm"){
    if (D.type == "population"){
      D <- make.D.sem(b.hat = b.hat, y = y, X = X, M = M)
      D <- D / n
    }
    if (D.type == "gradient"){
      D <- moments.lee2010(b.hat, y = y, X = X, M = M)$D 
    }
    Omega <- make.vmom.sem(b.hat = b.hat, y = y, X = X, M = M)
    var   <- solve(crossprod(D)) %*% t(D) %*% Omega %*% D %*% solve(crossprod(D)) / n
  }
  if (estimator == "ogmm"){
    if (D.type == "population"){
      D <- make.D.sem(b.hat = b.hat, y = y, X = X, M = M)
      D <- D / n
    }
    if (D.type == "gradient"){
      D <- moments.lee2010(b.hat, y = y, X = X, M = M)$D 
    }
    Psi <- object$Psi
    var <- solve(t(D) %*% Psi %*% D) / n
  }
  return(var)
}


## ----summary.gmm.sem---------------------------------------------------------------------------------------------------------------------
summary.gmm.sem <- function(object,
                            D = c("population", "gradient"),
                             table = TRUE, 
                             digits = max(3, .Options$digits - 3),
                             ...){
  D.type <- match.arg(D)
  n       <- nrow(object$X)
  df      <- n - length(object$coefficients)
  b       <- object$coefficients
  std.err <- sqrt(diag(vcov(object, D = D.type)))
  z       <- b / std.err
  p       <- 2 * pt(-abs(z), df = df)
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std.Error", "t-value", "Pr(>|t|)")
  result <- structure(
    list(
      CoefTable = CoefTable, 
      digits    = digits, 
      call      = object$call),
    class = 'summary.gmm.sem'
  )
  return(result)
}

print.summary.gmm.sem <- function(x, 
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


## ----testsemgmm, cache = TRUE------------------------------------------------------------------------------------------------------------
semgmm <- sem.gmm(y ~ x1 + x2 + x3, data = data, 
                  M = M, estimator = "gmm")
summary(semgmm, D = "population")
summary(semgmm, D = "gradient")

semOgmm <- sem.gmm(y ~ x1 + x2 + x3, data = data, 
                  M = M, estimator = "ogmm")
summary(semOgmm, D = "population")
summary(semOgmm, D = "gradient")

sem_mm    <- GMerrorsar(y ~ x1 + x2 + x3, 
                        data = data,
                        listw = mat2listw(M, style = "W"),
                        verbose = FALSE,
                        se.lambda = TRUE)
summary(sem_mm)

library("sphet")
sem.spreg    <- spreg(y ~ x1 + x2 + x3, 
                        data = data,
                        listw = mat2listw(M, style = "W"),
                        model = "error")
summary(sem.spreg)


## ----dgp-sac-----------------------------------------------------------------------------------------------------------------------------
# Generate Data Generating Process
set.seed(666)
n      <- 529
rho    <- 0.6
lambda <- 0.6
W.nb2  <- cell2nb(sqrt(n), sqrt(n))
W      <- nb2mat(W.nb2)
M.nb2  <- cell2nb(sqrt(n), sqrt(n), type = "queen")
M      <- nb2mat(M.nb2)

# Exogenous variables
x1     <- rnorm(n)
x2     <- rnorm(n)
x3     <- rnorm(n)

# DGP parameters
b0      <- 0 ; b1 <- -1; b2 <- 0; b3 <- 1
sigma2  <- 1
epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma2))

# Simulate the dependent variable
u <- solve(diag(n) -  lambda * M) %*% epsilon
y <- solve(diag(n) -  rho * W) %*% (b0 + b1 * x1 + b2 * x2 + b3 * x3 + u)

# Data as data.frame
data <- as.data.frame(cbind(y, x1, x2, x3))
names(data) <- c("y", "x1", "x2", "x3")


## ----moments.kp--------------------------------------------------------------------------------------------------------------------------
# This function provides the KP's 2010 moments conditions
moments.kp <- function(lambda, u, M, A1, A2){
  n    <- length(u)
  I    <- diag(n)
  R    <- I - lambda * M
  epsi <- R %*% u
  m.q1 <- crossprod(epsi, A1) %*% epsi
  m.q2 <- crossprod(epsi, A2) %*% epsi
  m    <- rbind(m.q1, m.q2)
  return(m)
}


## ----Q.min-------------------------------------------------------------------------------------------------------------------------------
# This function generates the function Q(theta) to be minimized.
Q.min <- function(lambda, u, M, A1, A2, Upsilon, verbose){
  m <- moments.kp(lambda = lambda, u = u, M = M, A1 = A1, A2 = A2)
  Q <- crossprod(m, Upsilon) %*% m / n
  if (verbose) cat("function:", Q, "lambda:", lambda, "\n")
  return(Q)
}


## ----make.Upsilon------------------------------------------------------------------------------------------------------------------------
# This function generates the VC matrix for 
make.Upsilon <- function(lambda, HH.i, H, Z, A1, A2, u, M, 
                         step = c("first", "second")){
  step   <- match.arg(step)
  n      <- nrow(H)
  A1s    <- A1 + t(A1)
  A2s    <- A2 + t(A2)
  R      <- diag(n) - lambda * M
  Zs     <- R %*% Z
  us     <- R %*% u
  alpha1 <- - (1 / n) * (crossprod(Zs, A1s) %*% us)
  alpha2 <- - (1 / n) * (crossprod(Zs, A2s) %*% us)
  if (step == "first"){
    # Based on 2SLS residuals
    HZ     <- crossprod(H, Z) / n
    P      <- HH.i %*% HZ %*% solve(t(HZ) %*% HH.i %*% HZ)
    R.inv  <- solve(t(R))
    a1     <- crossprod(R.inv, H) %*% P %*% alpha1
    a2     <- crossprod(R.inv, H) %*% P %*% alpha2
  } else {
    # Based on G2SLS residuals
    HZ     <- crossprod(H, Zs) / n
    P      <- HH.i %*% HZ %*% solve(t(HZ) %*% HH.i %*% HZ)
    a1     <- H %*% P %*% alpha1
    a2     <- H %*% P %*% alpha2
  }
  epsi   <- R %*% u
  Sigma  <- diag(drop(epsi^2))
  Psi.11 <- (1 / (2 * n)) * tr(A1s %*% Sigma %*% A1s %*% Sigma) + 
            (1 / n) * crossprod(a1, Sigma) %*% a1
  Psi.12 <- (1 / (2 * n)) * tr(A1s %*% Sigma %*% A2s %*% Sigma) + 
            (1 / n) * crossprod(a1, Sigma) %*% a2
  Psi.22 <- (1 / (2 * n)) * tr(A2s %*% Sigma %*% A2s %*% Sigma) + 
            (1 / n) * crossprod(a2, Sigma) %*% a2
  Upsi <- matrix(0, nrow = 2, ncol = 2)
  Upsi[1, 1] <- Psi.11
  Upsi[1, 2] <- Upsi[2, 1] <- Psi.12
  Upsi[2, 2] <- Psi.22
  out <- list(Upsi = Upsi, a1 = a1, a2 = a2, epsi = epsi, P = P, Sigma = Sigma)
  return(out)
}


## ----make.G------------------------------------------------------------------------------------------------------------------------------
# This function generates matrix G
make.G <- function(u, M, A1, A2){
  n       <- length(u)
  u.l     <- M %*% u
  G       <- matrix(0, 2, 2)
  G[1, 1] <- 2  * t(u.l) %*% (A1 + t(A1)) %*% u
  G[1, 2] <- -1 * t(u.l) %*% A1 %*% u.l
  G[2, 1] <- t(u.l) %*% (A2 + t(A2)) %*% u
  G[2, 2] <- -1 * t(u.l) %*% A2 %*% u.l
  G <- G / n
  return(G)
}  


## ----gstsls.sac--------------------------------------------------------------------------------------------------------------------------
gstsls.sac <- function(formula, data, W, M = NULL, instruments = 2,
                       verbose = FALSE){
    # Model Frame
    callT    <- match.call(expand.dots = TRUE)
    mf       <- callT
    m        <- match(c("formula", "data"), names(mf), 0L)
    mf       <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf       <- eval(mf, parent.frame())
    
    # Get variables and globals
    y  <- model.response(mf)
    X  <- model.matrix(formula, mf)
    n  <- nrow(X)
    
    # Check spatial weight matrices
    if (!is.null(M)){
      if(identical(M, W)) {
        M <- W
        has.M <- FALSE
      } else {
        has.M <- TRUE
      }
    } else {
      M <- W
      has.M <- FALSE
    }

    #---------------------------- #
    # Step 1a: S2SLS estimator
    #---------------------------- #
    # Make instruments
    WXI <- make.H(W = W, X = X, l = instruments)
    H   <- cbind(X, WXI)
    if(has.M){
      H  <- cbind(H, M %*% X[, -1], M %*% WXI)
    }
  
    # Compute S2SLSE
    Wy          <- W %*% y
    Z           <- cbind(X, Wy)
    colnames(Z) <- c(colnames(X), "Wy")
    HH          <- crossprod(H)
    HH.i        <- solve(HH)
    PH          <- H %*% tcrossprod(HH.i, H)
    Z_hat       <- PH %*% Z
    b_2sls      <- solve(crossprod(Z_hat)) %*% crossprod(Z_hat, y)
    y_hat       <- Z %*% b_2sls
    u_2sls      <- y - y_hat
    
    #---------------------------------------- #
    # Step 1b: Initial GMM estimator of lambda
    #---------------------------------------- #
    # Initial value for optimization 
    init.1b <- crossprod(M %*% u_2sls, u_2sls) / crossprod(u_2sls)
    
    # Optimization
    Upsilon.1b <- diag(2)
    MM         <- crossprod(M)
    A1.1b      <- MM - diag(diag(MM))
    A2.1b      <- M
    opt.1b <- nlminb(start     = init.1b, 
                     objective = Q.min,
                     u         = u_2sls, 
                     M         = M, 
                     A1        = A1.1b, 
                     A2        = A2.1b, 
                     Upsilon   = Upsilon.1b, 
                     verbose   = verbose,
                     lower     = -0.9 + .Machine$double.eps, 
                     upper      = 0.9 - .Machine$double.eps)
    
    #------------------------------------------- #
    # Step 1c: Efficient GMM estimator of lambda
    #------------------------------------------- #
    
    # Initial value for optimization 
    init.1c <- opt.1b$par
    
    # Optimization
    Upsilon.1c <- make.Upsilon(lambda = init.1c, 
                               HH.i = HH.i,
                               H = H, 
                               Z = Z, 
                               A1 = A1.1b, 
                               A2 = A2.1b, 
                               u = u_2sls, 
                               M = M, 
                               step = "first")$Upsi
    opt.1c <- nlminb(start     = init.1c, 
                     objective = Q.min,
                     u         = u_2sls, 
                     M         = M, 
                     A1        = A1.1b, 
                     A2        = A2.1b, 
                     Upsilon   = solve(Upsilon.1c), 
                     verbose   = verbose,
                     lower     = -0.9 + .Machine$double.eps, 
                     upper      = 0.9 -  .Machine$double.eps)
    
    #---------------------------------------- #
    # Step 2a: GS2SLS
    #---------------------------------------- #
    lambda.hat <- opt.1c$par
    ys         <- y - lambda.hat * M %*% y
    Zs         <- Z - lambda.hat * M %*% Z
    
    # Compute S2SLSE
    Zs_hat  <- PH %*% Zs
    b_g2sls <- solve(crossprod(Zs_hat)) %*% crossprod(Zs_hat, ys)
    y_hat   <- Z %*% b_g2sls
    u_g2sls <- y - y_hat
    
    #---------------------------------------- #
    # Step 2b: Efficient GMM
    #---------------------------------------- #
    # Initial value for optimization 
    init.2b  <- lambda.hat
    
    # Optimization
    Upsilon.2b <- make.Upsilon(lambda = init.2b, 
                               HH.i = HH.i,
                               H = H,
                               Z = Z, 
                               A1 = A1.1b, 
                               A2 = A2.1b, 
                               u = u_g2sls, 
                               M = M,
                               step = "second")$Upsi
    opt.2b <- nlminb(start     = init.2b, 
                     objective = Q.min,
                     u         = u_g2sls, 
                     M         = M, 
                     A1        = A1.1b, 
                     A2        = A2.1b, 
                     Upsilon   = solve(Upsilon.2b), 
                     verbose   = verbose,
                     lower     = -0.9 + .Machine$double.eps, 
                     upper      = 0.9 -  .Machine$double.eps)
    
    b.hat.f <- c(b_g2sls, opt.2b$par)
    names(b.hat.f) <- c(rownames(b_g2sls), "Wu")
    
    # Generate matrix G
    G <- make.G(u = u_g2sls, M = M, A1 = A1.1b, A2 = A2.1b)
    
    # Save results
    results <- structure(
      list(
        coefficients = b.hat.f, 
        call         = callT, 
        X            = X, 
        H            = H,
        HH.i         = HH.i,
        Z            = Z, 
        y            = y,
        W            = W, 
        M            = M, 
        A1           = A1.1b, 
        A2           = A2.1b, 
        G            = G 
      ), 
      class = 'mygs2sls'
    )
}


## ----vcov.mygs2sls-----------------------------------------------------------------------------------------------------------------------
vcov.mygs2sls <- function(object, ...){
  theta.hat  <- object$coefficients
  Z          <- object$Z
  k          <- ncol(Z)
  lambda.hat <- theta.hat[k + 1]
  delta.hat  <- theta.hat[1:k]
  H    <- object$H
  HH.i <- object$HH.i
  A1   <- object$A1
  A2   <- object$A2
  M    <- object$M
  y    <- object$y
  p    <- ncol(H) 
  u_g2sls <- y - Z %*% delta.hat
  fitL <- make.Upsilon(lambda = lambda.hat,
                       HH.i = HH.i,
                       H = H, 
                       Z = Z, 
                       A1 = A1, 
                       A2 = A2, 
                       u = u_g2sls, 
                       M = M,
                       step = "second")
  ## Make Omega
  Psi.rr      <- (1 / n) * crossprod(H, fitL$Sigma) %*% H
  Psi.rd      <- (1 / n) * crossprod(H, fitL$Sigma) %*% cbind(fitL$a1, fitL$a2)
  Psi.f10     <- cbind(Psi.rr, Psi.rd)
  Psi         <- fitL$Upsi 
  Psi.f20     <- cbind(t(Psi.rd), Psi)
  Psi.0       <- rbind(Psi.f10, Psi.f20)
  Psi.inv     <- solve(fitL$Upsi)
  G           <- object$G
  J           <- G %*% rbind(1, 2 * lambda.hat)
  bread.lower <- solve(t(J) %*% Psi.inv %*% J) %*% t(J) %*% Psi.inv
  bread       <- matrix(0, nrow = (k + 1), ncol = (p + 2))
  bread[1:k, 1:p] <- t(fitL$P)
  bread[k + 1, (p + 1):(p + 2)] <- bread.lower
  Omega <- bread %*% Psi.0 %*% t(bread) / n
  return(Omega)
}


## ----summary.mygs2sls--------------------------------------------------------------------------------------------------------------------
summary.mygs2sls <- function(object,
                             table = TRUE, 
                             digits = max(3, .Options$digits - 3), 
                             ...){
  n         <- nrow(object$Z)
  df        <- n - (ncol(object$Z) + 1)
  b         <- object$coefficients
  std.err   <- sqrt(diag(vcov(object)))
  z         <- b / std.err
  p         <- 2 * pt(-abs(z), df = df)
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std.Error", "t-value", "Pr(>|t|)")
  result <- structure(
    list(
      CoefTable = CoefTable, 
      digits    = digits, 
      call      = object$call),
    class = 'summary.mygs2sls'
  )
  return(result)
}

print.summary.mygs2sls <- function(x, 
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


## ----test-gstsls.sac---------------------------------------------------------------------------------------------------------------------
out <- gstsls.sac(y ~ x1 + x2 + x3, 
           data = data, 
           M = M, 
           W = W)
summary(out)


## ----test1-sac---------------------------------------------------------------------------------------------------------------------------
test1 <- spreg(y ~ x1 + x2 + x3, 
               data = data, 
               listw = mat2listw(W, style = "W"),
               listw2 = mat2listw(M, style = "W"),
               model = "sarar", 
               het = TRUE, 
               step1.c = TRUE)
summary(test1)

