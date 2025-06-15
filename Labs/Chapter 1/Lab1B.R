##############################################=
# Spatial Econometrics
# Laboratory: More on Moran's I
# By: Mauricio Sarrias
##############################################=

## 0: Clean directory and load packages ====
rm(list = ls(all = TRUE)) 
library("spdep")
set.seed(1) # In order to get the same results 

# 1. Creating Spatial Weights for a Grid Layout ----

# Create a rook type neighbor object for a 4 by 4 square grid layout (16 spatial units)
rook4x4 <- cell2nb(4, 4)
summary(rook4x4)

# 2. Converting Spatial Weights to a Matrix ----
help("nb2mat")
W <- nb2mat(rook4x4, 
            style = "W") #Normalized W matrix
W

# 3. Creating a Spatial Lag Using Matrix Multiplication ----
x <- rnorm(16, mean = 0, sd = 1)
Wx <- W %*% x  # Spatial lag by hand
sWx <- lag.listw(nb2listw(rook4x4), x) # Spatial Lag using command

Wx[1:4] == sWx[1:4]

# 5. Computing Moran's I Using Matrix Algebra ----

# Recall that the Moran's I is
#
# I = z' W z / (z'z) * (n / S_0)
#
z <- x - mean(x)
zz  <- crossprod(z)
Wz  <- lag.listw(nb2listw(rook4x4), z)
zWz <- crossprod(z, Wz) # t(z) %*% Wz
n   <- length(x)
S0  <- sum(as.vector(W))

# See also the following
spweights.constants(nb2listw(rook4x4))

mi1  <- zWz / zz
mi2  <- (n / S0) * zWz / zz
mi1
mi2

summary(lm(Wz ~ z - 1)) # In form of regression. 
coef(lm(Wz ~ z - 1))


# 4. Theoretical Moments of Moran's I ====
EI <- -1 / (n - 1)

wc <- spweights.constants(nb2listw(rook4x4))
S02 <- wc$S0*wc$S0
EI2 <-  (wc$nn*wc$S1 - wc$n*wc$S2 + 3*S02) / (S02*(wc$nn - 1))
VI <- EI2 - EI ^ 2
stat <- (mi1 - EI) / sqrt(VI)
stat

moran.test(x, nb2listw(rook4x4), randomisation = FALSE)

# Try with matrices
S0 <- sum(W)
c <- n / S0
Ws <- ((W + t(W)) /2) 
I <- c * (t(z) %*% Ws %*% z) / (t(z) %*% z)
I

tr <- function(x) sum(diag(x))

P <- diag(n) - (1 / n)* matrix(1, nrow = n, ncol = n)
A <- P %*% Ws %*% P
B <- P
A_s <- (A + t(A)) /2
B_s <- (B + t(B)) / 2
I2 <- c^2 * (sum(diag(A))^2 + sum(diag(A %*% A_s))) / (sum(diag(B))^2 + sum(diag(B %*% B_s))) 

I2 <- c^2 * ((tr(A)^2 + tr(A %*% t(A)) + tr(A %*% A)) / (tr(B)^2 + tr(B %*% t(B)) + tr(B %*% B)))
I2
EI2

I2 <- c^2 * (tr(A)^2 + tr(A %*% A_s)) / (tr(B)^2 + tr(B %*% B_s))
EI2


# 5. My own function under normality ====

my_moran_nor <- function(x, listw){
  # x     : a vector
  # listw : a list of neighboors of class "listw"
  
  require("spdep")        # just to be sure that spdep is loaded when required
  n <- length(listw$neighbours) # number of spatial units
  if (n != length(x)) stop("x and W does not have the same number of units")
  
  wc <- spweights.constants(listw) # Important constants from W
  S0  <- wc$S0
  S02 <- (S0) ^ 2
  
  # Compute Morans' I
  z   <- x - mean(x)
  zz  <- crossprod(z)
  Wz  <- lag.listw(listw, z)
  zWz <- crossprod(z, Wz)
  mI  <- (n / S0) * (zWz / zz)
  
  # Theoretical moments under normality
  EI <- (-1) / wc$n1
  EI2 <- (wc$nn * wc$S1 - wc$n * wc$S2 + 3 * S02) / (S02 * (wc$nn - 1))
  VI <- EI2 - EI ^ 2
  out <- list("I" = mI, "EI" = EI, "VI" = VI)
  out
}

my_moran_nor(x = x, listw = nb2listw(rook4x4))
moran.test(x, nb2listw(rook4x4), randomisation = FALSE, alternative = "two.sided")

# Two side test
moran_out <- my_moran_nor(x, nb2listw(rook4x4))
2 * pnorm(abs((moran_out$I - moran_out$EI) / sqrt(moran_out$VI)), lower.tail = FALSE) 


# 7. Monte Carlo Moran's I ----

S <- 999
res <- numeric(length = S + 1)
for (i in 1:S) {
  # sample takes a sample of the specified size from the elements of x
  # using either with or without replacement
    res[i] <- my_moran_nor(sample(x), nb2listw(rook4x4))$I
}
res[S + 1] <- my_moran_nor(x, nb2listw(rook4x4))$I
statistic <- res[S + 1]

hist(res[1:S], xlab = "I", col = "grey", main  = "Distribution of Morans' I")
abline(v = res[S + 1], col = "red", lwd = 2)
p <- (1 + sum(res[1:S] >= res[S + 1])) / (1 + S)
p

moran.mc(x, nb2listw(rook4x4), nsim = 999)

# other example,
rook5x10 <- cell2nb(5, 10)
y <- rnorm(5 * 10)
Iw <- invIrM(rook5x10, 0.9)
yc <- c(Iw %*% y)

res2 <- numeric(length = S + 1)
for (i in 1:S) {
  # sample takes a sample of the specified size from the elemenets of x
  # using either with or withour replacement
  res2[i] <- my_moran_nor(sample(yc), nb2listw(rook5x10))$I
}
res2[S + 1] <- my_moran_nor(yc, nb2listw(rook5x10))$I
statistic2 <- res2[S + 1]

hist(res2[1:S], xlab = "I", col = "grey", main  = "Distribution of Morans' I", xlim = c(-0.3, 1))
abline(v = res2[S + 1], col = "red", lwd = 2)
p <- (1 + sum(res2[1:S] >= res2[S + 1])) / (1 + S)
p

moran.mc(yc, nb2listw(rook5x10), nsim = 999)

moran.plot(yc,nb2listw(rook5x10), pch = 19,
           xlab = expression(x),
           ylab = expression(Wx),
           col = "red")
text(c(3,-3, -3, 3), c(-0.5, -0.5, -2, -3),  labels =c("Quadrant I", "Quadrant II",
                                                "Quadrant III", "Quadrant IV"))


# 8. Simulation Morans' I
sim_moran_nor <- function(nrow, ncol, S, pvalue = 0.05){
              morans <- vector(mode = "numeric", length = S)
              reject <- vector(mode = "numeric", length = S)
              W <- nb2listw(cell2nb(nrow, ncol))
              n <- nrow * ncol
              for (s in 1:S) {
                y <- rnorm(n)
                mi <- moran.test(y, W, alternative = "two.sided", randomisation = FALSE)
                morans[s] <- mi$statistic
                if (mi$p.value < pvalue) reject[s] <- 1 else reject[s] <- 0
              }
              rejfreq <- sum(reject) / S
              out <- list(morans = morans, rej = rejfreq)
}

m_1000 <- sim_moran_nor(5, 5, 1000)
m_1000$morans
m_1000$rej

plot(density(m_1000$morans), main = "Moran's I under Null-NORMAL", lwd = 2, col = 2)

# More

m_10    <- sim_moran_nor(5, 2, 1000)
m_25    <- sim_moran_nor(5, 5, 1000)
m_100   <- sim_moran_nor(10, 10, 1000)
m_1000  <- sim_moran_nor(10, 100, 1000)

plot(density(m_10$morans), main = "Moran's I under Null-NORMAL", lwd = 2, col = 1,
     xlim = c(-5, 5), ylim = c(0, 0.45))
lines(density(m_25$morans),  lwd = 2, col = 2)
lines(density(m_100$morans), lwd = 2, col = 3)
lines(density(m_1000$morans), lwd = 2, col = 4)

