getSummary.Sarlm <- function(obj, alpha = 0.05, ...){
  smry <- summary(obj)
  coef <- smry$Coef
  rho <- c(smry$rho, smry$rho.se, smry$rho/smry$rho.se, 
                2 * (1 - pnorm(abs(smry$rho/smry$rho.se))))
  lambda <- c(smry$lambda, smry$lambda.se, smry$lambda/smry$lambda.se, 
              2 * (1 - pnorm(abs(smry$lambda/smry$lambda.se))))
  coef <- rbind(coef, rho, lambda)
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha/2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha/2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  length(obj$residuals)
  sumstat <- c(logLik = logLik(obj), deviance = deviance(obj), AIC = AIC(obj), BIC = BIC(obj), N = N, 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
               Nagelkerke = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}

getSummary.Stsls <- function(obj, alpha = 0.05, ...){
  smry <- summary(obj)
  coef <- smry$Coef
  rownames(coef)[1] <- "rho"
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha/2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha/2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  length(obj$residuals)
  sumstat <- c(logLik = NA, deviance = NA, AIC = NA, BIC = NA, N = N, 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
               Nagelkerke = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}

getSummary.Gmsar <- function(obj, alpha = 0.05, ...){
  smry <- summary(obj)
  coef <- smry$Coef
  lambda <- c(obj$lambda, obj$lambda.se, obj$lambda/obj$lambda.se, 
           2 * (1 - pnorm(abs(obj$lambda/obj$lambda.se))))
  coef <- rbind(coef, lambda)
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha/2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha/2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  length(obj$residuals)
  sumstat <- c(logLik = NA, deviance = NA, AIC = NA, BIC = NA, N = N, 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
               Nagelkerke = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}