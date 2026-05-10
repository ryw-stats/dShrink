setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ar1_cor <- function(p, rho) {
  #AR1 correlation matrix
  exponent <-
    abs(matrix(
      1:p - 1,
      nrow = p,
      ncol = p,
      byrow = TRUE
    ) -
      (1:p - 1))
  rho ^ exponent
}

simone <- function(sim) {
  set.seed(sim)
  
  x <- mvrnorm(m, rep(0, p), Sig1)
  x_aug <- cbind(1, x)
  y <- c(x_aug %*% (beta + b) + rnorm(m, 0, 1))
  z <- x_aug[, 1:(q + 1)]
  Sig_z1 <- t(z) %*% z / m
  M_z1 <- solve(Sig_z1)
  alpha1 <- c(M_z1 %*% t(z) %*% y / m)
  Sig1a <- cov(c(y - z %*% alpha1) * (z %*% M_z1)) / (m - q - 1)
  
  
  x <- mvrnorm(n, rep(0, p), Sig0)
  x_aug <- cbind(1, x)
  y <- c(x_aug %*% beta + rnorm(n, 0, 1))
  
  Sig_x0_aug <- t(x_aug) %*% x_aug / n
  M_x0 <- solve(Sig_x0_aug)
  beta0 <- c(M_x0 %*% t(x_aug) %*% y / n)
  z <- cbind(1, x[, 1:q])
  Sig_z0 <- t(z) %*% z / n
  M_z0 <- solve(Sig_z0)
  alpha0 <- c(M_z0 %*% t(z) %*% y / n)
  Sig0b <- cov(c(y - x_aug %*% beta0) * (x_aug %*% M_x0)) / (n - p - 1)
  Sig0a <- cov(c(y - z %*% alpha0) * (z %*% M_z0)) / (n - q - 1)
  Sig0ba <- cov(c(y - x_aug %*% beta0) * (x_aug %*% M_x0), c(y - z %*% alpha0) * (z %*% M_z0)) / (n - p - 1)

  dat0 <- data.frame(y, x)
  
  pCML <- rep(NA, p + 1)
  try(pCML <- PCML(beta0, dat0, alpha1), silent = T)
  
  H0 <- Sig_x0_aug
  
  t(pCML - beta) %*% H0 %*% (pCML - beta)
}

p <- 30
q <- floor(2 / 3 * p)
n <- 300
m <- 1000

set.seed(0)
sig_x0 <- sqrt(rgamma(p, 2, scale = 0.5))
Sig0 <- diag(sig_x0) %*% ar1_cor(p, 0.4) %*% diag(sig_x0)
sig_x <- sqrt(rgamma(p, 2, scale = 0.5))
Sigx <- diag(sig_x) %*% ar1_cor(p, 0.2) %*% diag(sig_x)
sim <- 500

set.seed(0)
beta <- rnorm(p + 1, 0, sqrt(0.5 / (p + 1))) # runif(p, -0.2, 0.5) #rnorm(p, 0, sqrt(0.5/p)) #sqrt(0.5/p)
A_seq <- seq(0, 1, 0.1)

set.seed(0)
eta <- rbinom(p + 1, 1, 0.5) * rnorm(p + 1, sqrt(1 / (p + 1)), sqrt(0.5 / (p + 1))) #runif(p, - sqrt(3), sqrt(3)) 
Res_seq <- array(NA, c(sim, length(A_seq)))

library(parallel)
library(MASS)


for (i in 1:length(A_seq)) {
  b <-  A_seq[i] * eta
  Sig1 <- (1 - A_seq[i]) * Sig0 + A_seq[i] * Sigx
  
  cl <- makeCluster(60)
  clusterEvalQ(cl, c(library(MASS), library(stats), library(expm), library(abind), 
                     library(matrixStats), source("func_PCML.R")))
  clusterExport(cl, c("beta", "Sig0", "Sig1", "p", "n", "m", "b", "q"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, i] <- Res
}

save(Res_seq, file = paste(paste("Error", p, n, m, "PCML", sep = "-"), ".Rdata", sep = ""))
