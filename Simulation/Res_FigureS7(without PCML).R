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
  y <- c(x_aug %*% (beta + b) + rnorm(m, 0, sqrt(2)))
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
  
  A <- solve(t(z) %*% z) %*% t(z) %*% x_aug
  
  beta_cls <- beta0 - solve(t(x_aug) %*% x_aug) %*% t(A) %*% solve(A %*% solve(t(x_aug) %*% x_aug) %*% t(A)) %*% (A %*% beta0 - alpha1)
  
  eps <- y - x_aug %*% beta0
  w <- c(1 - (q - 1) * sum(eps^2) / (n - p + 1) / (t(x_aug %*% (beta0 - beta_cls)) %*% (x_aug %*% (beta0 - beta_cls))))
  Han <- w * beta0 + (1 - w) * beta_cls
  
  H0 <- Sig_x0_aug
  
  nsample <- m
  
  form1 <- 'y ~ x1'
  for (i in 2:q) {
    form1 <- paste(form1, " + x", i, sep = "")
  }
  
  varname <- "(Intercept)"
  for (i in 1:p) {
    varname <- c(varname, paste("x", i, sep = ""))
  }
  
  model <- list()
  model[[1]] <- list(form = form1, info = data.frame(var = varname[1:(q + 1)], bet = alpha1))
  dat0 <- data.frame(y, x)
  colnames(dat0)[2:(p + 1)] <- varname[-1]
  form0 <- 'y ~ x1'
  for (i in 2:p) {
    form0 <- paste(form0, " + x", i, sep = "")
  }
  
  
  CML <- rep(NA, p + 1)
  try(CML <- fxnCC_LinReg(q + 1, p + 1, y, x[, 1:q], x[, (q + 1):p], alpha1, beta0, n, tol = 1e-6, maxIter = 400, factor = 1)$gammaHat, silent = T)
  
  
  gim <- rep(NA, p + 1)
  try(gim <- gim(form0, 'gaussian', dat0, model, nsample, ref = dat0)$coefficients, silent = T)
  
  Q_meta <- Sig0ba %*% solve(Sig0a + Sig1a)
  meta <- beta0 - Q_meta %*% (alpha0 - alpha1)
  
  
  Sig_d <- Sig0a + Sig1a
  Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
  P_c <- Sig_c
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda <- max(sum(diag(P_c %*% H0)) - 2 * Re(eigen(HR0 %*% P_c %*% HR0)$values[1]), 0)
  r1 <- meta
  v_c <- t(r1) %*% H0 %*% r1
  Prop <- max(0, c(1 - lambda / v_c)) * r1
  
  P_d <- Sig_d
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r2 <- Q_meta %*% (alpha0 - alpha1)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda <- max(sum(diag(P_d %*% H_d)) - 2 * Re(eigen(HR_d %*% P_d %*% HR_d)$values[1]), 0)
  v_d <- t(r2) %*% H0 %*% r2
  Prop <- Prop +  max(0, c(1 - lambda / v_d)) * r2
  
  loss_beta0 <- t(beta0 - beta) %*% H0 %*% (beta0 - beta)
  #loss_meta <- t(meta - beta) %*% H0 %*% (meta - beta)
  loss_Han <- t(Han - beta) %*% H0 %*% (Han - beta)
  loss_gim <- t(gim - beta) %*% H0 %*% (gim - beta)
  loss_CML <- t(CML - beta) %*% H0 %*% (CML - beta)
  loss_Prop <- t(Prop - beta) %*% H0 %*% (Prop - beta)
  c(loss_Prop, loss_beta0, loss_gim, loss_CML, loss_Han) 
}

p <- 10
q <- floor(2 / 3 * p)
n <- 100
m <- 500

set.seed(0)
sig_x0 <- sqrt(rgamma(p, 2, scale = 0.5))
Sig0 <- diag(sig_x0) %*% ar1_cor(p, 0.4) %*% diag(sig_x0)
sig_x <- sqrt(rgamma(p, 2, scale = 0.5))
Sigx <- diag(sig_x) %*% ar1_cor(p, 0.2) %*% diag(sig_x)
sim <- 5000

set.seed(0)
beta <- rnorm(p + 1, 0, sqrt(0.5 / (p + 1))) # runif(p, -0.2, 0.5) #rnorm(p, 0, sqrt(0.5/p)) #sqrt(0.5/p)
A_seq <- seq(0, 1, 0.1)

set.seed(0)
eta <- rbinom(p + 1, 1, 0.5) * rnorm(p + 1, sqrt(1 / (p + 1)), sqrt(0.5 / (p + 1))) #runif(p, - sqrt(3), sqrt(3)) 
Res_seq <- array(NA, c(5, sim, length(A_seq)))

library(parallel)
library(MASS)


for (i in 1:length(A_seq)) {
  b <-  A_seq[i] * eta
  Sig1 <- (1 - A_seq[i]) * Sig0 + A_seq[i] * Sigx
  
  cl <- makeCluster(20)
  clusterEvalQ(cl, c(library(MASS), library(MetaIntegration), library(gim)))
  clusterExport(cl, c("beta", "Sig0", "Sig1", "p", "n", "m", "b", "q"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, , i] <- Res
}

save(Res_seq, file = paste(paste("Error", p, n, m, "Cmpr", "hVar", sep = "-"), ".Rdata", sep = ""))