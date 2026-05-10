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
  
  x1 <- mvrnorm(m, rep(0, p), Sig1)
  x_aug <- cbind(1, x1)
  y1 <- c(x_aug %*% (mu + b) + rnorm(m, 0, 1))
  Sig_x1_aug <- t(x_aug) %*% x_aug / m
  M_x1 <- solve(Sig_x1_aug)
  alpha1 <- c(M_x1 %*% t(x_aug) %*% y1 / m)
  Sig1a <- cov(c(y1 - x_aug %*% alpha1) * (x_aug %*% M_x1)) / (m - p - 1)
  
  beta <- mu - b
  x0 <- mvrnorm(n, rep(0, p), Sig0)
  x_aug <- cbind(1, x0)
  y0 <- c(x_aug %*% beta + rnorm(n, 0, 1))
  
  Sig_x0_aug <- t(x_aug) %*% x_aug / n
  M_x0 <- solve(Sig_x0_aug)
  beta0 <- c(M_x0 %*% t(x_aug) %*% y0 / n)
  alpha0 <- beta0
  Sig0b <- cov(c(y0 - x_aug %*% beta0) * (x_aug %*% M_x0)) / (n - p - 1)
  Sig0a <- Sig0b
  Sig0ba <- Sig0b
  
  Q_meta <- Sig0ba %*% solve(Sig0a + Sig1a)
  meta <- beta0 - Q_meta %*% (alpha0 - alpha1)
  
  H0 <- Sig_x0_aug

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
  loss_Prop <- t(Prop - beta) %*% H0 %*% (Prop - beta)
  
  n.vec <- c(n, m)
  I.til <- 1:floor(n / 3)
  
  X <- rbind(scale(x0, scale = F), scale(x1, scale = F))
  y <- c(scale(y0, scale = F), scale(y1, scale = F))
  TL <- try(Trans.lasso(X, y, n.vec, I.til, l1=T)$beta.hat, silent = T)
  if (mode(TL) == "character") {
    TLasso <- c(mean(y0), rep(0, p))
  } else {
    TLasso <- c(mean(y0), TL)
  }
  TGLM <- c(glmtrans(target = list(x = x0, y = y0), source = list(list(x = x1, y = y1)), family = "gaussian")$beta)
  
  
  loss_TLasso <- t(TLasso - beta) %*% H0 %*% (TLasso - beta)
  loss_TGLM <- t(TGLM - beta) %*% H0 %*% (TGLM - beta)
  
  lambda_seq <- seq(0,7,length.out=10)
  ISEDI <- LDW_meta(cbind(1, x1), cbind(1, x0), y1, y0, lambda_seq, family="gaussian")$estimates_opt[, 2]
  
  loss_ISEDI <- t(ISEDI - beta) %*% H0 %*% (ISEDI - beta)
  
  x_aug <- cbind(1, x1)
  V1 <- t(x_aug) %*% x_aug
  
  s2_1 <- sum((y1 - x_aug %*% alpha1)^2) / (m - p - 1)
  
  ev1 <- eigen(V1)
  
  x_aug <- cbind(1, x0)
  V0 <- t(x_aug) %*% x_aug
  
  s2_0 <- sum((y0 - x_aug %*% beta0)^2) / (n - p - 1)
  
  ev0 <- eigen(V0)
  
  V0_sqrt <- ev0$vectors %*% diag(sqrt(ev0$values)) %*% t(ev0$vectors)
  
  evM <- eigen(V0_sqrt %*% solve(V1) %*% V0_sqrt)
  
  gamma <- alpha1 - beta0
  
  ev_tmp <- eigen(gamma %*% t(gamma) - s2_1 * solve(V1) - s2_0 *  solve(V0))
  
  Theta <- s2_1 * solve(V1) + ev_tmp$vectors %*% diag(pmax(0, ev_tmp$values)) %*% t(ev_tmp$vectors)
  
  kappa2 <- diag(t(evM$vectors) %*% V0_sqrt %*% Theta %*% V0_sqrt %*% evM$vectors)
  nu <- evM$values
  
  mse_DELR <- function(lambda) {
    sum((s2_1 * (1 + lambda * nu)^2 + lambda^2 * kappa2) / (1 + lambda + lambda * nu)^2)
  }
  
  lambda_DELR <- nlminb(0, mse_DELR, lower = 0)$par
  
  W_DELR <- solve(V1 + lambda_DELR * V0 + lambda_DELR * V1) %*% (V1 + lambda_DELR * V0)
  
  DELR <- alpha1 + W_DELR %*% (beta0 - alpha1)
  
  loss_DELR <- t(DELR - beta) %*% H0 %*% (DELR - beta)
  
  c(loss_Prop, loss_beta0, loss_TLasso, loss_TGLM, loss_ISEDI, loss_DELR) 
}

p <- 25
n <- 100
m <- 200

set.seed(0)
sig_x0 <- sqrt(rgamma(p, 2, scale = 0.5))
Sig0 <- diag(sig_x0) %*% ar1_cor(p, 0.4) %*% diag(sig_x0)
sig_x <- sqrt(rgamma(p, 2, scale = 0.5))
Sigx <- diag(sig_x) %*% ar1_cor(p, 0.2) %*% diag(sig_x)
sim <- 5000

A_seq <- seq(0, 1, 0.1)

set.seed(0)
mu <-  rbinom(p + 1, 1, 0.8) * rnorm(p + 1, 0, sqrt(2 / (p + 1)))
eta <- rbinom(p + 1, 1, 0.5) * rnorm(p + 1, 0, sqrt(2 / (p + 1))) #runif(p, - sqrt(3), sqrt(3)) 
Res_seq <- array(NA, c(6, sim, length(A_seq)))

library(parallel)
library(MASS)


for (i in 1:length(A_seq)) {
  b <-  A_seq[i] * eta
  Sig1 <- (1 - A_seq[i]) * Sig0 + A_seq[i] * Sigx
  
  cl <- makeCluster(50)
  clusterEvalQ(cl, c(library(MASS), library(glmtrans), library(ISEDI), source("TransLasso-functions.R")))
  clusterExport(cl, c("mu", "Sig0", "Sig1", "p", "n", "m", "b", "q"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, , i] <- Res
}

save(Res_seq, file = paste(paste("Error", p, n, m, "Cmpr", "sameVariable", sep = "-"), ".Rdata", sep = ""))