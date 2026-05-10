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
  
  H0 <- solve(Sig0b)
  
  Sig_d <- Sig0a + Sig1a
  Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
  
  Pi_c <- diag(1, p + 1) #- A %*% solve(t(A) %*% solve(Sig_c) %*% A) %*% t(A) %*% solve(Sig_c)
  P_c <- Pi_c %*% Sig_c %*% t(Pi_c)
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda_c <- max(sum(diag(P_c %*% H0)) - 2 * Re(eigen(HR0 %*% P_c %*% HR0)$values[1]), 0)
  r1 <- meta
  v_c <- t(r1) %*% H0 %*% r1
  Prop <- max(0, c(1 - lambda_c / v_c)) * r1
  
  Pi_d <- diag(1, p + 1) #- matrix(1, q + 1, q + 1) %*% solve(Sig_d) / sum(solve(Sig_d))
  P_d <- Pi_d %*% Sig_d %*% t(Pi_d)
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r2 <- Q_meta %*% (alpha0 - alpha1)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda_d <- max(sum(diag(P_d %*% H_d)) - 2 * Re(eigen(HR_d %*% P_d %*% HR_d)$values[1]), 0)
  v_d <- t(r2) %*% H0 %*% r2
  Prop <- Prop + max(0, c(1 - lambda_d / v_d)) * r2
  
  set.seed(0)
  BT <- 10
  B <- 500
  LB <- matrix(NA, p + 1, BT)
  UB <- matrix(NA, p + 1, BT)
  
  Lik1 <- function(s) {
    S <- Sig_c + s * diag(1, p + 1)
    t(meta) %*% ginv(S) %*% meta + log(det(S))
  }
  
  Lik2 <- function(s) {
    S <- Sig_d + s * diag(1, p + 1)
    t(alpha0 - alpha1) %*% ginv(S) %*% (alpha0 - alpha1) + log(det(S))
  }
  
  
  s2_c <- nlminb(1, Lik1, lower = 0)$par
  s2_d <- nlminb(1, Lik2, lower = 0)$par
  
  meta_pos <- s2_c * ginv(Sig_c + s2_c * diag(1, p + 1)) %*% meta
  Sig_c_pos <- s2_c * diag(1, p + 1) - s2_c^{2} * ginv(Sig_c + s2_c * diag(1, p + 1))
  
  d_pos <- s2_d * ginv(Sig_d + s2_d * diag(1, p + 1)) %*% (alpha0 - alpha1)
  Sig_d_pos <- s2_d * diag(1, p + 1) - s2_d^{2} * ginv(Sig_d + s2_d * diag(1, p + 1))
  
  for (t in 1:BT) {
    Boots <- c()
    meta_tmp <- mvrnorm(1, meta_pos, Sig_c_pos)
    d_tmp <- mvrnorm(1, d_pos, Sig_d_pos)
    for (j in 1:B) {
      meta_b <- mvrnorm(1, meta_tmp, Sig_c)
      d_b <- mvrnorm(1, d_tmp, Sig_d)
      
      r1_b <- meta_b
      v_c_b <- t(r1_b) %*% H0 %*% r1_b
      Prop_b <- max(0, c(1 - lambda_c / v_c_b)) * r1_b
      
      H_d <- t(Q_meta) %*% H0 %*% Q_meta
      r2_b <- Q_meta %*% d_b
      v_d_b <- t(r2_b) %*% H0 %*% r2_b
      Prop_b <- Prop_b + max(0, c(1 - lambda_d / v_d_b)) * r2_b
      
      Boots <- cbind(Boots, Prop_b - meta_tmp - Q_meta %*% d_tmp)
    }
    LB[, t] <- apply(Boots, 1, function(x) quantile(x, 0.025))
    UB[, t] <- apply(Boots, 1, function(x) quantile(x, 0.975))
  }
  
  cri1 <- sqrt(diag(Sig0b)) * qnorm(0.975)
  cri2 <- cbind(apply(LB, 1, min), apply(UB, 1, max))
  
  
  TGLM_est <- c(glmtrans(target = list(x = x0, y = y0), source = list(list(x = x1, y = y1)), family = "gaussian")$beta)
  
  TGLM <- glmtrans_inf(target = list(x = x0, y = y0), source = list(list(x = x1, y = y1)), beta.hat = TGLM_est, family = "gaussian")
  
  lambda_seq <- seq(0,7,length.out=10)
  ISEDI <- LDW_meta(cbind(1, x1), cbind(1, x0), y1, y0, lambda_seq, family="gaussian")$estimates_opt[, 3:4]
  
  length_Prop <- mean(cri2[, 2] - cri2[, 1])
  length_beta0 <- mean(2 * cri1)
  length_TGLM <- mean(TGLM$CI[, 3] - TGLM$CI[, 2])
  length_ISEDI <- mean(ISEDI[, 2] - ISEDI[, 1])
  
  cover_Prop <- mean(((Prop - beta) >= cri2[, 1]) & ((Prop - beta) <= cri2[, 2]))
  cover_beta0 <- mean(abs(beta0 - beta) <= cri1)
  cover_TGLM <- mean(TGLM$CI[, 3] >= beta & TGLM$CI[, 2] <= beta)
  cover_ISEDI <- mean(ISEDI[, 2] >= beta & ISEDI[, 1] <= beta)
  
  c(length_Prop, length_beta0, length_TGLM, length_ISEDI, cover_Prop, cover_beta0, cover_TGLM, cover_ISEDI) 
}

p <- 25
n <- 100
m <- 200

set.seed(0)
sig_x0 <- sqrt(rgamma(p, 2, scale = 0.5))
Sig0 <- diag(sig_x0) %*% ar1_cor(p, 0.4) %*% diag(sig_x0)
sig_x <- sqrt(rgamma(p, 2, scale = 0.5))
Sigx <- diag(sig_x) %*% ar1_cor(p, 0.2) %*% diag(sig_x)
sim <- 500

A_seq <- seq(0, 1, 0.1)

set.seed(0)
mu <-  rbinom(p + 1, 1, 0.8) * rnorm(p + 1, 0, sqrt(2 / (p + 1)))
eta <- rbinom(p + 1, 1, 0.5) * rnorm(p + 1, 0, sqrt(2 / (p + 1))) #runif(p, - sqrt(3), sqrt(3)) 
Res_seq <- array(NA, c(8, sim, length(A_seq)))

library(parallel)
library(MASS)


for (i in 1:length(A_seq)) {
  b <-  A_seq[i] * eta
  Sig1 <- (1 - A_seq[i]) * Sig0 + A_seq[i] * Sigx
  
  cl <- makeCluster(50)
  clusterEvalQ(cl, c(library(MASS), library(glmtrans), library(ISEDI)))
  clusterExport(cl, c("mu", "Sig0", "Sig1", "p", "n", "m", "b", "q"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, , i] <- Res
}

save(Res_seq, file = paste(paste("Error", p, n, m, "Cmpr", "sameVariable", "CI", sep = "-"), ".Rdata", sep = ""))