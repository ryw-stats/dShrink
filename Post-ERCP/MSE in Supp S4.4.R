setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

DATA <- medicaldata::indo_rct
DATA <- DATA[, -c(1, 2, 33)]
DATA[, c(1, 4)] <- DATA[, c(4, 1)]
colnames(DATA)[c(1, 4)] <- colnames(DATA)[c(4, 1)]
DATA[, c(2, 30)] <- DATA[, c(30, 2)]
colnames(DATA)[c(2, 30)] <- c("treatment", "risk")
for (i in 1:ncol(DATA)) {
  DATA[, i] <- as.numeric(unlist(DATA[, i]))
}
DATA <- DATA - 1
DATA$gender <- 1 - DATA$gender 
DATA$young <- as.numeric(DATA$age < 50)

X <- as.matrix(DATA[, -c(1, 2, 31)])

RX <- cor(X)

r <- sqrt(1 - 1 / diag(solve(RX)) / diag(RX))

X_aug <- cbind(1, X[, r <= 0.95])

library(MASS)
library(sphunif)
alpha1 <- c(5/57, 6/47,#young
            16/214, 19/212,#female
            20/246, 20/217,#pdstent
            4/64, 6/52,#difcan
            4/35, 3/40#bsphinc
)

var_lst <- c("young", "gender", "pdstent", "difcan", "bsphinc")

y_tilde <- 2 * DATA$treatment * DATA$outcome - 2 * (1 - DATA$treatment) * DATA$outcome

p <- ncol(X_aug)
loss <- function(b) {
  mean((y_tilde - 2 * plogis(X_aug %*% b) + 1)^2)
}

beta <- nlminb(rep(0, p), loss, control = list(iter = 1000))$par

N <- nrow(X_aug)

# p_lin <- X_aug %*% beta
# h_or <- c(8 * plogis(p_lin)^2 * (1 - plogis(p_lin))^2)
# 
# H_or <- t(X_aug) %*% (h_or * X_aug) / N
# Var_beta <- solve(H_or) %*% cov(- 4 * c(plogis(p_lin) * (1 - plogis(p_lin)) * (y_tilde - 2 * plogis(p_lin) + 1)) * X_aug) %*% solve(H_or) / N
# 1 - pchisq(t(beta[-1]) %*% solve(Var_beta[-1, -1]) %*% beta[-1], p - 1)

Res <- c()

simone <- function(sim) {
  set.seed(sim)
  ind <- sample(1:N, N, replace = T)
  sim_DATA <- DATA[ind, ]
  X <- X_aug[ind, ]
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  q <- length(alpha1)
  
  
  alpha0 <- c()
  for (var in var_lst) {
    sel <- sim_DATA[, colnames(sim_DATA) == var]
    for (tr in 0:1) {
      alpha0 <- c(alpha0, mean(sim_DATA$outcome[sel == 1 & sim_DATA$treatment == tr]))
    }
  }
  
  IF_alpha <- c()
  for (var in var_lst) {
    sel <- DATA[ind, colnames(DATA) == var]
    for (tr in 0:1) {
      v1 <- as.numeric(sel == 1 & DATA$treatment[ind] == tr)
      v2 <- as.numeric((sel == 1 & DATA$treatment[ind] == tr) * DATA$outcome[ind])
      IF_alpha <- cbind(IF_alpha, (v2 - mean(v2)) / mean(v1) - mean(v2) / mean(v1)^2 * (v1 - mean(v1)))
    }
  }
  
  Sig0a <- cov(IF_alpha) / n
  
  y <- y_tilde[ind]
  
  #lambda_min <- cv.glmnet(X[, -1], y, family = "gaussian", lambda = seq(0.5, 1, 0.1) * sqrt(log(p) / n))$lambda.min
  
  loss <- function(b) {
    mean((y - 2 * plogis(X %*% b) + 1)^2) + 0.1 * sum(b^2) / n
  }
  
  beta0 <- nlminb(rep(0, p), loss, lower = -5, upper = 5, control = list(iter = 1000))$par
  
  p_lin <- X %*% beta0
  h0b <- c(8 * plogis(p_lin)^2 * (1 - plogis(p_lin))^2)
  H0 <- t(X) %*% (h0b * X) / n + diag(0.2, p) / n
  
  
  IF_beta <- - (- 4 * c(plogis(p_lin) * (1 - plogis(p_lin)) * (y_tilde - 2 * plogis(p_lin) + 1)) * X) %*% ginv(H0)
  
  Sig0b <- cov(IF_beta) / n
  Sig0ba <-  cov(IF_beta, IF_alpha) / n
  
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  
  Q_meta <- Sig0ba %*% ginv(Sig0a)
  
  meta_dag <- beta0 - Q_meta %*% (alpha0 - alpha1)
  
  Sig_d <- Sig0a
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
  A <- cbind(Sig_c %*% rep(1, p))
  Pi <- diag(1, p) #- A %*% ginv(t(A) %*% ginv(Sig_c) %*% A) %*% t(A) %*% ginv(Sig_c)
  P_c <- Pi %*% Sig_c %*% t(Pi)
  
  
  lambda <- max(sum(diag(P_c %*% H0)) - 2 * Re(eigen(HR0 %*% P_c %*% HR0)$values[1]), 0)
  mu_meta <- (diag(1, p) - Pi) %*% meta_dag
  r1 <- meta_dag - mu_meta
  v_c <- t(r1) %*% H0 %*% r1
  Prop <- mu_meta + max(1 - lambda / v_c, 0) * r1
  
  
  A <- rep(1, q)
  Pi <- diag(1, q) #- A %*% ginv(t(A) %*% ginv(Sig_d) %*% A) %*% t(A) %*% ginv(Sig_d)
  P_d <- Pi %*% Sig_d %*% t(Pi)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda <- max(sum(diag(P_d %*% H_d)) - 2 * Re(eigen(HR_d %*% P_d %*% HR_d)$values[1]), 0)
  mu_d <- (diag(1, q) - Pi) %*% (alpha0 - alpha1)
  r2 <- Q_meta %*% (alpha0 - alpha1 - mu_d)
  v_d <- t(r2) %*% H0 %*% r2
  Prop <- Prop + (Q_meta %*% mu_d + max(1 - lambda / v_d, 0) * r2)
  
  
  c(4 * mean((plogis(X_aug %*% beta) - plogis(X_aug %*% beta0))^2), 
    4 * mean((plogis(X_aug %*% beta) - plogis(X_aug %*% Prop))^2)
    )
}

library(parallel)

cl <- makeCluster(20)
clusterEvalQ(cl, c(library("MASS"), library(sphunif)))
clusterExport(cl, c("y_tilde", "X_aug", "DATA", "N", "alpha1", "var_lst", "beta"))
Res <- parSapply(cl, 1:5000, simone, simplify = "array")
stopCluster(cl)

ER <- rowMeans(Res) * N
ER
save(ER, file = "ITE_PRED.Rdata")