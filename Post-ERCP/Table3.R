##confidence interval does not work well because sig0b is nearly singular, PCA is required
##the overcoverage does not come from the sample size, the randomness of PCA, the small/large PC
##it is due to the conservative nature of CATE
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

DATA <- medicaldata::indo_rct
DATA <- DATA[, -c(1, 2, 33)] #remove useless variables
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

y <- DATA$outcome
t <- DATA$treatment

n <- length(y)


g1 <- DATA$age <= 25
g2 <- DATA$age > 25 & DATA$age <= 35
g3 <- DATA$age > 35 & DATA$age <= 45
g4 <- DATA$age > 45 & DATA$age <= 55
g5 <- DATA$age > 55 & DATA$age <= 65
g6 <- DATA$age > 65 


beta0 <- 
  c(mean(y[t == 1]) - mean(y[t == 0]),
    mean(y[t == 1 & g1]) - mean(y[t == 0 & g1]),
    mean(y[t == 1 & g2]) - mean(y[t == 0 & g2]),
    mean(y[t == 1 & g3]) - mean(y[t == 0 & g3]),
    mean(y[t == 1 & g4]) - mean(y[t == 0 & g4]),
    mean(y[t == 1 & g5]) - mean(y[t == 0 & g5]),
    mean(y[t == 1 & g6]) - mean(y[t == 0 & g6]))

IF_beta <- c()

IF_tmp <- matrix(NA, n, 2)

for (tr in 0:1) {
  v1 <- as.numeric(t == tr)
  v2 <- as.numeric((t == tr) * y)
  IF_tmp[, tr + 1] <- (v2 - mean(v2)) / mean(v1) - mean(v2) / mean(v1)^2 * (v1 - mean(v1))
}

IF_beta <- cbind(IF_beta, IF_tmp[, 2] - IF_tmp[, 1])

for (tr in 0:1) {
  v1 <- as.numeric(t == tr & g1)
  v2 <- as.numeric((t == tr & g1) * y)
  IF_tmp[, tr + 1] <- (v2 - mean(v2)) / mean(v1) - mean(v2) / mean(v1)^2 * (v1 - mean(v1))
}

IF_beta <- cbind(IF_beta, IF_tmp[, 2] - IF_tmp[, 1])

for (tr in 0:1) {
  v1 <- as.numeric(t == tr & g2)
  v2 <- as.numeric((t == tr & g2) * y)
  IF_tmp[, tr + 1] <- (v2 - mean(v2)) / mean(v1) - mean(v2) / mean(v1)^2 * (v1 - mean(v1))
}

IF_beta <- cbind(IF_beta, IF_tmp[, 2] - IF_tmp[, 1])

for (tr in 0:1) {
  v1 <- as.numeric(t == tr & g3)
  v2 <- as.numeric((t == tr & g3) * y)
  IF_tmp[, tr + 1] <- (v2 - mean(v2)) / mean(v1) - mean(v2) / mean(v1)^2 * (v1 - mean(v1))
}

IF_beta <- cbind(IF_beta, IF_tmp[, 2] - IF_tmp[, 1])

for (tr in 0:1) {
  v1 <- as.numeric(t == tr & g4)
  v2 <- as.numeric((t == tr & g4) * y)
  IF_tmp[, tr + 1] <- (v2 - mean(v2)) / mean(v1) - mean(v2) / mean(v1)^2 * (v1 - mean(v1))
}

IF_beta <- cbind(IF_beta, IF_tmp[, 2] - IF_tmp[, 1])

for (tr in 0:1) {
  v1 <- as.numeric(t == tr & g5)
  v2 <- as.numeric((t == tr & g5) * y)
  IF_tmp[, tr + 1] <- (v2 - mean(v2)) / mean(v1) - mean(v2) / mean(v1)^2 * (v1 - mean(v1))
}

IF_beta <- cbind(IF_beta, IF_tmp[, 2] - IF_tmp[, 1])

for (tr in 0:1) {
  v1 <- as.numeric(t == tr & g6)
  v2 <- as.numeric((t == tr & g6) * y)
  IF_tmp[, tr + 1] <- (v2 - mean(v2)) / mean(v1) - mean(v2) / mean(v1)^2 * (v1 - mean(v1))
}

IF_beta <- cbind(IF_beta, IF_tmp[, 2] - IF_tmp[, 1])

library(MASS)
library(sphunif)
alpha1 <- c(5/57, 6/47,#young
            16/214, 19/212,#female
            20/246, 20/217,#pdstent
            4/64, 6/52,#difcan
            4/35, 3/40#bsphinc
)

var_lst <- c("young", "gender", "pdstent", "difcan", "bsphinc")

p <- length(beta0)
q <- length(alpha1)


alpha0 <- c()
for (var in var_lst) {
  sel <- DATA[, colnames(DATA) == var]
  for (tr in 0:1) {
    alpha0 <- c(alpha0, mean(DATA$outcome[sel == 1 & DATA$treatment == tr]))
  }
}

IF_alpha <- c()
for (var in var_lst) {
  sel <- DATA[, colnames(DATA) == var]
  for (tr in 0:1) {
    v1 <- as.numeric(sel == 1 & DATA$treatment == tr)
    v2 <- as.numeric((sel == 1 & DATA$treatment == tr) * DATA$outcome)
    IF_alpha <- cbind(IF_alpha, (v2 - mean(v2)) / mean(v1) - mean(v2) / mean(v1)^2 * (v1 - mean(v1)))
  }
}

Sig0a <- cov(IF_alpha) / n

Sig0b <- cov(IF_beta) / n
Sig0ba <-  cov(IF_beta, IF_alpha) / n

tmp_ev <- eigen(Sig0b)

H0 <- ginv(Sig0b)

Q_meta <- Sig0ba %*% ginv(Sig0a)

meta_dag <- beta0 - Q_meta %*% (alpha0 - alpha1)

eval_loss <- c()

Sig_d <- Sig0a
H_d <- t(Q_meta) %*% H0 %*% Q_meta
Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
A <- rep(1, p)
Pi_c <- diag(1, p) #- A %*% ginv(t(A) %*% ginv(Sig_c) %*% A) %*% t(A) %*% ginv(Sig_c)
P_c <- Pi_c %*% Sig_c %*% t(Pi_c)

PC_H0 <- eigen(H0)
HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
lambda_c <- max(sum(diag(P_c %*% H0)) - 2 * Re(eigen(HR0 %*% P_c %*% HR0)$values[1]), 0)
mu_meta <- (diag(1, p) - Pi_c) %*% meta_dag
r1 <- meta_dag - mu_meta
v_c <- t(r1) %*% H0 %*% r1
Prop <- mu_meta + max(1 - lambda_c / v_c, 0) * r1


A <- cbind(rep(1, q))
Pi_d <- diag(1, q) #- A %*% ginv(t(A) %*% ginv(Sig_d) %*% A) %*% t(A) %*% ginv(Sig_d)
P_d <- Pi_d %*% Sig_d %*% t(Pi_d)
P_b <- Q_meta %*% P_d %*% t(Q_meta)


PC_H_d <- eigen(H_d)
HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
lambda_d <- max(sum(diag(P_d %*% H_d)) - 2 * Re(eigen(HR_d %*% P_d %*% HR_d)$values[1]), 0)
mu_d <- (diag(1, q) - Pi_d) %*% (alpha0 - alpha1)
r2 <- Q_meta %*% (alpha0 - alpha1 - mu_d)
v_d <- t(r2) %*% H0 %*% r2
Prop <- Prop + (Q_meta %*% mu_d + max(1 - lambda_d / v_d, 0) * r2)


set.seed(0)
BT <- 10
B <- 500
LB <- matrix(NA, p, BT)
UB <- matrix(NA, p, BT)

Lik1 <- function(s) {
  S <- Sig_c + s * diag(1, p)
  t(meta_dag) %*% ginv(S) %*% meta_dag + log(det(S))
}

Lik2 <- function(s) {
  S <- Sig_d + s * diag(1, q)
  t(alpha0 - alpha1) %*% ginv(S) %*% (alpha0 - alpha1) + log(det(S))
}


s2_c <- nlminb(1, Lik1, lower = 0)$par
s2_d <- nlminb(1, Lik2, lower = 0)$par

meta_pos <- s2_c * ginv(Sig_c + s2_c * diag(1, p)) %*% meta_dag
Sig_c_pos <- s2_c * diag(1, p) - s2_c^{2} * ginv(Sig_c + s2_c * diag(1, p))

d_pos <- s2_d * ginv(Sig_d + s2_d * diag(1, q)) %*% (alpha0 - alpha1)
Sig_d_pos <- s2_d * diag(1, q) - s2_d^{2} * ginv(Sig_d + s2_d * diag(1, q))

for (t in 1:BT) {
  Boots <- c()
  meta_tmp <- mvrnorm(1, meta_pos, Sig_c_pos)
  d_tmp <- mvrnorm(1, d_pos, Sig_d_pos)
  for (j in 1:B) {
    meta_b <- mvrnorm(1, meta_tmp, Sig_c)
    d_b <- mvrnorm(1, d_tmp, Sig_d)
    
    mu_meta_b <- (diag(1, p) - Pi_c) %*% meta_b
    r1_b <- meta_b - mu_meta_b
    v_c_b <- t(r1_b) %*% H0 %*% r1_b
    Prop_b <- mu_meta_b + max(0, c(1 - lambda_c / v_c_b)) * r1_b
    
    mu_d_b <- (diag(1, q) - Pi_d) %*% d_b
    H_d <- t(Q_meta) %*% H0 %*% Q_meta
    r2_b <- Q_meta %*% (d_b - mu_d)
    v_d_b <- t(r2_b) %*% H0 %*% r2_b
    Prop_b <- Prop_b + (Q_meta %*% mu_d_b + max(0, c(1 - lambda_d / v_d_b)) * r2_b)
    
    Boots <- cbind(Boots, Prop_b - meta_tmp - Q_meta %*% d_tmp)
  }
  LB[, t] <- apply(Boots, 1, function(x) quantile(x, 0.025))
  UB[, t] <- apply(Boots, 1, function(x) quantile(x, 0.975))
}

sd <- sqrt(diag(Sig0b))

EST_tar <- beta0

EST_Prop <- Prop

CI_tar <- cbind(beta0 - sd * qnorm(0.975), beta0 + sd * qnorm(0.975))

CI_Prop <- cbind(Prop - apply(UB, 1, max), Prop - apply(LB, 1, min))

library(plyr)

round_any(c(EST_tar[1], CI_tar[1, 1], CI_tar[1, 2], EST_Prop[1], CI_Prop[1, 1], CI_Prop[1, 2]), accuracy =  0.01, f = ceiling)

round_any(cbind(EST_tar[-1], CI_tar[-1, 1], CI_tar[-1, 2], CI_tar[-1, 2] - CI_tar[-1, 1], 
                EST_Prop[-1], CI_Prop[-1, 1], CI_Prop[-1, 2], CI_Prop[-1, 2] - CI_Prop[-1, 1]), accuracy =  0.01, f = ceiling)

