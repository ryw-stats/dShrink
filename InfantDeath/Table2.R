setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(MASS)

load("Infant_Death_2015.Rdata")

dat.White <- df[df$MWhite == 1, - c(6:21)]

dat.Black <- df[df$MBlack == 1, - c(6:21)]

dat.Target <- df[df$MFilipino + df$MViet == 1, - c(6:21)]

n0 <- nrow(dat.Target)

fit <- glm(Y ~ ., family = "binomial", data = dat.Target)

beta0 <- fit$coefficients 

beta_T <- beta0[-1]
p <- length(beta_T)

X0 <- as.matrix(cbind(1, dat.Target[, -1]))

alpha0 <- beta0[-1]

p_0 <- c(plogis(X0 %*% beta0))
G <- t(X0) %*% (p_0 * (1 - p_0) * X0) / n0

IF_beta <- ((dat.Target$Y - p_0) * X0) %*% ginv(G)

Sig0b <- cov(IF_beta)[-1, -1] / n0
Sig0ba <-  Sig0b
Sig0a <-  Sig0b

beta0 <- beta0[-1]

H0 <- ginv(Sig0b)

H0 <- H0 + 0.2 * max(diag(H0)) * diag(1, p)

n1_White <- nrow(dat.White)

n1_Black <- nrow(dat.Black)

fit <- glm(Y ~ ., family = "binomial", data = dat.White)

alpha1_White <- fit$coefficients 

fit <- glm(Y ~ ., family = "binomial", data = dat.Black)

alpha1_Black <- fit$coefficients 

X1 <- as.matrix(cbind(1, dat.White[, -1]))

p_1 <- c(plogis(X1 %*% alpha1_White))
G <- t(X1) %*% (p_1 * (1 - p_1) * X1) / n1_White

IF_alpha1 <- ((dat.White$Y - p_1) * X1) %*% ginv(G)

Sig1a_White <- cov(IF_alpha1)[-1, -1] / n1_White

alpha1_White <- alpha1_White[-1]

X1 <- as.matrix(cbind(1, dat.Black[, -1]))

p_1 <- c(plogis(X1 %*% alpha1_Black))
G <- t(X1) %*% (p_1 * (1 - p_1) * X1) / n1_Black

IF_alpha1 <- ((dat.Black$Y - p_1) * X1) %*% ginv(G)

Sig1a_Black <- cov(IF_alpha1)[- 1, - 1] / n1_Black

alpha1_Black <- alpha1_Black[- 1]

alpha1 <- ginv(ginv(Sig1a_White, tol = 1e-4) + ginv(Sig1a_Black, tol = 1e-4)) %*% 
  (ginv(Sig1a_White, tol = 1e-4) %*% alpha1_White + ginv(Sig1a_Black, tol = 1e-4) %*% alpha1_Black)

Sig1a <- ginv(ginv(Sig1a_White, tol = 1e-4) + ginv(Sig1a_Black, tol = 1e-4))

Q_meta <- Sig0ba %*% ginv(Sig0a + Sig1a)

meta <- beta0 - Q_meta %*% (alpha0 - alpha1)

Sig_d <- Sig0a + Sig1a
H_d <- t(Q_meta) %*% H0 %*% Q_meta
Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)


PC_H0 <- eigen(H0)
HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
lambda_c <- max(sum(diag(Sig_c %*% H0)) - 2 * eigen(HR0 %*% Sig_c %*% HR0, symmetric = T)$values[1], 0)
r1 <- meta
v_c <- t(r1) %*% H0 %*% r1
Prop <- max(1 - lambda_c / v_c, 0) * r1


PC_H_d <- eigen(H_d, symmetric = T)
HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
lambda_d <- max(sum(diag(Sig_d %*% H_d)) - 2 * eigen(HR_d %*% Sig_d %*% HR_d, symmetric = T)$values[1], 0)
r2 <- Q_meta %*% (alpha0 - alpha1)
v_d <- t(r2) %*% H0 %*% r2
Prop <- Prop +  max(1 - lambda_d / v_d, 0) * r2

set.seed(0)
BT <- 10
B <- 500
LB <- matrix(NA, p, BT)
UB <- matrix(NA, p, BT)

Lik1 <- function(s) {
  S <- Sig_c + s * diag(1, p)
  t(meta) %*% ginv(S) %*% meta + log(det(S))
}

Lik2 <- function(s) {
  S <- Sig_d + s * diag(1, p)
  t(alpha0 - alpha1) %*% ginv(S) %*% (alpha0 - alpha1) + log(det(S))
}


s2_c <- nlminb(1, Lik1, lower = 0)$par
s2_d <- nlminb(1, Lik2, lower = 0)$par

meta_pos <- s2_c * ginv(Sig_c + s2_c * diag(1, p)) %*% meta
Sig_c_pos <- s2_c * diag(1, p) - s2_c^{2} * ginv(Sig_c + s2_c * diag(1, p))

d_pos <- s2_d * ginv(Sig_d + s2_d * diag(1, p)) %*% (alpha0 - alpha1)
Sig_d_pos <- s2_d * diag(1, p) - s2_d^{2} * ginv(Sig_d + s2_d * diag(1, p))

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

library(glmtrans)

# TGLM_est <- c(glmtrans(target = list(x = dat.Target[, -1], y = dat.Target$Y), 
#                        source = list(list(x = dat.White[, -1], y = dat.White$Y),
#                                      list(x = dat.Black[, -1], y = dat.Black$Y)),
#                        family = "binomial")$beta)
# 
# TGLM <- glmtrans_inf(target = list(x = dat.Target[, -1], y = dat.Target$Y), 
#                      source = list(list(x = dat.White[, -1], y = dat.White$Y),
#                                    list(x = dat.Black[, -1], y = dat.Black$Y)),
#                      beta.hat = TGLM_est, family = "binomial", cores = 2)


Res_tar <- cbind(beta0, beta0 - cri1, beta0 + cri1)

Res_ds <- cbind(Prop, Prop - cri2[, 2], Prop - cri2[, 1])

rownames(Res_ds) <- rownames(Res_tar)

Res <- list(Res_tar, Res_ds)

save(Res, file = "RealData_Infant.Rdata")


library(plyr)

Res_tar <- Res[[1]]

Res_ds <- Res[[2]]

cbind(round_any(cbind(Res_tar[c(15, 13, 12, 16, 17), ], Res_tar[c(15, 13, 12, 16, 17), 3] - Res_tar[c(15, 13, 12, 16, 17), 2]), accuracy =  0.01, f = ceiling),
round_any(cbind(Res_ds[c(15, 13, 12, 16, 17), ], Res_ds[c(15, 13, 12, 16, 17), 3] - Res_ds[c(15, 13, 12, 16, 17), 2]), accuracy =  0.01, f = ceiling))

round_any(mean(Res_ds[, 3] - Res_ds[, 2]), accuracy =  0.001, f = ceiling)
round_any(mean(Res_tar[, 3] - Res_tar[, 2]), accuracy =  0.001, f = ceiling)

