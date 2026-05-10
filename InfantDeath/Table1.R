setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(MASS)

load("Infant_Death_2015.Rdata")

dat.White <- df[df$MWhite == 1, - c(6:21)]

dat.Black <- df[df$MBlack == 1, - c(6:21)]

dat.Target <- df[df$MFilipino + df$MViet == 1, - c(6:21)]


n_White <- nrow(dat.White)
n_Black <- nrow(dat.Black)
n_Target <- nrow(dat.Target)

n_case_White <- sum(dat.White$Y)
n_case_Black <- sum(dat.Black$Y)
n_case_Target <- sum(dat.Target$Y)

w <- rep(0, n_Target)

w[dat.Target$Y == 1] <- 1
w[dat.Target$Y == 0] <- 5 * n_case_Target / (n_Target - n_case_Target)

fit <- glm(Y ~ ., family = "binomial", data = dat.Target, weight = w)

beta0 <- fit$coefficients 

beta_T <- beta0
p <- length(beta_T) - 1

simone <- function(sim) {
  set.seed(sim)
  
  tmp.Target <- rbind(dat.Target[dat.Target$Y == 1, ][sample(1:n_case_Target, n_case_Target, replace = T), ],
                      dat.Target[dat.Target$Y == 0, ][sample(1:(n_Target - n_case_Target), 5 * n_case_Target, replace = T), ])
  tmp.White <- rbind(dat.White[dat.White$Y == 1, ][sample(1:n_case_White, n_case_White, replace = T), ],
                      dat.White[dat.White$Y == 0, ][sample(1:(n_White - n_case_White), 5 * n_case_White, replace = T), ])
  tmp.Black <- rbind(dat.Black[dat.Black$Y == 1, ][sample(1:n_case_Black, n_case_Black, replace = T), ],
                      dat.Black[dat.Black$Y == 0, ][sample(1:(n_Black - n_case_Black), 5 * n_case_Black, replace = T), ])
  n0 <- nrow(tmp.Target)
  
  n1_White <- nrow(tmp.White)
  
  n1_Black <- nrow(tmp.Black)
  
  fit <- glm(Y ~ ., family = "binomial", data = tmp.White)
  
  alpha1_White <- fit$coefficients 
  
  fit <- glm(Y ~ ., family = "binomial", data = tmp.Black)
  
  alpha1_Black <- fit$coefficients 
  
  fit <- logistf(Y ~ ., data = tmp.Target)
  
  beta0 <- fit$coefficients 
  
  X1 <- as.matrix(cbind(1, tmp.White[, -1]))
  
  p_1 <- c(plogis(X1 %*% alpha1_White))
  G <- t(X1) %*% (p_1 * (1 - p_1) * X1) / n1_White
  
  IF_alpha1 <- ((tmp.White$Y - p_1) * X1) %*% ginv(G)
  
  Sig1a_White <- cov(IF_alpha1) / n1_White
  
  alpha1_White <- alpha1_White
  
  X1 <- as.matrix(cbind(1, tmp.Black[, -1]))
  
  p_1 <- c(plogis(X1 %*% alpha1_Black))
  G <- t(X1) %*% (p_1 * (1 - p_1) * X1) / n1_Black
  
  IF_alpha1 <- ((tmp.Black$Y - p_1) * X1) %*% ginv(G)
  
  Sig1a_Black <- cov(IF_alpha1) / n1_Black
  
  alpha1_Black <- alpha1_Black
  
  alpha1 <- ginv(ginv(Sig1a_White, tol = 1e-4) + ginv(Sig1a_Black, tol = 1e-4)) %*% 
            (ginv(Sig1a_White, tol = 1e-4) %*% alpha1_White + ginv(Sig1a_Black, tol = 1e-4) %*% alpha1_Black)
  
  Sig1a <- ginv(ginv(Sig1a_White, tol = 1e-4) + ginv(Sig1a_Black, tol = 1e-4))
  
  X0 <- as.matrix(cbind(1, tmp.Target[, -1]))
  
  alpha0 <- beta0
  
  p_0 <- c(plogis(X0 %*% beta0))
  G <- t(X0) %*% (p_0 * (1 - p_0) * X0) / n0
  
  Sig_S <- G * n0
  
  IF_beta <- ((tmp.Target$Y - p_0) * X0) %*% ginv(G)
  
  Sig0b <- cov(IF_beta) / n0
  Sig0ba <-  Sig0b
  Sig0a <-  Sig0b
  
  H0 <- t(X0) %*% ((p_0 * (1 - p_0))^2 * X0) / n0

  #############DELogis second order optimal weight
  
  fit <- glm(Y ~ ., family = "binomial", data = rbind(tmp.White, tmp.Black))
  
  beta_pool <- fit$coefficients 
  
  
  X1 <- as.matrix(cbind(1, rbind(tmp.White[, -1], tmp.Black[, -1])))
  
  p_1 <- c(plogis(X1 %*% beta_pool))
  Sig_B_pool <- t(X1) %*% (p_1 * (1 - p_1) * X1)
  
  
  gamma_pool <- beta_pool - beta0
  
  n1_pool <- n1_White + n1_Black
  
  W_pool <- ginv(gamma_pool %*% t(gamma_pool) + ginv(Sig_B_pool, tol = 1e-4 * n1_pool) + ginv(Sig_S, tol = 1e-4 * n0)) %*% (gamma_pool %*% t(gamma_pool) + ginv(Sig_B_pool, tol = 1e-4 * n1_pool))
  
  DELogis <- beta_pool - W_pool %*% gamma_pool
  
  ##########calibration
  
  Q_meta <- Sig0ba %*% ginv(Sig0a + Sig1a)
  
  meta <- beta0 - Q_meta %*% (alpha0 - alpha1)
  
  ###########dShrink
  
  Sig_d <- Sig0a + Sig1a
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  
  lambda_c <- max(sum(diag(Sig_c %*% H0)) - 2 * Re(eigen(HR0 %*% Sig_c %*% HR0)$values[1]), 0)
  r1 <- meta
  v_c <- t(r1) %*% H0 %*% r1
  Prop <- max(0, c(1 - lambda_c / v_c)) * r1
  
  
  PC_H_d <- eigen(H_d, symmetric = T)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda_d <- max(sum(diag(Sig_d %*% H_d)) - 2 * eigen(HR_d %*%  Sig_d %*% HR_d, symmetric = T)$values[1], 0)
  r2 <- Q_meta %*% (alpha0 - alpha1)
  v_d <- t(r2) %*% H0 %*% r2
  Prop <- Prop + max(1 - lambda_d / v_d, 0) * r2
  
  
  H0 <- ginv(Sig0b)
  
  H0 <- H0 + 0.2 * max(diag(H0)) * diag(1, p + 1)
  
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  
  lambda_c <- max(sum(diag(Sig_c %*% H0)) - 2 * Re(eigen(HR0 %*% Sig_c %*% HR0)$values[1]), 0)
  r1 <- meta
  v_c <- t(r1) %*% H0 %*% r1
  Prop_Inf <- max(0, c(1 - lambda_c / v_c)) * r1
  
  
  PC_H_d <- eigen(H_d, symmetric = T)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda_d <- max(sum(diag(Sig_d %*% H_d)) - 2 * eigen(HR_d %*%  Sig_d %*% HR_d, symmetric = T)$values[1], 0)
  r2 <- Q_meta %*% (alpha0 - alpha1)
  v_d <- t(r2) %*% H0 %*% r2
  Prop_Inf <- Prop_Inf + max(1 - lambda_d / v_d, 0) * r2
  
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
      Prop_b <- max(0, c(1 - lambda_c / v_c)) * r1_b
      
      r2_b <- Q_meta %*% d_b
      v_d_b <- t(r2_b) %*% H0 %*% r2_b
      Prop_b <- Prop_b + max(0, c(1 - lambda_d / v_d_b)) * r2_b
      
      Boots <- cbind(Boots, Prop_b - meta_tmp - Q_meta %*% d_tmp)
    }
    LB[, t] <- apply(Boots, 1, function(x) quantile(x, 0.025))
    UB[, t] <- apply(Boots, 1, function(x) quantile(x, 0.975))
  }
  
  cri1 <- sqrt(diag(Sig0b)) * qnorm(0.975)
  cri2 <- sqrt(diag(Sig0b - Sig0ba %*% ginv(Sig0a + Sig1a) %*% t(Sig0ba))) * qnorm(0.975)
  cri3 <- cbind(apply(LB, 1, min), apply(UB, 1, max))
  
  TGLM_est <- NA#c(glmtrans(target = list(x = tmp.Target[, -1], y = tmp.Target$Y), 
                         # source = list(list(x = tmp.White[, -1], y = tmp.White$Y),
                         #               list(x = tmp.Black[, -1], y = tmp.Black$Y)),
                         # family = "binomial")$beta)
  
  TGLM <- NA#glmtrans_inf(target = list(x = tmp.Target[, -1], y = tmp.Target$Y), 
                       # source = list(list(x = tmp.White[, -1], y = tmp.White$Y),
                       #               list(x = tmp.Black[, -1], y = tmp.Black$Y)),
                       # beta.hat = TGLM_est, family = "binomial")
  
  lambda_seq <- seq(0,7,length.out=10)
  X0 <- as.matrix(cbind(1, tmp.Target[, -1]))
  X1 <- as.matrix(cbind(1, rbind(tmp.White, tmp.Black)[, -1]))
  Y0 <- tmp.Target$Y
  Y1 <- c(rbind(tmp.White, tmp.Black)[, 1])
  ISEDI <- LDW_meta(X1, X0, Y1, Y0, lambda_seq, family="binomial")$estimates_opt
  
  ISEDI_est <- ISEDI$estimate
  
  
  length_beta0 <- mean(2 * cri1)
  length_meta <- mean(2 * cri2)
  length_TGLM <- mean(TGLM$CI[, 3] - TGLM$CI[, 2])
  length_ISEDI <- mean(ISEDI$upper_CI[-1] - ISEDI$lower_CI[-1])
  length_Prop <- mean(cri3[, 2] - cri3[, 1])
  
  cover_beta0 <- mean(abs(beta0 - beta_T) <= cri1)
  cover_meta <- mean(abs(meta - beta_T) <= cri2)
  cover_TGLM <- mean(TGLM$CI[, 3] >= beta_T & TGLM$CI[, 2] <= beta_T)
  cover_ISEDI <- mean((beta_T >= ISEDI$lower_CI) & (beta_T <= ISEDI$upper_CI))
  cover_Prop <- mean(((Prop_Inf - beta_T) >= cri3[, 1]) & ((Prop_Inf - beta_T) <= cri3[, 2]))
  
  y_ref <- dat.Target$Y
  x_ref <- as.matrix(cbind(1, dat.Target[, -1]))
  
  w <- w / sum(w)
  
  c(sum((plogis(x_ref %*% beta0) - plogis(x_ref %*% beta_T))^2 * w),
    sum((plogis(x_ref %*% meta) - plogis(x_ref %*% beta_T))^2 * w),
    sum((plogis(x_ref %*% DELogis) - plogis(x_ref %*% beta_T))^2 * w),
    sum((plogis(x_ref %*% TGLM_est) - plogis(x_ref %*% beta_T))^2 * w),
    sum((plogis(x_ref %*% ISEDI_est) - plogis(x_ref %*% beta_T))^2 * w),
    sum((plogis(x_ref %*% Prop) - plogis(x_ref %*% beta_T))^2 * w),
    length_beta0,
    length_meta,
    length_TGLM,
    length_ISEDI,
    length_Prop,
    cover_beta0,
    cover_meta,
    cover_TGLM,
    cover_ISEDI,
    cover_Prop)
}


library(parallel)

cl <- makeCluster(20)
clusterEvalQ(cl, c(library(MASS), library(logistf), library(glmtrans), library(ISEDI)))
clusterExport(cl, c("dat.Target", "dat.White", "dat.Black", "n_Target", "n_White", "n_Black", "n_case_Target", "n_case_White", "n_case_Black", "w", "p", "q", "beta_T"))
Res1 <- parSapply(cl, 1:200, simone, simplify = "array")
stopCluster(cl)

save(Res, file = "Res_Infant.Rdata")


library(ggplot2)
Res <- data.frame(cbind(rep(c("Target", "Calibration", "dShrink"), c(200, 200, 200)), c(t(Res[c(1, 2, 6), ]))))
colnames(Res) <- c("Method", "MSE")
Res$MSE <- as.numeric(Res$MSE)
ggplot(Res, aes(x=factor(Method, levels = c("Target", "Calibration", "dShrink")), y=MSE, fill=Method)) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot() + 
  ylab("Squared Error") +
  xlab("") +
  theme_light() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.title = element_blank(),                   
        legend.text = element_text(size = 20, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=20),
        legend.key.size=unit(1.5,'cm'))

ggsave(paste(paste("Boxplot", "Infant", "WhiteS", sep = "_"), ".pdf", sep = ""),
       path = ".", device = "pdf", width = 12, height = 9, units = "in")
