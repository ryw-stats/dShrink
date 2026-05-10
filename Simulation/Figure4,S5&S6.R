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
  Sig1a <- cov(c(y - z %*% alpha1) * (z %*% M_z1)) / m
  
  
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
  Sig0b <- cov(c(y - x_aug %*% beta0) * (x_aug %*% M_x0)) / n
  Sig0a <- cov(c(y - z %*% alpha0) * (z %*% M_z0)) / n
  Sig0ba <- cov(c(y - x_aug %*% beta0) * (x_aug %*% M_x0), c(y - z %*% alpha0) * (z %*% M_z0)) / n
  
  Q_meta <- Sig0ba %*% solve(Sig0a + Sig1a)
  meta <- beta0 - Q_meta %*% (alpha0 - alpha1)
  
  sigma_comp <- sqrt(diag(Sig0a) + diag(Sig1a))
  
  test <- (abs(alpha0 - alpha1) < qnorm(0.975) * sigma_comp)
  comp_test <- beta0
  if (sum(test) > 0) {
    comp_test <- comp_test - Sig0ba[, test] %*% solve(Sig0a[test, test] + Sig1a[test, test]) %*% (alpha0[test] - alpha1[test])
  }
  
  test <- c(t(alpha0 - alpha1) %*% solve(Sig0a + Sig1a) %*% (alpha0 - alpha1)  < qchisq(0.95, q))
  
  global_test <- beta0 + (meta - beta0) * test
  
  H0 <- Sig_x0_aug
  Prop <- c()
  
  Sig_d <- Sig0a + Sig1a
  Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
  P_c <- Sig_c
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda <- max(sum(diag(P_c %*% H0)) - 2 * Re(eigen(HR0 %*% P_c %*% HR0)$values[1]), 0)
  r1 <- meta
  v_c <- t(r1) %*% H0 %*% r1
  Prop <- cbind(Prop, max(0, c(1 - lambda / v_c)) * r1)
  
  P_d <- Sig_d
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r2 <- Q_meta %*% (alpha0 - alpha1)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda <- max(sum(diag(P_d %*% H_d)) - 2 * Re(eigen(HR_d %*% P_d %*% HR_d)$values[1]), 0)
  v_d <- t(r2) %*% H0 %*% r2
  Prop[, 1] <- Prop[, 1] +  max(0, c(1 - lambda / v_d)) * r2
  
  A <- cbind(rep(1, p + 1), Q_meta %*% rep(1, q + 1))
  Pi <- diag(1, p + 1) - A %*% solve(t(A) %*% solve(Sig_c) %*% A) %*% t(A) %*% solve(Sig_c)
  P_c <- Pi %*% Sig_c %*% t(Pi)
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda <- max(sum(diag(P_c %*% H0)) - 2 * Re(eigen(HR0 %*% P_c %*% HR0)$values[1]), 0)
  mu_meta <- (diag(1, p + 1) - Pi) %*% meta
  r1 <- meta - mu_meta
  v_c <- t(r1) %*% H0 %*% r1
  Prop <- cbind(Prop, mu_meta + max(0, c(1 - lambda / v_c)) * r1)
  
  Pi <- diag(1, q + 1) - matrix(1, q + 1, q + 1) %*% solve(Sig_d) / sum(solve(Sig_d))
  P_d <- Pi %*% Sig_d %*% t(Pi)
  mu_d <- (diag(1, q + 1) - Pi) %*% (alpha0 - alpha1)
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r2 <- Q_meta %*% (alpha0 - alpha1 - mu_d)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda <- max(sum(diag(P_d %*% H_d)) - 2 * Re(eigen(HR_d %*% P_d %*% HR_d)$values[1]), 0)
  v_d <- t(r2) %*% H0 %*% r2
  Prop[, 2] <- Prop[, 2] + (Q_meta %*% mu_d + max(0, c(1 - lambda / v_d)) * r2)
  
  #############without external variances
  Sig_d <- Sig0a
  Q_meta <- Sig0ba %*% solve(Sig_d)
  meta_dag <- beta0 - Q_meta %*% (alpha0 - alpha1)
  
  Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
  P_c <- Sig_c
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda <- max(sum(diag(P_c %*% H0)) - 2 * Re(eigen(HR0 %*% P_c %*% HR0)$values[1]), 0)
  r1 <- meta_dag
  v_c <- t(r1) %*% H0 %*% r1
  Prop <- cbind(Prop, max(0, c(1 - lambda / v_c)) * r1)
  
  P_d <- Sig_d
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r2 <- Q_meta %*% (alpha0 - alpha1)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda <- max(sum(diag(P_d %*% H_d)) - 2 * Re(eigen(HR_d %*% P_d %*% HR_d)$values[1]), 0)
  v_d <- t(r2) %*% H0 %*% r2
  Prop[, 3] <- Prop[, 3] + max(0, c(1 - lambda / v_d)) * r2
  
  A <- cbind(rep(1, p + 1), Q_meta %*% rep(1, q + 1))
  Pi <- diag(1, p + 1) - A %*% solve(t(A) %*% solve(Sig_c) %*% A) %*% t(A) %*% solve(Sig_c)
  P_c <- Pi %*% Sig_c %*% t(Pi)
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda <- max(sum(diag(P_c %*% H0)) - 2 * Re(eigen(HR0 %*% P_c %*% HR0)$values[1]), 0)
  mu_meta <- (diag(1, p + 1) - Pi) %*% meta_dag
  r1 <- meta_dag - mu_meta
  v_c <- t(r1) %*% H0 %*% r1
  Prop <- cbind(Prop, mu_meta + max(0, c(1 - lambda / v_c)) * r1)
  
  Pi <- diag(1, q + 1) - matrix(1, q + 1, q + 1) %*% solve(Sig_d) / sum(solve(Sig_d))
  P_d <- Pi %*% Sig_d %*% t(Pi)
  mu_d <- (diag(1, q + 1) - Pi) %*% (alpha0 - alpha1)
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r2 <- Q_meta %*% (alpha0 - alpha1 - mu_d)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda <- max(sum(diag(P_d %*% H_d)) - 2 * Re(eigen(HR_d %*% P_d %*% HR_d)$values[1]), 0)
  v_d <- t(r2) %*% H0 %*% r2
  Prop[, 4] <- Prop[, 4] + (Q_meta %*% mu_d + max(0, c(1 - lambda / v_d)) * r2)
  
  #################Joint risk minimization
  Sig_d <- Sig0a + Sig1a
  Q_meta <- Sig0ba %*% solve(Sig_d)
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r1 <- meta
  r2 <- Q_meta %*% (alpha1 - alpha0)
  r12 <- cbind(r1, r2)
  B <- t(r12) %*% H0 %*% r12
  lambda_star <- c(sum(diag(Sig_c %*% H0)), - sum(diag(Sig_d %*% H_d)))
  weight <- solve(B) %*% lambda_star
  JointS <- r12 %*% (c(1, -1) - weight)
  
  meta_or <- beta0 - Q_meta_or %*% (alpha0 - alpha1)
  r1 <- meta_or
  r2 <- Q_meta_or %*% (alpha1 - alpha0)
  r12 <- cbind(r1, r2)
  comb_or <- r12 %*% (c(1, -1) - w_or)
  
  loss_Prop <- c()
  
  loss_beta0 <- t(beta0 - beta) %*% H0 %*% (beta0 - beta)
  loss_meta <- t(meta - beta) %*% H0 %*% (meta - beta)
  loss_comp_test <- t(comp_test - beta) %*% H0 %*% (comp_test - beta)
  loss_global_test <- t(global_test - beta) %*% H0 %*% (global_test - beta)
  loss_Prop[1] <- t(Prop[, 1] - beta) %*% H0 %*% (Prop[, 1] - beta)
  loss_Prop[2] <- t(Prop[, 2] - beta) %*% H0 %*% (Prop[, 2] - beta)
  loss_Prop[3] <- t(Prop[, 3] - beta) %*% H0 %*% (Prop[, 3] - beta)
  loss_Prop[4] <- t(Prop[, 4] - beta) %*% H0 %*% (Prop[, 4] - beta)
  loss_JointS <- t(JointS - beta) %*% H0 %*% (JointS - beta)
  loss_or <- t(comb_or - beta) %*% H0 %*% (comb_or - beta)
  
  c(loss_Prop, loss_beta0, loss_meta, loss_comp_test, loss_global_test, loss_JointS, loss_or) 
}

p <- 30
q <- floor(2 / 3 * p)
n <- 300
m <- 1000

set.seed(0)
sig_x0 <- sqrt(rgamma(p, 2, scale = 0.5))
Sig0 <- diag(sig_x0) %*% ar1_cor(p, 0.4) %*% diag(sig_x0)
sig_x <-sqrt(rgamma(p, 2, scale = 0.5))
Sigx <- diag(sig_x) %*% ar1_cor(p, 0.2) %*% diag(sig_x)
sim <- 5000

set.seed(0)
beta <- rnorm(p + 1, 0, sqrt(0.5 / (p + 1))) # runif(p, -0.2, 0.5) #rnorm(p, 0, sqrt(0.5/p)) #sqrt(0.5/p)
A_seq <- seq(0, 1, 0.1)

set.seed(0)
eta <- rbinom(p + 1, 1, 0.5) * rnorm(p + 1, sqrt(1 / (p + 1)), sqrt(0.5 / (p + 1))) #runif(p, - sqrt(3), sqrt(3)) 
Res_seq <- matrix(NA, 10, length(A_seq))

library(parallel)
library(MASS)


for (i in 1:length(A_seq)) {
  b <-  A_seq[i] * eta
  Sig1 <- (1 - A_seq[i]) * Sig0 + A_seq[i] * Sigx
  
  Sig1_aug <- rbind(c(1, rep(0, p)), cbind(0, Sig1))
  
  alpha1 <-   solve(Sig1_aug[1:(q + 1), 1:(q + 1)]) %*% Sig1_aug[1:(q + 1), ] %*% (beta + b)
  Sig1a <-  c(t((beta + b)[(q + 2):(p + 1)]) %*% (Sig1_aug[(q + 1):p, (q + 2):(p + 1)] - Sig1_aug[(q + 2):(p + 1), 1:(q + 1)] %*% 
              solve(Sig1_aug[1:(q + 1), 1:(q + 1)]) %*% Sig1_aug[1:(q + 1), (q + 2):(p + 1)]) %*%
               (beta + b)[(q + 2):(p + 1)] + 1) * solve(Sig1_aug[1:(q + 1), 1:(q + 1)]) / m
  
  
  Sig0_aug <- rbind(c(1, rep(0, p)), cbind(0, Sig0))
  beta0 <- beta
  alpha0 <- solve(Sig0_aug[1:(q + 1), 1:(q + 1)]) %*% Sig0_aug[1:(q + 1), ] %*% beta
  Sig0b <- solve(Sig0_aug) / n
  Sig0ba <- solve(Sig0_aug) %*% Sig0_aug[, 1:(q + 1)] %*% solve(Sig0_aug[1:(q + 1), 1:(q + 1)]) / n
  Sig0a <- c(t((beta + b)[(q + 2):(p + 1)]) %*% (Sig0_aug[(q + 2):(p + 1), (q + 2):(p + 1)] - Sig0_aug[(q + 1):p, 1:(q + 1)] %*% 
             solve(Sig0_aug[1:(q + 1), 1:(q + 1)]) %*% Sig0_aug[1:(q + 1), (q + 2):(p + 1)]) %*%
               (beta + b)[(q + 2):(p + 1)] + 1) * solve(Sig0_aug[1:(q + 1), 1:(q + 1)]) / n
  
  Sig_d <- Sig0a + Sig1a
  
  Q_meta_or <- Sig0ba %*% solve(Sig_d)
  meta <- beta0 - Q_meta_or %*% (alpha0 - alpha1)
  
  Sig_c <- Sig0b - Q_meta_or %*% (Sig_d) %*% t(Q_meta_or)
  H0 <- Sig0_aug
  r1 <- meta
  
  H_d <- t(Q_meta_or) %*% H0 %*% Q_meta_or
  r2 <- Q_meta_or %*% (alpha1 - alpha0)
  r12 <- cbind(r1, r2)
  B <- t(r12) %*% H0 %*% r12
  B <- B + diag(c(sum(diag(Sig_c %*% H0)), sum(diag(Sig_d %*% H_d))))
  w_or <- solve(B) %*% c(sum(diag(Sig_c %*% H0)), - sum(diag(Sig_d %*% H_d)))
  
  ##############compute the estimators
  
  cl <- makeCluster(10)
  clusterEvalQ(cl, c(library(MASS)))
  clusterExport(cl, c("beta", "Sig0", "Sig1", "p", "n", "m", "b", "q", "Q_meta_or",  "w_or"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, i] <- apply(Res, 1, mean)
}

library(latex2exp)

library(ggplot2)

Error <- data.frame(rep(A_seq, 5),
                    as.vector(t(Res_seq[c(5, 6, 8, 9, 1), ])),
                    rep(c("T", "C", "gt", "SURE", "ds"), rep(length(A_seq), 5)))

colnames(Error) <- c("t", "Error", "label")

Error$label <- factor(Error$label, levels= c("T", "C", "gt", "SURE", "ds"))


color <- c("T" = "green4", "C" = "blue", "gt" =  "darkgoldenrod1", "SURE" = "red", "ds" = "purple")

ggplot(Error) + geom_line(aes(t, Error, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, Error, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "Calibration", "Pretest", "SURE", "dShrink")) +
  coord_cartesian(ylim=c(0.05, 0.15)) +
  scale_shape_manual(values = c(1, 0, 15, 16, 17), labels = c("Target Pop", "Calibration", "Pretest", "SURE", "dShrink")) +
  labs(y = "MSE", x = "t") +
  theme_light() +
  theme(axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        legend.title = element_blank(),                   
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=30),
        legend.key.size=unit(1.5,'cm'))

ggsave(paste(paste("Error", p, n, m, "ds", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 8, units = "in")




Error <- data.frame(rep(A_seq, 3),
                    as.vector(t(Res_seq[c(5, 1, 2), ])),
                    rep(c("T", "ds", "mds"), rep(length(A_seq), 3)))

colnames(Error) <- c("t", "Error", "label")

Error$label <- factor(Error$label, levels= c("T", "ds", "mds"))


color <- c("T" = "green4", "ds" = "purple", "mds" = "red")

ggplot(Error) + geom_line(aes(t, Error, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, Error, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "dShink", "m-dShrink")) +
  coord_cartesian(ylim=c(0.05, 0.15)) +
  scale_shape_manual(values = c(1, 17, 15), labels = c("Target Pop", "dShink", "m-dShrink")) +
  labs(y = "MSE", x = "t") +
  theme_light() +
  theme(axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        legend.title = element_blank(),                   
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=30),
        legend.key.size=unit(1.5,'cm'))

ggsave(paste(paste("Error", p, n, m, "mds", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 8, units = "in")


Error <- data.frame(rep(A_seq, 3),
                    as.vector(t(Res_seq[c(5, 1, 3), ])),
                    rep(c("T", "ds", "cds"), rep(length(A_seq), 3)))

colnames(Error) <- c("t", "Error", "label")

Error$label <- factor(Error$label, levels= c("T", "ds", "cds"))


color <- c("T" = "green4", "ds" = "purple", "cds" = "red")

ggplot(Error) + geom_line(aes(t, Error, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, Error, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "dShink", "c-dShrink")) +
  coord_cartesian(ylim=c(0.05, 0.15)) +
  scale_shape_manual(values = c(1, 17, 15), labels = c("Target Pop", "dShink", "c-dShrink")) +
  labs(y = "MSE", x = "t") +
  theme_light() +
  theme(axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        legend.title = element_blank(),                   
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=30),
        legend.key.size=unit(1.5,'cm'))

ggsave(paste(paste("Error", p, n, m, "cds", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 8, units = "in")

