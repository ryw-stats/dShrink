
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
  
  dat1 <- data.frame(y, x)
  
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
  
  Pi_d <- diag(1, q + 1) #- matrix(1, q + 1, q + 1) %*% solve(Sig_d) / sum(solve(Sig_d))
  P_d <- Pi_d %*% Sig_d %*% t(Pi_d)
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r2 <- Q_meta %*% (alpha0 - alpha1)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda_d <- max(sum(diag(P_d %*% H_d)) - 2 * Re(eigen(HR_d %*% P_d %*% HR_d)$values[1]), 0)
  v_d <- t(r2) %*% H0 %*% r2
  Prop <- Prop + max(0, c(1 - lambda_d / v_d)) * r2
  
  
  ###bootstrap
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
    S <- Sig_d + s * diag(1, q + 1)
    t(alpha0 - alpha1) %*% ginv(S) %*% (alpha0 - alpha1) + log(det(S))
  }
  
  
  s2_c <- nlminb(1, Lik1, lower = 0)$par
  s2_d <- nlminb(1, Lik2, lower = 0)$par
  
  meta_pos <- s2_c * ginv(Sig_c + s2_c * diag(1, p + 1)) %*% meta
  Sig_c_pos <- s2_c * diag(1, p + 1) - s2_c^{2} * ginv(Sig_c + s2_c * diag(1, p + 1))
  
  d_pos <- s2_d * ginv(Sig_d + s2_d * diag(1, q + 1)) %*% (alpha0 - alpha1)
  Sig_d_pos <- s2_d * diag(1, q + 1) - s2_d^{2} * ginv(Sig_d + s2_d * diag(1, q + 1))
  
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
  
  cover <- c(abs(beta0 - beta) <= cri1, ((Prop - beta) >= cri2[, 1]) & ((Prop - beta) <= cri2[, 2]))
  c(cover, 2 * cri1, cri2[, 2] - cri2[, 1]) 
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
sim <- 5000

set.seed(0)
beta <- rnorm(p + 1, 0, sqrt(0.5/ (p + 1))) # runif(p, -0.2, 0.5) #rnorm(p, 0, sqrt(0.5/p)) #sqrt(0.5/p)
A_seq <- seq(0, 1, 0.1)

set.seed(0)
eta <- rbinom(p + 1, 1, 0.5) * rnorm(p + 1, sqrt(1 / (p + 1)), sqrt(0.5 / (p + 1))) #runif(p, - sqrt(3), sqrt(3)) 
Res_seq <- matrix(NA, 4 * (p + 1), length(A_seq))

library(parallel)
library(MASS)


for (i in 1:length(A_seq)) {
  b <-  A_seq[i] * eta
  Sig1 <- (1 - A_seq[i]) * Sig0 + A_seq[i] * Sigx
  
  cl <- makeCluster(20)
  clusterEvalQ(cl, c(library(MASS), library(sphunif)))
  clusterExport(cl, c("beta", "Sig0", "Sig1", "p", "n", "m", "b", "q"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, i] <- apply(Res, 1, mean)
}


library(latex2exp)
library(ggplot2)



CI <- data.frame(rep(A_seq, 2),
                 c(colMeans(Res_seq[1:(p + 1), ]), colMeans(Res_seq[(p + 2):(2*p + 2), ])),
                 rep(c("T", "ds"), rep(length(A_seq), 2)))

colnames(CI) <- c("t", "cover", "label")

CI$label <- factor(CI$label, levels= c("T", "ds"))


color <- c("T" = "green4", "ds" = "purple")

p1 <- ggplot(CI) + geom_line(aes(t, cover, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, cover, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "dShrink")) +
  coord_cartesian(ylim=c(0.8, 1)) +
  scale_shape_manual(values = c(1, 17), labels = c("Target Pop","dShrink")) +
  labs(y = "Coverage Rate", x = "t") +
  geom_hline(yintercept = 0.95, linewidth = 1.2, linetype = 3, color = "darkgrey") +
  theme_light() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=20),
        legend.key.size=unit(1.5,'cm'))

CI <- data.frame(rep(A_seq, 1), 10 * c(colMeans(Res_seq[(2 * p + 3):(3 * p + 3), ]),
                                       colMeans(Res_seq[(3 * p + 4):(4 * p + 4), ])), rep(c("T", "ds"), rep(length(A_seq), 2)))

colnames(CI) <- c("t", "width", "label")

CI$label <- factor(CI$label, levels= c("T", "ds"))


p2 <- ggplot(CI) + geom_line(aes(t, width, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, width, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "dShrink")) +
  scale_shape_manual(values = c(1, 17), labels = c("Target Pop","dShrink")) +
  labs(y = "Width", x = "t") +
  theme_light() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=20),
        legend.key.size=unit(1.5,'cm'))

library(ggpubr)
ggarrange(p1, p2, common.legend = T)

ggsave(paste(paste("CI", p, n, m, "ds", "comp", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 6, units = "in")
