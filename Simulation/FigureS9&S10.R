setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

simone <- function(sim) {
  set.seed(sim)
  
  x <- matrix(rnorm(m * p), m)
  y <- plogis(x %*% (mu + t * eta) / sqrt(sum((mu + t * eta)^2))) %*% t(rep(1, 3)) + 0.5 * rt(m, 10) %*% t(rep(1, 3)) + 
       0.5 * matrix((2 * rbinom(3 * m, 1, 0.5) - 1) * rexp(3 * m, sqrt(2)), m)
  
  GEE <- function(par) {
    sum((y - plogis(x %*% par) %*% t(rep(1, 3)))^{2})
  }
  
  init1 <- nlminb(rep(0, p), GEE)$par
  
  e1 <- y - plogis(x %*% init1) %*% t(rep(1, 3))
  
  V1 <- solve(cov(e1))
  
  GEE <- function(par) {
    g <- colMeans(c((y - plogis(x %*% par) %*% t(rep(1, 3))) %*% V1 %*% rep(1, 3)) * (c(plogis(x %*% init1) * (1 - plogis(x %*% init1))) * x))
    sum(g^2)
  }
  
  alpha1 <- nlminb(rep(0, p), GEE)$par
  
  G <- solve(sum(V1) * (t(x) %*% (c(plogis(x %*% alpha1) * (1 - plogis(x %*% alpha1)))^{2} * x)) / m)
  
  Sig1a <- cov(c((y - plogis(x %*% alpha1) %*% t(rep(1, 3))) %*% V1 %*% rep(1, 3)) * (c(plogis(x %*% alpha1) * (1 - plogis(x %*% alpha1))) * x) %*% G) / m
  
  beta <- (mu - t * eta) / sqrt(sum((mu - t * eta)^2))
  
  x <- matrix(rnorm(n * p), n)
  y <- plogis(x %*% beta) %*% t(rep(1, 3)) + 0.5 * rt(n, 10) %*% t(rep(1, 3)) + 
    0.5 * matrix((2 * rbinom(3 * n, 1, 0.5) - 1) * rexp(3 * n, sqrt(2)), n)
  
  GEE <- function(par) {
    sum((y - plogis(x %*% par) %*% t(rep(1, 3)))^{2})
  }
  
  init0 <- nlminb(rep(0, p), GEE)$par
  
  e0 <- y - plogis(x %*% init0) %*% t(rep(1, 3))
  
  V0 <- solve(cov(e0))
  
  GEE <- function(par) {
    g <- colMeans(c((y - plogis(x %*% par) %*% t(rep(1, 3))) %*% V0 %*% rep(1, 3)) * (c(plogis(x %*% init0) * (1 - plogis(x %*% init0))) * x))
    sum(g^2)
  }
  
  beta0 <- nlminb(rep(0, p), GEE)$par
  
  G <- solve(sum(V1) * (t(x) %*% (c(plogis(x %*% beta0) * (1 - plogis(x %*% beta0)))^{2} * x)) / n)
  
  Sig0b <- cov(c((y - plogis(x %*% beta0) %*% t(rep(1, 3))) %*% V0 %*% rep(1, 3)) * (c(plogis(x %*% beta0) * (1 - plogis(x %*% beta0))) * x) %*% G) / n
  
  alpha0 <- beta0
  
  Sig0a <- Sig0b
  Sig0ba <- Sig0b
  
  
  Q_meta <- Sig0ba %*% solve(Sig0a + Sig1a)
  meta <- beta0 - Q_meta %*% (alpha0 - alpha1)
  
  test <- c(t(alpha0 - alpha1) %*% solve(Sig0a + Sig1a) %*% (alpha0 - alpha1)  < qchisq(0.95, p))
  
  global_test <- beta0 + (meta - beta0) * test
  
  H0 <- diag(1, p)
  
  Sig_d <- Sig0a + Sig1a
  Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda <- max(sum(diag(Sig_c %*% H0)) - 2 * Re(eigen(HR0 %*% Sig_c %*% HR0)$values[1]), 0)
  r1 <- meta
  v_c <- t(r1) %*% H0 %*% r1
  Prop <- max(0, c(1 - lambda / v_c)) * r1
  
  Sig_d <- Sig_d
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r2 <- Q_meta %*% (alpha0 - alpha1)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda <- max(sum(diag(Sig_d %*% H_d)) - 2 * Re(eigen(HR_d %*% Sig_d %*% HR_d)$values[1]), 0)
  v_d <- t(r2) %*% H0 %*% r2
  Prop <- Prop +  max(0, c(1 - lambda / v_d)) * r2
  
  
  Sig_d <- Sig0a + Sig1a
  Q_meta <- Sig0ba %*% solve(Sig_d)
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r1 <- meta
  r2 <- Q_meta %*% (alpha1 - alpha0)
  r12 <- cbind(r1, r2)
  B <- t(r12) %*% H0 %*% r12
  lambda_star <- c(sum(diag(Sig_c %*% H0)), - sum(diag(Sig_d %*% H_d)))
  weight <- ginv(B) %*% lambda_star
  JointS <- r12 %*% (c(1, -1) - weight)
  
  loss_beta0 <- t(beta0 - beta) %*% H0 %*% (beta0 - beta)
  loss_meta <- t(meta - beta) %*% H0 %*% (meta - beta)
  loss_global_test <- t(global_test - beta) %*% H0 %*% (global_test - beta)
  loss_Prop <- t(Prop - beta) %*% H0 %*% (Prop - beta)
  loss_JointS <- t(JointS - beta) %*% H0 %*% (JointS - beta)
  
  
  ##Inference
  H0 <- solve(Sig0b)
  
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda_c <- max(sum(diag(Sig_c %*% H0)) - 2 * Re(eigen(HR0 %*% Sig_c %*% HR0)$values[1]), 0)
  r1 <- meta
  v_c <- t(r1) %*% H0 %*% r1
  Prop_Inf <- max(0, c(1 - lambda_c / v_c)) * r1
  
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r2 <- Q_meta %*% (alpha0 - alpha1)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda_d <- max(sum(diag(Sig_d %*% H_d)) - 2 * Re(eigen(HR_d %*% Sig_d %*% HR_d)$values[1]), 0)
  v_d <- t(r2) %*% H0 %*% r2
  Prop_Inf <- Prop_Inf + max(0, c(1 - lambda_d / v_d)) * r2
  
  
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
  
  length_Prop <- mean(cri2[, 2] - cri2[, 1])
  length_beta0 <- mean(2 * cri1)
  
  cover_Prop <- mean(((Prop_Inf - beta) >= cri2[, 1]) & ((Prop_Inf - beta) <= cri2[, 2]))
  cover_beta0 <- mean(abs(beta0 - beta) <= cri1)
  
  
  c(loss_Prop, loss_beta0, loss_meta, loss_global_test, loss_JointS, length_Prop, length_beta0, cover_Prop, cover_beta0) 
}

p <- 6
n <- 300
m <- 600

sim <- 5000

set.seed(0)
mu <- rnorm(p, 0, sqrt(0.25 / p)) # runif(p, -0.2, 0.5) #rnorm(p, 0, sqrt(0.5/p)) #sqrt(0.5/p)
A_seq <- seq(0, 1, 0.1)

eta <- rnorm(p, 0, sqrt(0.25 / p)) #runif(p, - sqrt(3), sqrt(3)) 
Res_seq <- matrix(NA, 9, length(A_seq))

library(parallel)
library(MASS)


for (i in 1:length(A_seq)) {
  t <-  A_seq[i]
  
  cl <- makeCluster(20)
  clusterEvalQ(cl, c(library(MASS)))
  clusterExport(cl, c("mu", "p", "n", "m", "eta", "t"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, i] <- apply(Res, 1, mean)
}

library(latex2exp)

library(ggplot2)


Error <- data.frame(rep(A_seq, 5),
                    as.vector(t(Res_seq[c(2, 3, 4, 5, 1), ])),
                    rep(c("T", "C", "gt", "SURE", "ds"), rep(length(A_seq), 5)))

colnames(Error) <- c("t", "Error", "label")

Error$label <- factor(Error$label, levels= c("T", "C", "gt", "SURE", "ds"))


color <- c("T" = "green4", "C" = "blue", "gt" =  "darkgoldenrod1", "SURE" = "red", "ds" = "purple")

ggplot(Error) + geom_line(aes(t, Error, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, Error, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "Calibration", "Pretest", "SURE", "dShrink")) +
  coord_cartesian(ylim=c(0.05, 0.35)) +
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

ggsave(paste(paste("Error", p, n, m, "longitudinal", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 8, units = "in")


CI <- data.frame(rep(A_seq, 2),
                 c(Res_seq[9, ], Res_seq[8, ]),
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

CI <- data.frame(rep(A_seq, 1), c(Res_seq[7, ], Res_seq[6, ]), rep(c("T", "ds"), rep(length(A_seq), 2)))

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

ggsave(paste(paste("CI", p, n, m, "longitudinal", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 6, units = "in")


