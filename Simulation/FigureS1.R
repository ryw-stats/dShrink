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
  
  beta <- mu + b
  alpha1 <- rnorm(p, mu - b, sqrt(2 / m))
  Sig1a <- diag(2, p) / m
  
  beta0 <- rnorm(p, beta, sqrt(2 / n))
  alpha0 <- beta0
  Sig0b <- diag(2, p) / n
  Sig0a <- Sig0b
  Sig0ba <- Sig0b
  
  Q_meta <- Sig0ba %*% solve(Sig0a + Sig1a)
  meta <- beta0 - Q_meta %*% (alpha0 - alpha1)
  
  w_sep <- c()
  
  H0 <- diag(1, p)
  
  Sig_d <- Sig0a + Sig1a
  Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda <- max(sum(diag(Sig_c %*% H0)) - 2 * Re(eigen(HR0 %*% Sig_c %*% HR0)$values[1]), 0)
  r1 <- meta
  v_c <- t(r1) %*% H0 %*% r1

  w_sep[1] <- max(0, c(1 - lambda / v_c))
  
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  r2 <- Q_meta %*% (alpha0 - alpha1)
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda <- max(sum(diag(Sig_d %*% H_d)) - 2 * Re(eigen(HR_d %*% Sig_d %*% HR_d)$values[1]), 0)
  v_d <- t(r2) %*% H0 %*% r2

  w_sep[2] <- max(0, c(1 - lambda / v_d))
  
  Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
  r1 <- meta
  r2 <- Q_meta %*% (alpha1 - alpha0)
  r12 <- cbind(r1, r2)
  B <- t(r12) %*% H0 %*% r12
  lambda_star <- c(sum(diag(Sig_c %*% H0)), - sum(diag(Sig_d %*% H_d)))
  weight <- solve(B) %*% lambda_star
  
  c(w_sep - w_or, weight - w_or) 
}

p <- 5
n <- 300
m <- 300

sim <- 5000

mu <- sqrt(1 / p) * c(0.05, 0.02, 0.1, 0.1, 0.1) 
A_seq <- seq(0, 1, 0.1)

eta <-  c(0.2, 0.3, 0.3, 0.3, 0.3) * sqrt(1 / p)
Res_seq <- matrix(NA, 4, length(A_seq))

library(parallel)
library(MASS)

cond <- c()

for (i in 1:length(A_seq)) {
  b <-  A_seq[i] * eta 
  
  beta0 <- alpha0 <- mu + b
  alpha1 <- mu - b
  
  Sig0b <- Sig0ba <- diag(2 / n, p)
  Sig_d <- diag(2 / n + 2 / m, p)
  
  Q_meta_or <- Sig0ba %*% solve(Sig_d)
  meta <- beta0 - Q_meta_or %*% (alpha0 - alpha1)
  
  Sig_c <- Sig0b - Q_meta_or %*% (Sig_d) %*% t(Q_meta_or)
  H0 <- diag(1, p)
  r1 <- meta
  
  H_d <- t(Q_meta_or) %*% H0 %*% Q_meta_or
  r2 <- Q_meta_or %*% (alpha1 - alpha0)
  r12 <- cbind(r1, r2)
  B <- t(r12) %*% H0 %*% r12
  B <- B + diag(c(sum(diag(Sig_c %*% H0)), sum(diag(Sig_d %*% H_d))))
  
  ev <- eigen(B)$values
  cond[i] <- ev[1] / ev[2]
}

for (i in 1:length(A_seq)) {
  b <-  A_seq[i] * eta 
  
  beta0 <- alpha0 <- mu + b
  alpha1 <- mu - b
  
  Sig0b <- Sig0ba <- diag(2 / n, p)
  Sig_d <- diag(2 / n + 2 / m, p)
  
  Q_meta_or <- Sig0ba %*% solve(Sig_d)
  meta <- beta0 - Q_meta_or %*% (alpha0 - alpha1)
  
  Sig_c <- Sig0b - Q_meta_or %*% (Sig_d) %*% t(Q_meta_or)
  H0 <- diag(1, p)
  r1 <- meta
  
  H_d <- t(Q_meta_or) %*% H0 %*% Q_meta_or
  r2 <- Q_meta_or %*% (alpha1 - alpha0)
  r12 <- cbind(r1, r2)
  B <- t(r12) %*% H0 %*% r12
  B <- B + diag(c(sum(diag(Sig_c %*% H0)), sum(diag(Sig_d %*% H_d))))
  w_or <- solve(B) %*% c(sum(diag(Sig_c %*% H0)), - sum(diag(Sig_d %*% H_d)))
  
  w_s <- c(t(meta) %*% H0 %*% (meta) / (sum(diag(Sig_c %*% H0)) + t(meta) %*% H0 %*% (meta)),
           t(alpha1 - alpha0) %*% H_d %*% (alpha1 - alpha0) / (sum(diag(Sig_d %*% H_d)) + t(alpha1 - alpha0) %*% H_d %*% (alpha1 - alpha0)))
  
  ##############compute the estimators
  
  cl <- makeCluster(10)
  clusterEvalQ(cl, c(library(MASS)))
  clusterExport(cl, c("mu", "p", "n", "m", "b", "w_or", "w_s"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, i] <- apply(Res^2, 1, mean)
}

library(latex2exp)
library(ggplot2)

Error <- data.frame(rep(A_seq, 2),
                    as.vector(t(Res_seq[c(3, 1), ])),
                    rep(c("SURE", "ds"), rep(length(A_seq), 2)))

colnames(Error) <- c("t", "Error", "label")

Error$label <- factor(Error$label, levels= c("SURE","ds"))


color <- c("SURE" = "red", "ds" = "purple")

p1 <- ggplot(Error) + geom_line(aes(t, Error, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, Error, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c(TeX('$widehat(lambda)_{SURE,1}$'), TeX('$widehat(lambda)_{s,1}$'))) +
  scale_shape_manual(values = c(16, 17), labels = c(TeX('$widehat(lambda)_{SURE,1}$'), TeX('$widehat(lambda)_{s,1}$'))) +
  labs(y = "MSE", x = "t") +
  theme_light() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.title = element_blank(),                   
        legend.text = element_text(size = 20, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=20),
        legend.key.size=unit(1.5,'cm'))

Error <- data.frame(rep(A_seq, 2),
                    as.vector(t(Res_seq[c(4, 2), ])),
                    rep(c("SURE", "ds"), rep(length(A_seq), 2)))

colnames(Error) <- c("t", "Error", "label")

Error$label <- factor(Error$label, levels= c("SURE","ds"))


color <- c("SURE" = "red", "ds" = "purple")

p2 <- ggplot(Error) + geom_line(aes(t, Error, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, Error, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c(TeX('$widehat(lambda)_{SURE,2}$'), TeX('$widehat(lambda)_{s,2}$'))) +
  scale_shape_manual(values = c(16, 17), labels = c(TeX('$widehat(lambda)_{SURE,2}$'), TeX('$widehat(lambda)_{s,2}$'))) +
  labs(y = "MSE", x = "t") +
  theme_light() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.title = element_blank(),                  
        legend.text = element_text(size = 20, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=20),
        legend.key.size=unit(1.5,'cm'))

library(ggpubr)
ggarrange(p1, p2)

ggsave(paste(paste("Weight", p, n, m, "wOPT", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 6, units = "in")

