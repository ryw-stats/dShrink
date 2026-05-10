setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


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
  
  
  Sig_c <- Sig0b - Q_meta %*% (Sig0a + Sig1a) %*% t(Q_meta)
  H0 <- diag(1, p)
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda <- max(sum(diag(Sig_c %*% H0)) - 2 * Re(eigen(HR0 %*% Sig_c %*% HR0)$values[1]), 0)
  v_c <- t(meta) %*% H0 %*% meta
  Prop <- max(0, c(1 - lambda / v_c)) * meta
  
  Sig_d <- Sig0a + Sig1a
  H_d <- t(Q_meta) %*% H0 %*% Q_meta
  PC_H_d <- eigen(H_d)
  HR_d <- PC_H_d$vectors %*% diag(sqrt(PC_H_d$values - PC_H_d$values * (PC_H_d$values < 0))) %*% t(PC_H_d$vectors)
  lambda <- max(sum(diag(Sig_d %*% H_d)) - 2 * Re(eigen(HR_d %*% Sig_d %*% HR_d)$values[1]), 0)
  v_d <- t(alpha0 - alpha1) %*% H_d %*% (alpha0 - alpha1)
  Prop <- Prop + Q_meta %*% (max(0, c(1 - lambda / v_d)) * (alpha0 - alpha1))
  
  
  r2 <- Q_meta %*% (alpha1 - alpha0)
  r12 <- cbind(meta, r2)
  B <- t(r12) %*% H0 %*% r12
  lambda_star <- c(sum(diag(Sig_c %*% H0)), - sum(diag(Sig_d %*% H_d)))
  weight <- solve(B) %*% lambda_star
  JointS <- r12 %*% (c(1, -1) - weight)
  
  r2 <- Q_meta %*% (alpha1 - alpha0)
  r12 <- cbind(meta, r2)
  B <- t(r12) %*% H0 %*% r12 + 0.1 * diag(1, 2)
  lambda_star <- c(sum(diag(Sig_c %*% H0)), - sum(diag(Sig_d %*% H_d)))
  weight <- solve(B) %*% lambda_star
  JointS1 <- r12 %*% (c(1, -1) - weight)
  
  r2 <- Q_meta %*% (alpha1 - alpha0)
  r12 <- cbind(meta, r2)
  B <- t(r12) %*% H0 %*% r12 + 0.005 * diag(1, 2)
  lambda_star <- c(sum(diag(Sig_c %*% H0)), - sum(diag(Sig_d %*% H_d)))
  weight <- solve(B) %*% lambda_star
  JointS2 <- r12 %*% (c(1, -1) - weight)
  
  r2 <- Q_meta %*% (alpha1 - alpha0)
  r12 <- cbind(meta, r2)
  B <- t(r12) %*% H0 %*% r12 + 0.0005 * diag(1, 2)
  lambda_star <- c(sum(diag(Sig_c %*% H0)), - sum(diag(Sig_d %*% H_d)))
  weight <- solve(B) %*% lambda_star
  JointS3 <- r12 %*% (c(1, -1) - weight)
  
  loss_beta0 <- t(beta0 - beta) %*% H0 %*% (beta0 - beta)
  loss_Prop <- t(Prop - beta) %*% H0 %*% (Prop - beta)
  loss_JointS <- t(JointS - beta) %*% H0 %*% (JointS - beta)
  loss_JointS1 <- t(JointS1 - beta) %*% H0 %*% (JointS1 - beta)
  loss_JointS2 <- t(JointS2 - beta) %*% H0 %*% (JointS2 - beta)
  loss_JointS3 <- t(JointS3 - beta) %*% H0 %*% (JointS3 - beta)
  

  c(loss_beta0, loss_Prop, loss_JointS, loss_JointS1, loss_JointS2, loss_JointS3) 
}

p <- 5
n <- 300
m <- 300

set.seed(0)
sim <- 5000

mu <- sqrt(1 / p) * c(0.05, 0.02, 0.1, 0.1, 0.1) 
A_seq <- seq(0, 1, 0.1)

eta <-  c(0.2, 0.3, 0.3, 0.3, 0.3) * sqrt(1 / p)
Res_seq <- matrix(NA, 6, length(A_seq))

library(parallel)
library(MASS)

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
  
  ##############compute the estimators
  
  cl <- makeCluster(10)
  clusterEvalQ(cl, c(library(MASS)))
  clusterExport(cl, c("mu", "p", "n", "m", "b", "w_or"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, i] <- apply(Res, 1, mean)
}

library(latex2exp)

# pdf(paste(p, "-", n, "-", m, "-DS",".pdf", sep = ""), width = 8, height = 6, onefile = FALSE)
# plot(A_seq, Res_seq[6, ], type = "l", col = 4, ylim = c(0.08, 0.125), lwd = 1.8, lty = 2, xlab = TeX('$delta$'), ylab = "Error", cex.lab = 1.2)
# points(A_seq, Res_seq[6, ], col = 4, cex = 1.2, pch = 0)
# lines(A_seq, Res_seq[5, ], type = "l", col = 3, lwd = 1.8, lty = 2)
# points(A_seq, Res_seq[5, ], col = 3, cex = 1.2, pch = 1)
# lines(A_seq, Res_seq[7, ], type = "l", col = 5, lwd = 1.8, lty = 2)
# points(A_seq, Res_seq[7, ], col = 5, cex = 1.2, pch = 2)
# lines(A_seq, Res_seq[8, ], type = "l", col = 2, lwd = 1.8, lty = 2)
# points(A_seq, Res_seq[8, ], col = 2, cex = 1.2, pch = 15)
# lines(A_seq, Res_seq[9, ], type = "l", col = 7, lwd = 1.8, lty = 2)
# points(A_seq, Res_seq[9, ], col = 7, cex = 1.2, pch = 18)
# lines(A_seq, Res_seq[10, ], type = "l", col = 8, lwd = 1.8, lty = 2)
# points(A_seq, Res_seq[10, ], col = 8, cex = 1.2, pch = 19)
# lines(A_seq, Res_seq[1, ], type = "l", col = 6, lwd = 1.8, lty = 2)
# points(A_seq, Res_seq[1, ], col = 6, cex = 1.2, pch = 17)
# dev.off()


library(ggplot2)


Error <- data.frame(rep(A_seq, 6),
                    as.vector(t(Res_seq)),
                    rep(c("T", "dShrink", "SURE", "ridge1", "ridge2", "ridge3"), rep(length(A_seq), 6)))

colnames(Error) <- c("t", "Error", "label")

ridge1 <- unname(TeX("\\bf{Ridge} $(r = 0.1)$"))
ridge2 <- unname(TeX("\\bf{Ridge} $(r = 5\\times 10^{-3})$"))
ridge3 <- unname(TeX("\\bf{Ridge} $(r = 5\\times 10^{-4})$"))

Error$label <- factor(Error$label, levels= c("T", "ridge1", "ridge2", "ridge3", "SURE", "dShrink"))


color <- c("T" = "green4", "ridge1" = "blue", "ridge2" = "deepskyblue", "ridge3" = "darkgoldenrod1", "SURE" = "red",  "dShrink" = "purple")

ggplot(Error) + geom_line(aes(t, Error, color = label), size = 1.5, linetype = 2) +
  geom_point(aes(t, Error, color = label, shape = label), size = 5) +
  scale_color_manual(values = color, labels = c("Target Pop", ridge1, ridge2, ridge3, "SURE", "dShrink")) +
  coord_cartesian(ylim=c(0.012, 0.05)) +
  scale_shape_manual(values = c(1, 0, 2, 15, 16, 17), labels = c("Target Pop", ridge1, ridge2, ridge3, "SURE", "dShrink")) +
  labs(y = "MSE", x = "t") +
  theme_light() +
  theme(legend.text.align = 0) +
  theme(axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        legend.title = element_blank(),                    
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=30),
        legend.key.size=unit(1.5,'cm'))

ggsave(paste(paste("Eg", p, "ridge", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 8, units = "in")

