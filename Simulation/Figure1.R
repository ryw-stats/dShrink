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
  
  sigma_comp <- sqrt(diag(Sig0a) + diag(Sig1a))
  
  test <- (abs(alpha0 - alpha1) < qnorm(0.975) * sigma_comp)
  comp_test <- beta0
  if (sum(test) > 0) {
    comp_test <- comp_test - Sig0ba[, test] %*% solve(Sig0a[test, test] + Sig1a[test, test]) %*% (alpha0[test] - alpha1[test])
  }
  
  test <- c(t(alpha0 - alpha1) %*% solve(Sig0a + Sig1a) %*% (alpha0 - alpha1)  < qchisq(0.95, p))
  
  global_test <- beta0 + (meta - beta0) * test
  
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

  loss_beta0 <- t(beta0 - beta) %*% H0 %*% (beta0 - beta)
  loss_meta <- t(meta - beta) %*% H0 %*% (meta - beta)
  loss_comp_test <- t(comp_test - beta) %*% H0 %*% (comp_test - beta)
  loss_global_test <- t(global_test - beta) %*% H0 %*% (global_test - beta)
  loss_Prop <- t(Prop - beta) %*% H0 %*% (Prop - beta)
  loss_JointS <- t(JointS - beta) %*% H0 %*% (JointS - beta)

  c(loss_beta0, loss_meta, loss_comp_test, loss_global_test, loss_Prop, loss_JointS) 
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
  
  
  cl <- makeCluster(10)
  clusterEvalQ(cl, c(library(MASS)))
  clusterExport(cl, c("mu", "p", "n", "m", "b"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, i] <- apply(Res, 1, mean)
}

library(latex2exp)

library(ggplot2)


Error <- data.frame(rep(A_seq, 5),
                    as.vector(t(Res_seq[c(1, 2, 4, 6, 5), ])),
                    rep(c("T", "C", "gt", "SURE", "ds"), rep(length(A_seq), 5)))

colnames(Error) <- c("t", "Error", "label")

Error$label <- factor(Error$label, levels= c("T", "C", "gt", "SURE", "ds"))


color <- c("T" = "green4", "C" = "blue", "gt" = "darkgoldenrod1", "SURE" = "red" , "ds" = "purple")

ggplot(Error) + geom_line(aes(t, Error, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, Error, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "Calibration", "Pretest", "SURE", "dShrink")) +
  coord_cartesian(ylim=c(0.012, 0.05)) +
  scale_shape_manual(values = c(1, 0, 15, 16, 17), labels = c("Target Pop", "Calibration", "Pretest", "SURE", "dShrink")) +
  labs(y = "Error", x = "t") +
  theme_light() +
  theme(axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        legend.title = element_blank(),                    
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=30),
        legend.key.size=unit(1.5,'cm'))

ggsave(paste(paste("Eg", p, "ds", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 8, units = "in")

