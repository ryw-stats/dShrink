##bias is not zero unless they are both zero
##componentwise selection works when the bias is sparse
##EB performs very well and is robust against prior misspecification and bias level when p is large (over 300)
##EB can be explained as a finite population least squares estimator
##double shrink estimator or linear Empirical bayes estimator
##linear EB estimator performs better when the sample size is large
##double shrink outperform EB when the sample size is around 400
##linear EB outperforms EB when the sample size is larger than 800
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
  
  H0 <- diag(1, p)
  
  Q_meta <- Sig0ba %*% solve(Sig0a)
  
  PC_H0 <- eigen(H0)
  HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
  lambda <- max(sum(diag(Sig0b %*% H0)) - 2 * Re(eigen(HR0 %*% Sig0b %*% HR0)$values[1]), 0)
  v_JS <- t(beta0 - alpha1) %*% H0 %*% (beta0 - alpha1)
  JS <- alpha1 + max(0, c(1 - lambda / v_JS)) * (beta0 - alpha1)
  
  
  w_delr <- sum((alpha1 - alpha0)^2) / (sum((alpha1 - alpha0)^2) + 2 * p / n)
  
  delr <- w_delr * beta0 + (1 - w_delr) * alpha1 
  
  Q_meta <- Sig0ba %*% solve(Sig0a + Sig1a)
  meta <- beta0 - Q_meta %*% (alpha0 - alpha1)
  
  Sig_c <- Sig0b - Q_meta %*% (Sig0a + Sig1a) %*% t(Q_meta)
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
  

  y0 <- U0 %*% beta + rnorm(n, 0, sqrt(2))
  y1 <- U1 %*% (mu - b) + rnorm(m, 0, sqrt(2))
  
  
  loss_beta0 <- t(beta0 - beta) %*% H0 %*% (beta0 - beta)
  loss_JS <- t(JS - beta) %*% H0 %*% (JS - beta)
  loss_delr <- t(delr - beta) %*% H0 %*% (delr - beta)
  loss_Prop <- t(Prop - beta) %*% H0 %*% (Prop - beta)

  c(loss_beta0, loss_JS, loss_delr, loss_Prop) 
}

p <- 5
n <- 300
m <- 300

library("mclust")
U0 <- randomOrthogonalMatrix(n, p) * sqrt(n)
U1 <- randomOrthogonalMatrix(m, p) * sqrt(m)

set.seed(0)
sim <- 5000

mu <- sqrt(1 / p) * c(0.05, 0.02, 0.1, 0.1, 0.1) 
A_seq <- seq(0, 1, 0.1)

eta <-  c(0.2, 0.3, 0.3, 0.3, 0.3) * sqrt(1 / p)
Res_seq <- matrix(NA, 4, length(A_seq))

library(parallel)
library(MASS)

for (i in 1:length(A_seq)) {
  b <-  A_seq[i] * eta 
  
  ##############compute the estimators
  
  cl <- makeCluster(10)
  clusterEvalQ(cl, c(library(MASS), library(ISEDI)))
  clusterExport(cl, c("mu", "p", "n", "m", "b", "U0", "U1"))
  Res <- parSapply(cl, 1:sim, simone)
  stopCluster(cl)
  Res_seq[, i] <- apply(Res, 1, mean)
}

library(latex2exp)



library(ggplot2)

Error <- data.frame(rep(A_seq, 4),
                    as.vector(t(Res_seq[, ])),
                    rep(c("T", "JS Combine", "Data Enriched", "ds"), rep(length(A_seq), 4)))

colnames(Error) <- c("t", "Error", "label")

Error$label <- factor(Error$label, levels= c("T", "JS Combine", "Data Enriched", "ds"))


color <- c("T" = "green4", "JS Combine" =  "red", "Data Enriched" = "deepskyblue3", "ds" = "purple")

ggplot(Error) + geom_line(aes(t, Error, color = label), linewidth = 1.2, linetype = 2) +
  geom_point(aes(t, Error, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "JS Combine", "Data Enriched", "dShrink")) +
  coord_cartesian(ylim=c(0.013, 0.04)) +
  scale_shape_manual(values = c(1, 15, 0, 17), labels = c("Target Pop", "JS Combine", "Data Enriched", "dShrink")) +
  labs(y = "MSE", x = "t") +
  theme_light() +
  theme(axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        legend.title = element_blank(),                    
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=30),
        legend.key.size=unit(1.5,'cm'))

ggsave(paste(paste("Eg", p, "Stein", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 8, units = "in")

