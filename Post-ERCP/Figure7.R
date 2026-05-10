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


X <- as.matrix(DATA[, -c(1, 2, 31)])

RX <- cor(X)

r <- sqrt(1 - 1 / diag(solve(RX)) / diag(RX))

X_aug <- cbind(1, X[, r <= 0.95])


library(MASS)
library(sphunif)
alpha1 <- c(5/57, 6/47,#young
            16/214, 19/212,#female
            20/246, 20/217,#pdstent
            4/64, 6/52,#difcan
            4/35, 3/40#bsphinc
)

var_lst <- c("young", "gender", "pdstent", "difcan", "bsphinc")

y_tilde <- 2 * DATA$treatment * DATA$outcome - 2 * (1 - DATA$treatment) * DATA$outcome

n <- nrow(X_aug)
p <- ncol(X_aug)
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


#lambda_min <- cv.glmnet(X[, -1], y_tilde, family = "gaussian", lambda = seq(0.5, 1, 0.1) * sqrt(log(p) / n))$lambda.min

loss <- function(b) {
  mean((y_tilde - 2 * plogis(X_aug %*% b) + 1)^2)
}

beta0 <- nlminb(rep(0, p), loss, control = list(iter = 1000))$par

p_lin <- X_aug %*% beta0
h0b <- c(8 * plogis(p_lin)^2 * (1 - plogis(p_lin))^2)
G <- t(X_aug) %*% (h0b * X_aug) / n + diag(0.1, p) / n

IF_beta <- - (- 4 * c(plogis(p_lin) * (1 - plogis(p_lin)) * (y_tilde - 2 * plogis(p_lin) + 1)) * X_aug) %*% ginv(G)

Sig0b <- cov(IF_beta) / n
Sig0ba <-  cov(IF_beta, IF_alpha) / n

tmp_ev <- eigen(Sig0b)

H0 <- ginv(Sig0b)

H0 <- H0 + 0.1 * diag(diag(Sig0b)^{-1})

Q_meta <- Sig0ba %*% ginv(Sig0a)

meta_dag <- beta0 - Q_meta %*% (alpha0 - alpha1)

eval_loss <- c()

Sig_d <- Sig0a
H_d <- t(Q_meta) %*% H0 %*% Q_meta
Sig_c <- Sig0b - Q_meta %*% Sig_d %*% t(Q_meta)
A <- cbind(Sig_c %*% rep(1, p))
Pi_c <- diag(1, p) #- A %*% ginv(t(A) %*% ginv(Sig_c) %*% A) %*% t(A) %*% ginv(Sig_c)
P_c <- Pi_c %*% Sig_c %*% t(Pi_c)

PC_H0 <- eigen(H0)
HR0 <- PC_H0$vectors %*% diag(sqrt(PC_H0$values - PC_H0$values * (PC_H0$values < 0))) %*% t(PC_H0$vectors)
lambda_c <- max(sum(diag(P_c %*% H0)) - 2 * Re(eigen(HR0 %*% P_c %*% HR0)$values[1]), 0)
mu_meta <- (diag(1, p) - Pi_c) %*% meta_dag
r1 <- meta_dag - mu_meta
v_c <- t(r1) %*% H0 %*% r1
Prop <- mu_meta + max(1 - lambda_c / v_c, 0) * r1


A <- cbind(Sig_d %*% rep(1, q))
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
LB <- matrix(NA, n, BT)
UB <- matrix(NA, n, BT)

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
    
    Boots <- cbind(Boots, X_aug %*% (Prop_b - meta_tmp - Q_meta %*% d_tmp))
  }
  LB[, t] <- apply(Boots, 1, function(x) quantile(x, 0.025))
  UB[, t] <- apply(Boots, 1, function(x) quantile(x, 0.975))
}

sd <- sqrt(diag(X_aug %*% Sig0b %*% t(X_aug)))

EST_tar <- 2 * plogis(X_aug %*% beta0) - 1

EST_Prop <- 2 * plogis(X_aug %*% Prop) - 1

CI_tar <- cbind(2 * plogis(X_aug %*% beta0 - sd * qnorm(0.975)) - 1, 2 * plogis(X_aug %*% beta0 + sd * qnorm(0.975)) - 1)

CI_Prop <- cbind(2 * plogis(X_aug %*% Prop - apply(UB, 1, max)) - 1, 2 * plogis(X_aug %*% Prop - apply(LB, 1, min)) - 1)

round(c(mean(CI_tar[, 2] - CI_tar[, 1]),
        mean((CI_Prop[, 2] - CI_Prop[, 1]))), 2)

library(ggplot2)

set.seed(0)
ind <- sample(1:n, 10)

dat <- data.frame(c(1:10 - 0.13, 1:10 + 0.13), Method = factor(c(rep("Target", 10), rep("dShrink", 10))),
                  c(EST_tar[ind], EST_Prop[ind]), rbind(CI_tar[ind, ], CI_Prop[ind, ]))
colnames(dat) <- c("index", "Method", "EST", "LB", "UB")

ggplot(dat) +
  geom_pointrange(aes(x=index, y = EST, ymin = LB, ymax = UB,  color=Method), shape = 4, size = 1, linewidth = 1) +
  theme_light() +
  labs(y = "CI", x = "Index") +
  scale_x_continuous(breaks = 1:10) + 
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey", size = 1) +
theme(axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      legend.title = element_blank(),                   
      legend.text = element_text(size = 20, face = "bold"),
      axis.title=element_text(size=20),
      legend.key.size=unit(1.5,'cm'))

ggsave(paste("CI_ITE", ".pdf", sep = ""),
       path = ".", device = "pdf", width = 10, height = 6, units = "in")

