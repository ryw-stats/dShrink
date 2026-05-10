setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(latex2exp)
library(ggplot2)

p <- 30
n <- 300
m <- 1000

load(paste(paste("Error", p, n, m, "PCML", "hVar", sep = "-"), ".Rdata", sep = ""))
Res <- colMeans(Res_seq)
load(paste(paste("Error", p, n, m, "Cmpr", "hVar", sep = "-"), ".Rdata", sep = ""))
Res_seq[Res_seq > 100] <- NA
Res <- rbind(apply(Res_seq, c(1, 3), function(x) mean(x, na.rm = T)), Res)
##ylim is c(0.12, 0.135) when p = 12; c(0.065, 0.14) when p = 30

A_seq <- seq(0, 1, 0.1)

Error <- data.frame(rep(A_seq, 6),
                    as.vector(t(Res[c(2, 4, 3, 6, 5, 1), ])),
                    rep(c("T", "CML", "GIM", "PCML", "JS", "ds"), rep(length(A_seq), 6)))

colnames(Error) <- c("t", "Error", "label")

Error$label <- factor(Error$label, levels= c("T", "CML", "GIM", "PCML", "JS", "ds"))


color <- c("T" = "green4", "CML" = "blue", "GIM" = "red", "PCML" = "darkgoldenrod1", "JS" = "deepskyblue3", "ds" = "purple")

ggplot(Error) + geom_line(aes(t, Error, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, Error, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "CML", "GIM", "PCML", "JS Comb", "dShrink")) +
  coord_cartesian(ylim=c(0.05, 0.15)) +
  scale_shape_manual(values = c(1, 0, 15, 16, 18, 17), labels = c("Target Pop", "CML", "GIM", "PCML", "JS Comb", "dShrink")) +
  labs(y = "MSE", x = "t") +
  theme_light() +
  theme(axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25, face = "bold"),
        legend.position = "top",
        axis.title=element_text(size=30),
        legend.key.size=unit(1.5,'cm'))

ggsave(paste(paste("Error", p, n, m, "Cmpr", "hVar", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 8, units = "in")
