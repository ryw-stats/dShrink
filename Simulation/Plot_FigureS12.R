setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(latex2exp)
library(ggplot2)

p <- 25
n <- 100
m <- 200


load(paste(paste("Error", p, n, m, "Cmpr", "sameVariable", "CI", sep = "-"), ".Rdata", sep = ""))



CI <- data.frame(rep(A_seq, 4),
                 c(t(apply(Res_seq[5:8, , ], c(1, 3), mean))),
                 rep(c("ds",  "T", "TransGLM", "ISEDI"), rep(length(A_seq), 4)))

colnames(CI) <- c("t", "cover", "label")

CI$label <- factor(CI$label, levels= c("T", "TransGLM", "ISEDI", "ds"))


color <- c("T" = "green4", "TransGLM" = "red", "ISEDI" = "darkgoldenrod1", "ds" = "purple")

p1 <- ggplot(CI) + geom_line(aes(t, cover, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, cover, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "TransGLM", "ISEDI", "dShrink")) +
  coord_cartesian(ylim=c(0.5, 1)) +
  scale_shape_manual(values = c(1, 15, 16, 17), labels =c("Target Pop", "TransGLM", "ISEDI", "dShrink")) +
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

CI <- data.frame(rep(A_seq, 4),
                 c(t(apply(Res_seq[1:4, , ], c(1, 3), mean))),
                 rep(c("ds",  "T", "TransGLM", "ISEDI"), rep(length(A_seq), 4)))

colnames(CI) <- c("t", "width", "label")

CI$label <- factor(CI$label, levels= c("T", "TransGLM", "ISEDI", "ds"))


p2 <- ggplot(CI) + geom_line(aes(t, width, color = label), size = 1.2, linetype = 2) +
  geom_point(aes(t, width, color = label, shape = label), size = 4) +
  scale_color_manual(values = color, labels = c("Target Pop", "TransGLM", "ISEDI", "dShrink")) +  
  coord_cartesian(ylim=c(0, 0.7)) +
  scale_shape_manual(values = c(1, 15, 16, 17), labels =c("Target Pop", "TransGLM", "ISEDI", "dShrink")) +
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

ggsave(paste(paste("CI", p, n, m, "ds", "comp", "sameVariable", sep = "-"), ".pdf", sep = ""),
       path = "fig", device = "pdf", width = 12, height = 6, units = "in")
