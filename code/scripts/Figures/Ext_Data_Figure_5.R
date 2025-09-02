# Xavier Castellanos-Girouard

# Date First Created: Sometime Summer 2022
# Date Last Modified: May 17 2024


#### Import Libraries ####

library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(cowplot)

#### Import data ####

Human_Est_vs_Lit_DF <- read.csv("../results/Kd_and_GI/Human_Kd_Est_vs_Lit.csv",
                                row.names = 1)

#Human_Est_vs_Lit_vs_Comp_DF <- read.csv("../results/Kd_and_GI/Human_Kd_Est_vs_Lit_vs_Comp.csv",
#                                        row.names = 1)

Yeast_Kd_Est_vs_PCA_Kd <- read.csv("../results/Kd_and_GI/Yeast_Kd_Est_vs_PCA_Kd.csv",
                                row.names = 1)


#### Make Figures ####

Estimated_Kd_vs_Literature_Kd_p <- 
  ggplot(data = Human_Est_vs_Lit_DF,
         mapping = aes(x = Kd,
                       y = Value)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "red") +
  ggpubr::stat_cor(method="pearson",
                   show.legend = FALSE,
                   digits = 2) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = 12), 
    axis.text.y = element_text(colour="black", size = 12),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.position="none") +
  scale_x_log10(breaks = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2),
                labels = c(expression(10^-9), expression(10^-8), expression(10^-7), 
                           expression(10^-6), expression(10^-5), expression(10^-4), 
                           expression(10^-3), expression(10^-2))) +
  scale_y_log10(breaks = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^-9), expression(10^-8), expression(10^-7), 
                           expression(10^-6), expression(10^-5), expression(10^-4), 
                           expression(10^-3))) +
  xlab("Estimated Kd (M)") +
  ylab("Literature Kd (M)")

#png("../graphs/Estimated_vs_Lit.png", units="mm", width=120, height=80, res=300)
#Estimated_Kd_vs_Literature_Kd_p
#dev.off()


Estimated_IntStoich_vs_Literature_Kd_p <- 
  ggplot(data = Human_Est_vs_Lit_DF,
         mapping = aes(x = Interaction_Stoichiometry,
                       y = Value)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "red") +
  ggpubr::stat_cor(method="pearson",
                   show.legend = FALSE,
                   digits = 2) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = 12), 
    axis.text.y = element_text(colour="black", size = 12),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.position="none") +
  scale_x_log10(breaks = c(10^-4, 10^-3, 10^-2, 10^-1, 10^0),
                labels = c(expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0))) +
  scale_y_log10(breaks = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^-9), expression(10^-8), expression(10^-7), 
                           expression(10^-6), expression(10^-5), 
                           expression(10^-4), expression(10^-3))) +
  xlab("Interaction Stoichiometry") +
  ylab("Literature Kd (M)")

#png("../graphs/Estimated_vs_Lit.png", units="mm", width=120, height=80, res=300)
#Estimated_IntStoich_vs_Literature_Kd_p


#### Benchmarking with Computed Kd values ####

#Estimated_computed_Kd_vs_Literature_Kd_p <- 
#  ggplot(data = Human_Est_vs_Lit_vs_Comp_DF,
#         mapping = aes(x = Computed_Kd,
#                       y = Value)) +
#  geom_point() +
#  geom_smooth(method = "lm",
#              se = FALSE,
#              color = "red") +
#  ggpubr::stat_cor(method="pearson",
#                   show.legend = FALSE,
#                   digits = 2) +
#  theme(
#    panel.grid.major = element_blank(), 
#    panel.grid.minor = element_blank(),
#    panel.background = element_blank(),
#    axis.text.x = element_text(colour="black", size = 12), 
#    axis.text.y = element_text(colour="black", size = 12),
#    axis.line.x = element_line(color = "black", linewidth = 0.5),
#    axis.line.y = element_line(color = "black", linewidth = 0.5),
#    axis.ticks = element_line(color = "black"),
#    axis.title.x.bottom = element_text(color = "black", size = 14),
#    axis.title.y.left = element_text(color = "black", size = 14),
#    legend.position="none") +
#  scale_x_log10(breaks = c(10^-25, 10^-20, 10^-15, 10^-10, 10^-5, 10^-1),
#                labels = c(expression(10^-25), expression(10^-20), 
#                           expression(10^-15), expression(10^-10), 
#                           expression(10^-5), expression(10^-1))) +
#  scale_y_log10(breaks = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
#                labels = c(expression(10^-9), expression(10^-8), expression(10^-7), 
#                           expression(10^-6), expression(10^-5), 
#                           expression(10^-4), expression(10^-3))) +
#  xlab("Computed Kd (M)") +
#  ylab("Literature Kd (M)")

#Estimated_computed_Kd_vs_Literature_Kd_p

#### Corroborate Kd values with that of DHFR PCA (Yeast) ####


Estimated_kd_vs_PCA_Kd_p <-
  ggplot(data = Yeast_Kd_Est_vs_PCA_Kd,
       mapping = aes(x = Kd,
                     y = PCA_Kd)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "red") +
  ggpubr::stat_cor(method="pearson",
                   show.legend = FALSE,
                   digits = 2,
                   label.x.npc = 0.6,
                   label.y.npc = 0.15) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = 12), 
    axis.text.y = element_text(colour="black", size = 12),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_x_log10(breaks = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2),
                labels = c(expression(10^-8), expression(10^-7), 
                           expression(10^-6), expression(10^-5), expression(10^-4), 
                           expression(10^-3), expression(10^-2))) +
  scale_y_log10(breaks = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4),
                labels = c(expression(10^-9), expression(10^-8), 
                           expression(10^-7), expression(10^-6), 
                           expression(10^-5), expression(10^-4))) +
  xlab("Mass Spectrometry Estimated Kd (M)") +
  ylab("PCA Estimated Kd (M)")

#Estimated_kd_vs_PCA_Kd_p

#### Combine Plots ####

combined_plot <- 
  plot_grid(Estimated_Kd_vs_Literature_Kd_p,
          Estimated_IntStoich_vs_Literature_Kd_p,
          Estimated_kd_vs_PCA_Kd_p,
          NULL,
          labels = c("A", "B", "C", ""),
          rel_widths = c(1, 1, 1, 1),
          ncol = 2)

png("../results/Figures/Extended_Data_Figures/SupplementaryFigure5.png", units="mm", width=240, height=160, res=300)
combined_plot
dev.off()
