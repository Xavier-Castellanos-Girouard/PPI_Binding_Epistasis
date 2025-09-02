# Xavier Castellanos-Girouard

# Date First Created: March 14 2024 
# Date Last Modified: May 27 2024


#### Import Libraries ####

library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggsignif)

#### Import Datasets ####

## Get KaKs (dN/dS) measurements for all paralog pairs
all_paralogs_kaks_DF <- 
  read.csv("../results/Paralog_Analyses/all_paralog_kaks.csv",
           row.names = 1)

## Kd/GI dataset with Divergence measurements
Kd_GI_AllPara_DF <- 
  read.csv("../results/Paralog_Analyses/Kd_GI_Divergence.csv",
           row.names = 1)

## Stoichiometry/GI dataset with Paralog Type
Stoich_GI_AllPara_DF <- 
  read.csv("../results/Paralog_Analyses/Kd_GI_ParalogType.csv")

#### GI Graphs ####


## Make labels that are appropriate for legend
Kd_GI_AllPara_DF$Label <- NA

Kd_GI_AllPara_DF$Label[Kd_GI_AllPara_DF$Type == "NoParalog"] <- "No Paralog"
Kd_GI_AllPara_DF$Label[Kd_GI_AllPara_DF$Type == "OneParalog"] <- "One Paralog"
Kd_GI_AllPara_DF$Label[Kd_GI_AllPara_DF$Type == "TwoParalog"] <- "Two Paralogs"
Kd_GI_AllPara_DF$Label[Kd_GI_AllPara_DF$Type == "IntParalog"] <- "Paralogous\nHeteromer"

Kd_GI_AllPara_DF$Label <-
  factor(Kd_GI_AllPara_DF$Label,
         levels = c("No Paralog", "One Paralog", "Two Paralogs", "Paralogous\nHeteromer"))

## Check distribution in GI score and Kd 
Kd_vs_GI_para_plot <-
  ggplot(data = Kd_GI_AllPara_DF,
       mapping = aes(x = Kd,
                     y = scores,
                     color = Label)) +
  geom_point() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = 12), 
    axis.text.y = element_text(colour="black", size = 12),
    axis.ticks = element_line(color="black"),
    axis.line.x.bottom=element_line(color="black", linewidth = 0.5),
    axis.line.y.left=element_line(color="black", linewidth = 0.5),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.text=element_text(size=14),
    legend.title = element_text(size = 14)) +
    #legend.position = "none") +
  scale_x_log10() +
  scale_y_continuous(breaks = c(-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5)) +
  scale_color_manual(values = c("black", 
                                "royalblue2", 
                                "green2",
                                "violet")) +
  ylab("GI Score") +
  xlab("Kd (M)") +
  labs(color="Paralog Type")
#Kd_vs_GI_para_plot

## Positive GI
Kd_PosGI_AllPara_DF <- Kd_GI_AllPara_DF[Kd_GI_AllPara_DF$scores>0,]
PosGI_dist_para_plot <- 
  ggplot(data = Kd_PosGI_AllPara_DF,
       mapping = aes(x = Label,
                     y = scores,
                     fill = Label)) +
  geom_violin(trim = TRUE,
              color = "black",
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_signif(comparisons = list(c(1,2), c(1,3)), # Only keep significant
              #position="dodge",
              tip_length = 0.02,
              y_position = c(0, 0.25),
              #step_increase = 0.05,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = 12), 
    axis.text.y = element_text(colour="black", size = 12),
    axis.ticks = element_line(color="black"),
    axis.line.x.bottom=element_line(color="black", linewidth = 0.5),
    axis.line.y.left=element_line(color="black", linewidth = 0.5),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.position = "none") +
  scale_fill_manual(values = c("grey50", 
                               "royalblue2", 
                               "green2",
                               "violet")) +
  ylab("Positive GI Score") +
  xlab("Paralog Type") +
  scale_x_discrete(labels = c("None", "One", "Two", "Paralogous\nHeteromer")) +
  scale_y_log10(limits = c(10^-4, 10^1),
                breaks = c(10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1),
                labels = c(expression(10^"-4"), expression(10^"-3"), 
                           expression(10^"-2"), expression(10^"-1"),
                           expression(10^0), expression(10^1)))
#PosGI_dist_para_plot


## Negative GI
Kd_NegGI_AllPara_DF <- Kd_GI_AllPara_DF[Kd_GI_AllPara_DF$scores<0,]
NegGI_dist_para_plot <- 
  ggplot(data = Kd_NegGI_AllPara_DF,
       mapping = aes(x = Label,
                     y = abs(scores),
                     fill = Label)) +
  geom_violin(trim = TRUE,
              color = "black",
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_signif(comparisons = list(c(1,2), c(2,3), c(3,4), c(1,3), c(2,4)), # Only keep significant
              #position="dodge",
              tip_length = 0.02,
              y_position = c(0, 0, 0, 0.25, 0.5),
              #step_increase = 0.05,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = 12), 
    axis.text.y = element_text(colour="black", size = 12),
    axis.ticks = element_line(color="black"),
    axis.line.x.bottom=element_line(color="black", linewidth = 0.5),
    axis.line.y.left=element_line(color="black", linewidth = 0.5),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.position = "none") +
  scale_fill_manual(values = c("grey50", 
                               "royalblue2", 
                               "green2",
                               "violet")) +
  ylab("Negative GI Score") +
  xlab("Paralog Type") +
  scale_x_discrete(labels = c("None", "One", "Two", "Paralogous\nHeteromer")) +
  scale_y_log10(limits = c(10^-4, 10^1),
                breaks = c(10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1),
                labels = c(expression("-"*10^"-4"), expression("-"*10^"-3"), 
                           expression("-"*10^"-2"), expression("-"*10^"-1"),
                           expression("-"*10^0), expression("-"*10^1)))
#NegGI_dist_para_plot

#### Stoichiometry Graphs ####

Stoich_GI_AllPara_DF$Type <-
  factor(Stoich_GI_AllPara_DF$Type,
         levels = c("NoParalog", "OneParalog", "TwoParalog", "IntParalog"))

# Stoichiometry map
Stoich_GI_Para_p <- 
  ggplot(data = Stoich_GI_AllPara_DF,
       mapping = aes(x = Interaction_Stoichiometry,
                     y = Abundance_Stoichiometry,
                     color = Type)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  theme_classic() +
  scale_color_manual(values = c("black", 
                                "royalblue2", 
                                "green2", 
                                "violet"))

#tiff("../graphs/Paralogs_vs_StoichMap.tiff", units = "mm", height = 120, width = 180, res = 300)
#Stoich_GI_Para_p
#dev.off()


# Interaction Stoichiometry 
IntStoich_Para_p <-
  ggplot(data = Stoich_GI_AllPara_DF,
       mapping = aes(x = Type,
                     y = Interaction_Stoichiometry,
                     fill = Type)) +
  geom_violin(trim = TRUE,
              draw_quantiles = c(0.25, 0.5, 0.75),
              color = "black") +
  geom_signif(comparisons = list(c(1,2), c(1,3), c(2,4), c(3,4)), # Only keep significant
              #position="dodge",
              tip_length = 0.01,
              step_increase = 0.05,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black"), 
    axis.text.y = element_text(colour="black"),
    axis.ticks = element_line(color="black"),
    axis.line.x.bottom=element_line(color="black"),
    axis.line.y.left=element_line(color="black"),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_y_log10(limits = c(10^-6, 10^7),
    breaks = c(10^-6, 10^-4, 10^-2, 10^0, 10^2, 10^4, 10^6),
    labels = c(expression(10^-6), expression(10^-4), expression(10^-2),
               expression(10^0), expression(10^2), expression(10^4), 
               expression(10^6))) +
  scale_x_discrete(labels = c("No Paralog", "One Paralog", "Two Paralogs", "Paralogous\nHeteromer")) +
  scale_fill_manual(values = c("grey50", 
                                "royalblue2", 
                                "green2", 
                                "violet")) +
  xlab("Paralog Type") +
  ylab("Interaction Stoichiometry")

#IntStoich_Para_p



# Abundance Stoichiometry
AbStoich_Para_p <-
  ggplot(data = Stoich_GI_AllPara_DF,
       mapping = aes(x = Type,
                     y = Abundance_Stoichiometry,
                     fill = Type)) +
  geom_violin(trim = TRUE,
              draw_quantiles = c(0.25, 0.5, 0.75),
              color = "black") +
  geom_signif(comparisons = list(c(1,2), c(2,3)), # Only keep significant
              tip_length = 0.01,
              step_increase = 0.05,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black"), 
    axis.text.y = element_text(colour="black"),
    axis.ticks = element_line(color="black"),
    axis.line.x.bottom=element_line(color="black"),
    axis.line.y.left=element_line(color="black"),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_y_log10(limits = c(10^-3, 10^3),
                breaks = c(10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3),
                labels = c(expression(10^-3), expression(10^-2), expression(10^-1),
                           expression(10^0), expression(10^1), expression(10^2),
                           expression(10^3))) +
  scale_x_discrete(labels = c("No Paralog", "One Paralog", "Two Paralogs", "Paralogous\nHeteromer")) +
  scale_fill_manual(values = c("grey50", 
                               "royalblue2", 
                               "green2", 
                               "violet")) +
  xlab("Paralog Type") +
  ylab("Abundance Stoichiometry")

#tiff("../graphs/Paralogs_AbundanceStoich.tiff", units = "mm", height = 120, width = 180, res = 300)
#AbStoich_Para_p
#dev.off()


#### Export ####

legend <- get_legend(Kd_vs_GI_para_plot)
Kd_vs_GI_para_plot <- Kd_vs_GI_para_plot + theme(legend.position = "none")

Combined_plot <-
  plot_grid(Kd_vs_GI_para_plot,
            legend,
            NegGI_dist_para_plot,
            PosGI_dist_para_plot,
            labels = c("a","","b", "c"),
            rel_widths = c(1, 1, 1, 1),
            ncol = 2)

png("../results/Figures/Extended_Data_Figures/SupplementaryFigure9.png", units="mm", width=240, height=180, res=300)
Combined_plot
dev.off()


