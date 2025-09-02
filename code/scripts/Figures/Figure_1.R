# Xavier Castellanos-Girouard
# 
#
# Date First Created: February 7 2024
# Date Last Modified: September 6 2024

#### Import libraries ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggsignif)

#### Import data ####

Yeast_Stoich_DF <- read.csv("../results/Stoichiometry_and_GI/Yeast_Interactome_Stoich.csv",
                            row.names = 1)

Yeast_Stoich_GI_DF <- read.csv("../results/Stoichiometry_and_GI/Yeast_Stoich_GI.csv",
                               row.names = 1)

Yeast_Stoich_GI_wReg_DF <- read.csv("../results/Stoichiometry_and_GI/Yeast_Stoich_GI_wRegions.csv",
                                    row.names = 1)

Human_Stoich_DF <- read.csv("../results/Stoichiometry_and_GI/Human_Interactome_Stoich.csv",
                            row.names = 1)

Human_Stoich_GI_DF <- read.csv("../results/Stoichiometry_and_GI/Human_Stoich_GI.csv",
                               row.names = 1)

Human_Stoich_GI_wReg_DF <- read.csv("../results/Stoichiometry_and_GI/Human_Stoich_GI_wRegions.csv",
                                    row.names = 1)


#### Set Text Sizes ####

Axis_text_size = 5
Axis_title_size = 6
plot_title_size = 7

Axis_tick_width = 0.2
Axis_tick_length = 0.5

geom_point_size = 0.5
geom_violin_linewidth = 0.2

xy_plotline_width = 0.2

stoichplot_dashed_linewidth = 0.2

geom_signif_linewidth = 0.3
geom_signif_starsize = 2.5

#### Yeast Standard Stoichiometry plot ####


Yeast_Stoichiometry_plot <-
  ggplot(data = Yeast_Stoich_DF,
         mapping = aes(x = Interaction_Stoichiometry,
                       y = Abundance_Stoichiometry)) +
  geom_point(alpha = 0.05,
             size = geom_point_size) +
  geom_hline(yintercept = 1, 
             linetype='dashed', 
             color = "grey60", 
             linewidth = stoichplot_dashed_linewidth) +
  geom_vline(xintercept = 1, 
             linetype='dashed', 
             color = "grey60", 
             linewidth = stoichplot_dashed_linewidth) +
  geom_density2d(bins = 7,
                 color = "red",
                 linewidth = 0.3) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(color = "black", linewidth = 0.5),
        #axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.5),
        legend.text = element_text(colour = "black", 
                                   size = Axis_text_size),
        legend.title = element_text(colour = "black",
                                    size = Axis_title_size),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white")) +
  xlab("Interaction Stoichiometry") +
  ylab("Abundance Stoichiometry") +
  scale_x_log10(breaks = c(10^-6, 10^-4, 10^-2, 10^0, 10^2, 10^4),
                labels = c(expression(10^"-6"), expression(10^"-4"), 
                           expression(10^"-2"), expression(10^0), 
                           expression(10^2), expression(10^4))) +
  scale_y_log10(breaks = c(10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-3"), expression(10^"-2"), 
                           expression(10^"-1"), expression(10^0),
                           expression(10^1), expression(10^2)))


#Yeast_Stoichiometry_plot



#### Yeast Stoichiometry StatPlot with Mean GIs ####


## Positive GI stat Plot
Yeast_PosGI_Stoichiometry_plot <-
  ggplot(data = Yeast_Stoich_GI_DF[(Yeast_Stoich_GI_DF$GI_scores > 0),],
         mapping = aes(x = Interaction_Stoichiometry, 
                       y = Abundance_Stoichiometry,
                       z = GI_scores)) +
  stat_summary_2d(bins = 30, fun = "mean") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(color ="black", linewidth = 0.5),
        #axis.line.y = element_line(color ="black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black",
                                  size = plot_title_size),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    linewidth=0.5),
        legend.text = element_text(colour = "black", 
                                   size = Axis_text_size),
        legend.title = element_text(colour = "black",
                                    size = Axis_title_size),
        #legend.text.align = 0,
        #legend.key = element_rect(fill = "white"),
        #legend.position = "none",
        legend.key = element_rect(fill = "black", 
                                  linewidth = 1),
        legend.key.size = unit(3, "mm"),
        legend.justification = c(0.5,0.75)) +
  scale_x_log10(expand = c(0, 0),
                breaks = c(10^-4, 10^-2, 10^0, 10^2, 10^4),
                labels = c(expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0), expression(10^2), 
                           expression(10^4)),
                limits = c(10^-5, 10^4)) +
  scale_y_log10(expand = c(0, 0),
                breaks = c(10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-3"), expression(10^"-2"), 
                           expression(10^"-1"), expression(10^0),
                           expression(10^1), expression(10^2)),
                limits = c(10^-3, 10^2.5)) +
  scale_fill_continuous(type = "viridis",
                        limit = c(0, .2),
                        oob = scales::squish,
                        breaks = c(0, 0.05, 0.1, 0.15, 0.2),
                        guide = guide_colorbar(frame.colour = "black",
                                               frame.linewidth = 0.25,
                                               ticks = TRUE,
                                               ticks.colour='black')) +
  xlab("Interaction Stoichiometry") +
  ylab("Abundance Stoichiometry") +
  labs(fill = '   Mean\nGI Score') +
  ggtitle("Positive GI")

#Yeast_PosGI_Stoichiometry_plot


## Negative GI stat Plot
Yeast_NegGI_Stoichiometry_plot <-
  ggplot(data = Yeast_Stoich_GI_DF[(Yeast_Stoich_GI_DF$GI_scores < 0),],
         mapping = aes(x = Interaction_Stoichiometry, 
                       y = Abundance_Stoichiometry,
                       z = abs(GI_scores))) +
  stat_summary_2d(bins = 30, fun = "mean") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(color ="black", linewidth = 0.5),
        #axis.line.y = element_line(color ="black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black",
                                  size = plot_title_size),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    linewidth=0.5),
        legend.text = element_text(colour = "black", 
                                   size = Axis_text_size),
        legend.title = element_text(colour = "black",
                                    size = Axis_title_size),
        #legend.text.align = 0,
        #legend.key = element_rect(fill = "white"),
        #legend.position = "none",
        legend.key = element_rect(fill = "black", 
                                  linewidth = 1),
        legend.key.size = unit(3, "mm"),
        legend.justification = c(0.5,0.75)) +
  scale_x_log10(expand = c(0, 0),
                breaks = c(10^-4, 10^-2, 10^0, 10^2, 10^4),
                labels = c(expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0), expression(10^2), 
                           expression(10^4)),
                limits = c(10^-5, 10^4)) +
  scale_y_log10(expand = c(0, 0),
                breaks = c(10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-3"), expression(10^"-2"), 
                           expression(10^"-1"), expression(10^0),
                           expression(10^1), expression(10^2)),
                limits = c(10^-3, 10^2.5)) +
  scale_fill_continuous(type = "viridis",
                        limit = c(0, .2),
                        oob = scales::squish,
                        breaks = c(0, 0.05, 0.1, 0.15, 0.2),
                        labels = c("|0.00|", "|0.05|", "|0.10|", "|0.15|", "|0.20|"),
                        guide = guide_colorbar(frame.colour = "black",
                                               frame.linewidth = 0.25,
                                               ticks = TRUE,
                                               ticks.colour='black')) +
  xlab("Interaction Stoichiometry") +
  ylab("Abundance Stoichiometry") +
  labs(fill = '   Mean\nGI Score') +
  ggtitle("Negative GI")

#Yeast_NegGI_Stoichiometry_plot


#### Yeast Stoichiometry Region vs GI Violin plot ####

## Set to factor to conserve order in graph
Yeast_Stoich_GI_wReg_DF$Region <- 
  factor(Yeast_Stoich_GI_wReg_DF$Region, 
         levels = 1:4)

PosGI_StoichReg_violin_plot <-
  ggplot(data = Yeast_Stoich_GI_wReg_DF[(Yeast_Stoich_GI_wReg_DF$GI_scores > 0),],
         mapping = aes(x = Region,
                       y = abs(GI_scores))) +
  geom_violin(trim = TRUE,
              draw_quantiles = c(0.25, 0.5, 0.75),
              fill = "yellow2",
              color = "black", 
              linewidth = geom_violin_linewidth) +
  geom_signif(comparisons = list(c(1,2), c(1,3)),
              tip_length = 0.01,
              step_increase = 0.07,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
    axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_y_log10(breaks = c(10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1),
                labels = c(expression(10^"-4"), expression(10^"-3"), 
                           expression(10^"-2"), expression(10^"-1"),
                           expression(10^0), expression(10^1)),
                limits = c(10^-4, 10^1)) +
  scale_x_discrete(labels = c("Reg. 1", "Reg. 2", "Reg. 3", "Reg. 4")) +
  xlab("Stoichiometry Region") +
  ylab("Positive GI Score")

#PosGI_StoichReg_violin_plot




NegGI_StoichReg_violin_plot <-
  ggplot(data = Yeast_Stoich_GI_wReg_DF[(Yeast_Stoich_GI_wReg_DF$GI_scores < 0),],
         mapping = aes(x = Region,
                       y = abs(GI_scores))) +
  geom_violin(trim = TRUE,
              draw_quantiles = c(0.25, 0.5, 0.75),
              fill = "blue",
              color = "black",
              alpha = 0.75,
              linewidth = geom_violin_linewidth) +
  geom_signif(comparisons = list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4)),
              tip_length = 0.01,
              #step_increase = 0.05,
              #y_position = c(0, 0.2, 0.6, 0, 0.4),
              y_position = c(0, 0.25, 0.75, 0, 0.5),
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
    axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_y_log10(breaks = c(10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1),
                labels = c(expression("-"*10^"-4"), expression("-"*10^"-3"), 
                           expression("-"*10^"-2"), expression("-"*10^"-1"),
                           expression("-"*10^0), expression("-"*10^1)),
                limits = c(10^-4, 10^1)) +
  scale_x_discrete(labels = c("Reg. 1", "Reg. 2", "Reg. 3", "Reg. 4")) +
  xlab("Stoichiometry Region") +
  ylab("Negative GI Score")

#NegGI_StoichReg_violin_plot


#### Human Standard Stoichiometry plot ####


## plot
Human_Stoichiometry_plot <-
  ggplot(data = Human_Stoich_DF,
         mapping = aes(x = Interaction_Stoichiometry,
                       y = Abundance_Stoichiometry)) +
  geom_point(alpha = 0.02,
             size = geom_point_size) +
  geom_hline(yintercept = 1, 
             linetype='dashed', 
             color = "grey60",
             linewidth = stoichplot_dashed_linewidth) +
  geom_vline(xintercept = 1, 
             linetype='dashed', 
             color = "grey60",
             linewidth = stoichplot_dashed_linewidth) +
  geom_density2d(bins = 7,
                 color = "red",
                 linewidth = 0.3) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    #axis.line.x = element_line(color ="black", linewidth = 0.5),
    #axis.line.y = element_line(color ="black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(colour = "black",
                                fill = NA,
                                linewidth = 0.5)) +
  xlab("Interaction Stoichiometry") +
  ylab("Abundance Stoichiometry") +
  scale_x_log10(limits = c(10^-8, 10^4),
                breaks = c(10^-8, 10^-6, 10^-4, 10^-2, 10^0, 10^2, 10^4),
                labels = c(expression(10^"-8"), expression(10^"-6"),
                           expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0), expression(10^2), 
                           expression(10^4))) +
  scale_y_log10(breaks = c(10^-4, 10^-2, 10^0, 10^2, 10^4, 10^6),
                labels = c(expression(10^"-4"), expression(10^"-2"),
                           expression(10^0), expression(10^2), 
                           expression(10^4), expression(10^6)))
#Human_Stoichiometry_plot


#### Human Stoichiometry StatPlot with Mean GIs ####

## Positive GI stat Plot
Human_PosGI_Stoichiometry_plot <-
  ggplot(data = Human_Stoich_GI_DF[(Human_Stoich_GI_DF$GI_Score > 0),],
         mapping = aes(x = Interaction_Stoichiometry, 
                       y = Abundance_Stoichiometry,
                       z = GI_Score)) +
  stat_summary_2d(bins = 18, fun = "mean") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(color ="black", linewidth = 0.5),
        #axis.line.y = element_line(color ="black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black",
                                  size = plot_title_size),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    linewidth=0.5),
        legend.text = element_text(colour = "black", 
                                   size = Axis_text_size),
        legend.title = element_text(colour = "black",
                                    size = Axis_title_size),
        #legend.text.align = 0,
        #legend.key = element_rect(fill = "white"),
        #legend.position = "none",
        legend.key = element_rect(fill = "black", linewidth = 1),
        legend.key.size = unit(3, "mm"),
        legend.justification = c(0.5,0.75)) +
  scale_x_log10(expand = c(0, 0),
                breaks = c(10^-6, 10^-4, 10^-2, 10^0, 10^2, 10^4),
                labels = c(expression(10^-"6"), expression(10^"-4"), 
                           expression(10^"-2"), expression(10^0), 
                           expression(10^2), expression(10^4)),
                limits = c(10^-6, 10^2)) +
  scale_y_log10(expand = c(0, 0),
                breaks = c(10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-3"), expression(10^"-2"), 
                           expression(10^"-1"), expression(10^0),
                           expression(10^1), expression(10^2)),
                limits = c(10^-3, 10^2.5)) +
  scale_fill_continuous(type = "viridis",
                        limit = c(0, 3),
                        oob = scales::squish,
                        breaks = c(0, 1, 2, 3),
                        guide = guide_colorbar(frame.colour = "black",
                                               frame.linewidth = 0.25,
                                               ticks = TRUE,
                                               ticks.colour='black')) +
  xlab("Interaction Stoichiometry") +
  ylab("Abundance Stoichiometry") +
  labs(fill = '   Mean\nGI Score') +
  ggtitle("Positive GI")

#Human_PosGI_Stoichiometry_plot


## Negative GI Stat Plot
Human_NegGI_Stoichiometry_plot <-
  ggplot(data = Human_Stoich_GI_DF[(Human_Stoich_GI_DF$GI_Score < 0),],
         mapping = aes(x = Interaction_Stoichiometry, 
                       y = Abundance_Stoichiometry,
                       z = abs(GI_Score))) +
  stat_summary_2d(bins = 18, fun = "mean") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.line.x = element_line(color ="black", linewidth = 0.5),
        #axis.line.y = element_line(color ="black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black",
                                  size = plot_title_size),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    linewidth=0.5),
        legend.text = element_text(colour = "black", 
                                   size = Axis_text_size),
        legend.title = element_text(colour = "black",
                                    size = Axis_title_size),
        #legend.text.align = 0,
        #legend.key = element_rect(fill = "white"),
        #legend.position = "none",
        legend.key = element_rect(fill = "black", linewidth = 0.5),
        legend.key.size = unit(3, "mm"),
        legend.justification = c(0.5,0.75)) +
  scale_x_log10(expand = c(0, 0),
                breaks = c(10^-6, 10^-4, 10^-2, 10^0, 10^2, 10^4),
                labels = c(expression(10^"-6"), expression(10^"-4"), 
                           expression(10^"-2"), expression(10^0), 
                           expression(10^2), expression(10^4)),
                limits = c(10^-6, 10^2)) +
  scale_y_log10(expand = c(0, 0),
                breaks = c(10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-3"), expression(10^"-2"), 
                           expression(10^"-1"), expression(10^0),
                           expression(10^1), expression(10^2)),
                limits = c(10^-3, 10^2.5)) +
  scale_fill_continuous(type = "viridis",
                        limit = c(0, 3),
                        oob = scales::squish,
                        breaks = c(0, 1, 2, 3),
                        labels = c("|0.0|", "|1.0|", "|2.0|", "|3.0|"),
                        guide = guide_colorbar(frame.colour = "black",
                                               frame.linewidth = 0.25,
                                               ticks = TRUE,
                                               ticks.colour='black')) +
  xlab("Interaction Stoichiometry") +
  ylab("Abundance Stoichiometry") +
  labs(fill = '   Mean\nGI Score') +
  ggtitle("Negative GI")

#Human_NegGI_Stoichiometry_plot



#### Human Stoichiometry Region vs GI Violin pot ####

## Set to factor to conserve order in graph
Human_Stoich_GI_wReg_DF$Region <- 
  factor(Human_Stoich_GI_wReg_DF$Region, 
         levels = 1:4)

Human_PosGI_StoichReg_violin_plot <-
  ggplot(data = Human_Stoich_GI_wReg_DF[(Human_Stoich_GI_wReg_DF$GI_Score > 0),],
         mapping = aes(x = Region,
                       y = abs(GI_Score))) +
  geom_violin(trim = TRUE,
              draw_quantiles = c(0.25, 0.5, 0.75),
              fill = "yellow2",
              color = "black",
              linewidth = geom_violin_linewidth) +
  geom_signif(comparisons = list(c(1,2), c(1,3)),
              tip_length = 0.01,
              step_increase = 0.09,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
    axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(3, "mm"),
    legend.justification = c(0.5,0.75)) +
  scale_y_log10(breaks = c(10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-2"), expression(10^"-1"), 
                           expression(10^0), expression(10^1),
                           expression(10^2)),
                limits = c(10^-2.2, 10^1.5)) +
  scale_x_discrete(labels = c("Reg. 1", "Reg. 2", "Reg. 3", "Reg. 4")) +
  xlab("Stoichiometry Region") +
  ylab("Positive GI Score")

#Human_PosGI_StoichReg_violin_plot



Human_NegGI_StoichReg_violin_plot <-
  ggplot(data = Human_Stoich_GI_wReg_DF[(Human_Stoich_GI_wReg_DF$GI_Score < 0),],
         mapping = aes(x = Region,
                       y = abs(GI_Score))) +
  geom_violin(trim = TRUE,
              draw_quantiles = c(0.25, 0.5, 0.75),
              fill = "blue",
              color = "black",
              alpha = 0.75,
              linewidth = geom_violin_linewidth) +
  geom_signif(comparisons = list(c(1,2), c(1,3)),
              tip_length = 0.01,
              step_increase = 0.06,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
    axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    legend.key.size = unit(3, "mm"),
    legend.justification = c(0.5,0.75)) +
  scale_y_log10(breaks = c(10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression("-"*10^"-2"), expression("-"*10^"-1"), 
                           expression("-"*10^0), expression("-"*10^1),
                           expression("-"*10^2)),
                limits = c(10^-2.2, 10^1.5)) +
  scale_x_discrete(labels = c("Reg. 1", "Reg. 2", "Reg. 3", "Reg. 4")) +
  xlab("Stoichiometry Region") +
  ylab("Negative GI Score")

Human_NegGI_StoichReg_violin_plot



#### Combined plot Display ####


## Stoichiometry Plots

TopDivisions <-
  plot_grid(NULL,
            NULL,
            NULL,
            NULL,
            labels = c("a", "b\tYeast", "c\tHuman", "d"),
            label_size = 8,
            label_fontface = "bold",
            rel_widths = c(1, 1, 1, 1),
            ncol = 4)

TopPlots <- 
  plot_grid(
    NULL,
    Yeast_Stoichiometry_plot,
    Human_Stoichiometry_plot,
    NULL,
    #labels = c("A", "    B", "    C"),
    labels = c(),
    label_size = 8,
    label_fontface = "bold",
    rel_widths = c(1, 1, 1, 1),
    ncol = 4)

#TopPlots

## GI Stoichiometry Plots

# Get legends
yeastGILegend <- get_legend(Yeast_NegGI_Stoichiometry_plot)
humanGILegend <- get_legend(Human_NegGI_Stoichiometry_plot)

# Remove legends from yeast plots
Yeast_NegGI_Stoichiometry_plot_noLeg <- Yeast_NegGI_Stoichiometry_plot + theme(legend.position = "none")
Yeast_PosGI_Stoichiometry_plot_noLeg <- Yeast_PosGI_Stoichiometry_plot + theme(legend.position = "none")

# Remove legends from human plots
Human_NegGI_Stoichiometry_plot_noLeg <- Human_NegGI_Stoichiometry_plot  + theme(legend.position = "none")
Human_PosGI_Stoichiometry_plot_noLeg <- Human_PosGI_Stoichiometry_plot + theme(legend.position = "none")

MidPlots <-
  plot_grid(
    plot_grid(Yeast_NegGI_Stoichiometry_plot_noLeg,
              Yeast_PosGI_Stoichiometry_plot_noLeg,
              yeastGILegend,
              NULL, 
              labels = c("e", "", "", ""),
              label_size = 8,
              label_fontface = "bold",
              rel_widths = c(1, 1, 0.25, 0.05),
              ncol = 4),
    plot_grid(Human_NegGI_Stoichiometry_plot_noLeg,
              Human_PosGI_Stoichiometry_plot_noLeg,
              humanGILegend,
              NULL,
              labels = c("f", "", "", ""),
              label_size = 8,
              label_fontface = "bold",
              rel_widths = c(1, 1, 0.25, 0.05),
              ncol = 4),
    labels = c("Yeast", "Human"),
    label_size = 8,
    label_fontface = "bold",
    label_x = 0.4,
    label_y = 1.1,
    rel_widths = c(1, 1),
    ncol = 2
  )


## Stoichiometry Region GI violin plots
BotPlots <- plot_grid(
  NegGI_StoichReg_violin_plot,
  PosGI_StoichReg_violin_plot,
  NULL,
  Human_NegGI_StoichReg_violin_plot,
  Human_PosGI_StoichReg_violin_plot,
  NULL,
  labels = c("g", "", "", "h", "", ""),
  label_size = 8,
  label_fontface = "bold",
  rel_widths = c(1, 1, 0.25, 1, 1, 0.25),
  ncol = 6)

# Yeast and Human Combined
CombinedPlot <-
  plot_grid(TopDivisions,
            TopPlots,
            NULL,
            MidPlots,
            BotPlots,
            labels = c("", "", "", ""),
            rel_heights = c(0.1, 1, 0.1, 1, 1),
            ncol = 1)

png("../results/Figures/Main_Figures/Figure_1.png", units="mm", width=180, height=120, res=300)
CombinedPlot
dev.off()

mm_to_inch = 0.0393701
svg("../results/Figures/Main_Figures/Figure_1.svg", width=180*mm_to_inch, height=120*mm_to_inch)
CombinedPlot
dev.off()
