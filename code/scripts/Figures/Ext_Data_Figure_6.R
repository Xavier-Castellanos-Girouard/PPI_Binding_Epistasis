# Xavier Castellanos-Girouard
# 
# Determining the Relationship between free energies (from Kd) and Epistasis
#
# Date First Created: November 9 2023
# Date Last Modified: December 2 2024

#### Import libraries ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(viridis)

#### Import Data ####

Yeast_Kd_GI_DF <- read.csv("../results/Kd_and_GI/Yeast_Kd_GI.csv", row.names = 1)
Human_Kd_GI_DF <- read.csv("../results/Kd_and_GI/Human_Kd_GI.csv", row.names = 1)

#### Set Text Sizes ####

Axis_text_size = 12
Axis_title_size = 14
plot_title_size = 16

Axis_tick_width = 0.5
Axis_tick_length = 1

xy_plotline_width = 0.5

geom_violin_linewidth = 0.5


#### Yeast 2d Density plots ####


## Binding
Yeast_Binding_density_2d <-
  ggplot(data = Yeast_Kd_GI_DF,
         mapping = aes(x = Kd, 
                       y = scores)) +
  stat_density_2d_filled(bins = 30,
                         contour_var = "count") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = xy_plotline_width),
    axis.line.y = element_line(color = "black", linewidth = xy_plotline_width),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    plot.title = element_text(size = plot_title_size,
                              hjust = 0.5), 
    legend.position = "none") +
  scale_x_log10(expand = c(0,0),
                limits = c(10^-9, 10^-1.3),
                breaks = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1),
                labels = c(expression(10^"-9"), expression(10^"-8"),
                           expression(10^"-7"), expression(10^"-6"),
                           expression(10^"-5"), expression(10^"-4"),
                           expression(10^"-3"), expression(10^"-2"),
                           expression(10^"-1"))) +
  scale_y_continuous(limits = c(-0.5, 0.2),
                     expand = c(0,0),
                     breaks = c(0.4, 0.3, 0.2, 0.1, 0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1),
                     labels = c(0.4, 0.3, 0.2, 0.1, 0, -0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1)) +
  #scale_y_log10() +
  ylab("GI Score") +
  xlab("Kd (M)") +
  scale_fill_viridis_d(option = "H",
                       begin = 0,
                       end = 0.6,
                       direction = 1) +
  ggtitle("Yeast")

#Yeast_Binding_density_2d

#### Yeast GI Kd Violinplot ####

### Bin Data

## Negative GI
Kd_NegGI_Network_DF <- 
  as.data.frame(Yeast_Kd_GI_DF) %>%
  dplyr::filter(scores < 0)

# Divide into deciles according to Interaction Stoichiometry
Kd_NegGI_Network_DF$Division <- 
  dplyr::ntile(x = Kd_NegGI_Network_DF$Kd, 
               n = 10)

Kd_NegGI_Network_DF$Division <- 
  factor(Kd_NegGI_Network_DF$Division,
         levels = 1:10)

## Positive GI
# Divide into deciles according to Kd
Kd_PosGI_Network_DF <- 
  as.data.frame(Yeast_Kd_GI_DF) %>%
  dplyr::filter(scores > 0)

# Divide into deciles according to Kd
Kd_PosGI_Network_DF$Division <- 
  dplyr::ntile(x = Kd_PosGI_Network_DF$Kd, 
               n = 10)

Kd_PosGI_Network_DF$Division <- 
  factor(Kd_PosGI_Network_DF$Division,
         levels = 1:10)

## Negative GI
Yeast_negGI_Kd_lin_plot <-
  ggplot(data = Kd_NegGI_Network_DF,
         mapping = aes(x = Division,
                       y = abs(scores))) +
  geom_violin(trim = TRUE,
              fill = "blue",
              alpha = 0.75,
             draw_quantiles = c(.25, .5, .75),
             linewidth = geom_violin_linewidth) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = xy_plotline_width),
    axis.line.y = element_line(color = "black", linewidth = xy_plotline_width),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position="none",
    legend.background = element_rect(fill='transparent'),
    legend.text = element_text(size = Axis_text_size),
    legend.title = element_blank(), #element_text(colour = "black", size = Axis_title_size),
    plot.title = element_text(size = plot_title_size,
                              hjust = 0.5)) +
  scale_y_log10(breaks = c(10^-4, 10^-3, 10^-2, 10^-1, 10^0),
                labels = c(expression("-"*10^"-4"), expression("-"*10^"-3"), 
                           expression("-"*10^"-2"), expression("-"*10^"-1"),
                           expression("-"*10^0)),
                limits = c(10^-4, 10^0)) +
  xlab("Quantile") +
  ylab("Negative GI Score") +
  ggtitle("Yeast")

#Yeast_negGI_Kd_lin_plot

## Positive GI
Yeast_posGI_Kd_lin_plot <-
  ggplot(data = Kd_PosGI_Network_DF,
         mapping = aes(x = Division,
                       y = scores)) +
  geom_violin(trim = TRUE,
             fill = "yellow2",
             draw_quantiles = c(.25, .5, .75),
             linewidth = geom_violin_linewidth) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = xy_plotline_width),
    axis.line.y = element_line(color = "black", linewidth = xy_plotline_width),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position="none",
    legend.background = element_rect(fill='transparent'),
    legend.text = element_text(size = Axis_text_size),
    legend.title = element_blank(), #element_text(colour = "black", size = Axis_title_size),
    plot.title = element_text(size = plot_title_size,
                              hjust = 0.5)) +
  scale_y_log10(breaks = c(10^-4, 10^-3, 10^-2, 10^-1, 10^0),
                labels = c(expression(10^"-4"), expression(10^"-3"), 
                           expression(10^"-2"), expression(10^"-1"),
                           expression(10^0)),
                limits = c(10^-4, 10^0)) +
  xlab("Quantile") +
  ylab("Positive GI Score") +
  ggtitle("Yeast")

#Yeast_posGI_Kd_lin_plot


#### Human 2d Density plots ####


## Binding
Human_Binding_density_2d <-
  ggplot(data = Human_Kd_GI_DF,
         mapping = aes(x = Kd, 
                       y = GI_Score)) +
  geom_density_2d_filled(bins = 20) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = xy_plotline_width),
    axis.line.y = element_line(color = "black", linewidth = xy_plotline_width),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    plot.title = element_text(size = plot_title_size,
                              hjust = 0.5),
    legend.position = "none") +
  scale_x_log10(expand = c(0,0),
                limits = c(10^-9, 10^-1.3),
                breaks = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1),
                labels = c(expression(10^"-9"), expression(10^"-8"),
                           expression(10^"-7"), expression(10^"-6"),
                           expression(10^"-5"), expression(10^"-4"),
                           expression(10^"-3"), expression(10^"-2"),
                           expression(10^"-1"))) +
  scale_y_continuous(expand = c(0,0),
                     breaks = c(-10, -7.5, -5, -2.5, 0, 2.5, 5)) +
  ylab("GI Score") +
  xlab("Kd (M)") +
  scale_fill_viridis_d(option = "H",
                       begin = 0,
                       end = 0.6,
                       direction = 1) +
  ggtitle("Human")

#Human_Binding_density_2d


#### Human GI Kd Violinplot ####

## Binning

# Negative GI
Human_Kd_NegGI_Network_DF <- Human_Kd_GI_DF[Human_Kd_GI_DF$GI_Score < 0,]

Human_Kd_NegGI_Network_DF$Quintile <- 
  Human_Kd_NegGI_Network_DF$Kd%>%
  ntile(n = 5)

Human_Kd_NegGI_Network_DF$Quintile <-
  factor(Human_Kd_NegGI_Network_DF$Quintile, levels = 1:5)

# Positive GI
Human_Kd_PosGI_Network_DF <- Human_Kd_GI_DF[Human_Kd_GI_DF$GI_Score > 0,]

Human_Kd_PosGI_Network_DF$Quintile <- 
  Human_Kd_PosGI_Network_DF$Kd %>%
  ntile(n = 5)

Human_Kd_PosGI_Network_DF$Quintile <-
  factor(Human_Kd_PosGI_Network_DF$Quintile,
         levels = 1:5)

## Negative GI
Human_negGI_Kd_lin_plot <- 
  ggplot(data = Human_Kd_NegGI_Network_DF,
         mapping = aes(x = Quintile,
                       y = abs(GI_Score))) +
  geom_violin(trim = TRUE,
              size = geom_violin_linewidth,
             fill = "blue",
             alpha = 0.75,
             draw_quantiles = c(.25, .5, .75)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = xy_plotline_width),
    axis.line.y = element_line(color = "black", linewidth = xy_plotline_width),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position="none",
    #legend.text = element_text(size = Axis_text_size),
    #legend.title = element_blank(), #element_text(colour = "black", size = Axis_title_size),
    plot.title = element_text(size = plot_title_size,
                              hjust = 0.5)) +
  scale_y_log10(breaks = c(10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression("-"*10^"-2"), expression("-"*10^"-1"), 
                           expression("-"*10^0), expression("-"*10^1),
                           expression("-"*10^2)),
                limits = c(10^-2.2, 10^1.5)) +
  xlab("Quantile") +
  ylab("Negative GI Score") +
  ggtitle("Human")

#Human_negGI_Kd_lin_plot


## Positive GI
Human_posGI_Kd_lin_plot <- 
  ggplot(data = Human_Kd_PosGI_Network_DF,
         mapping = aes(x = Quintile,
                       y = GI_Score)) +
  geom_violin(trim = TRUE,
             size = geom_violin_linewidth,
             fill = "yellow2",
             draw_quantiles = c(.25, .5, .75)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = xy_plotline_width),
    axis.line.y = element_line(color = "black", linewidth = xy_plotline_width),
    axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
    axis.text.y.left = element_text(color = "black", size = Axis_text_size),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position="none",
    #legend.text = element_text(size = Axis_text_size),
    #legend.title = element_blank(), #element_text(colour = "black", size = Axis_title_size),
    plot.title = element_text(size = plot_title_size,
                              hjust = 0.5)) +
  scale_y_log10(breaks = c(10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-2"), expression(10^"-1"), 
                           expression(10^0), expression(10^1),
                           expression(10^2)),
                limits = c(10^-2.2, 10^1.5)) +
  xlab("Quantile") +
  ylab("Positive GI Score") +
  ggtitle("Human")

#Human_posGI_Kd_lin_plot

#### Combine plots ####

TopPlot <-
  plot_grid(Yeast_Binding_density_2d,
            Yeast_negGI_Kd_lin_plot,
            Yeast_posGI_Kd_lin_plot,
            ncol = 3,
            labels = c("a", "b", "c"),
            label_size = 16,
            label_fontface = "bold")

BotPlot <-
  plot_grid(Human_Binding_density_2d,
            Human_negGI_Kd_lin_plot,
            Human_posGI_Kd_lin_plot,
            ncol = 3,
            labels = c("", "", ""),
            label_size = 16,
            label_fontface = "bold")

CombinedPlot <-
  plot_grid(TopPlot,
            BotPlot,
            ncol = 1)


png("../results/Figures/Extended_Data_Figures/SupplementaryFigure6.png", units="mm", width=360, height=180, res=300)
CombinedPlot
dev.off()

#mm_to_inch = 0.0393701
#svg(paste0(dir, "Figures/Supplementary_Figure_4/graphs/SupplementaryFigure2.svg"), width=89*mm_to_inch, height=120*mm_to_inch)
#CombinedPlot
#dev.off()
