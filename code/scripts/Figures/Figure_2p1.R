# Xavier Castellanos-Girouard
# 
# Determining the Relationship between free energies (from Kd) and Epistasis
#
# Date First Created: November 9 2023
# Date Last Modified: November 25 2024

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

Yeast_Kd_GI_Deciles_DF <- read.csv("../results/Kd_and_GI/Yeast_Mean_GI_KdDeciles.csv", row.names = 1)
Human_Kd_GI_Quintiles_DF <- read.csv("../results/Kd_and_GI/Human_Mean_GI_KdQuintiles.csv", row.names = 1)

#### Set Text Sizes ####

Axis_text_size = 5
Axis_title_size = 6
plot_title_size = 7

Axis_tick_width = 0.2
Axis_tick_length = 0.5

xy_plotline_width = 0.2

geomline_linewidth = 0.5
geompoint_size = 0.5


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

#### Yeast GI Kd Lineplot ####

## Negative is a sigmoidal function of general form (1/(1 + e^x))
## Specific form is (a/(1 + e^b*(x-d)) + c)
Yeast_NegGI_Model <- nls(GI_mean ~ (a/(1 + exp(b*(log10_Kd_mean-d)))+c), data = Yeast_Kd_GI_Deciles_DF[Yeast_Kd_GI_Deciles_DF$GI_Score_Type == "Negative",], start = list(b = 8, a = 0.3, c = 0.1, d = -6))
#Yeast_NegGI_Model <- nls(GI_mean ~ (1/(1 + exp(b*(log10_Kd_mean-d)))+c), data = Yeast_Kd_GI_Deciles_DF[Yeast_Kd_GI_Deciles_DF$GI_Score_Type == "Negative",], start = list(b = 0.26, c = -0.035, d = -12))

## Positive is a linear equation a*x+b
#Yeast_PosGI_Model <- nls(GI_mean ~ a*log10_Kd_mean+b, data = Yeast_Kd_GI_Deciles_DF[Yeast_Kd_GI_Deciles_DF$GI_Score_Type == "Positive",], start = list(b = 1, a = 1))
Yeast_PosGI_Model <- nls(GI_mean ~ (a/(1 + exp(b*(log10_Kd_mean-d)))+c), data = Yeast_Kd_GI_Deciles_DF[Yeast_Kd_GI_Deciles_DF$GI_Score_Type == "Positive",], start = list(b = 3.2, a = 0.32, c = 0.12, d = -6))


## Get parameters for negative GI in Yeast
Yeast_NegGI_Model_sum <- summary(Yeast_NegGI_Model)
Y_NegGI_a <- Yeast_NegGI_Model_sum$parameters['a','Estimate']
Y_NegGI_b <- Yeast_NegGI_Model_sum$parameters['b','Estimate']
Y_NegGI_c <- Yeast_NegGI_Model_sum$parameters['c','Estimate']
Y_NegGI_d <- Yeast_NegGI_Model_sum$parameters['d','Estimate']

## Get parameters for positive GI in Yeast
Yeast_PosGI_Model_sum <- summary(Yeast_PosGI_Model)
Y_PosGI_a <- Yeast_PosGI_Model_sum$parameters['a','Estimate']
Y_PosGI_b <- Yeast_PosGI_Model_sum$parameters['b','Estimate']
Y_PosGI_c <- Yeast_PosGI_Model_sum$parameters['c','Estimate']
Y_PosGI_d <- Yeast_PosGI_Model_sum$parameters['d','Estimate']

## Both positive and negative
Yeast_GI_Kd_lin_plot <-
  ggplot(data = Yeast_Kd_GI_Deciles_DF,
         mapping = aes(x = Kd_mean,
                       y = GI_mean,
                       color = GI_Score_Type)) +
  geom_point(size = geompoint_size) +
  ggpubr::stat_cor(method="spearman",
                   show.legend = FALSE,
                   cor.coef.name = "rho",
                   alternative = "less",
                   label.x = -6.5,
                   size = 2) +
  geom_function(data = Yeast_Kd_GI_Deciles_DF[Yeast_Kd_GI_Deciles_DF$GI_Score_Type == "Negative",],
                fun = function(x) (Y_NegGI_a/(1 + exp(Y_NegGI_b*(log10(x)-(Y_NegGI_d))))) + Y_NegGI_c,
                linewidth = geomline_linewidth) +
  geom_function(data = Yeast_Kd_GI_Deciles_DF[Yeast_Kd_GI_Deciles_DF$GI_Score_Type == "Positive",],
                fun = function(x) Y_PosGI_a*log10(x) + Y_PosGI_b,
                linewidth = geomline_linewidth) +
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
    legend.position="bottom",
    legend.background = element_rect(fill='transparent'),
    legend.text = element_text(size = Axis_text_size),
    legend.title = element_blank(), #element_text(colour = "black", size = Axis_title_size),
    plot.title = element_text(size = plot_title_size,
                              hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05),
                     limits = c(0, 0.25)) +
  scale_x_log10(limits = c(10^-8, 10^-3),
                breaks = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^"-8"),
                           expression(10^"-7"), expression(10^"-6"),
                           expression(10^"-5"), expression(10^"-4"),
                           expression(10^"-3"))) +
  scale_color_manual(labels = c("Negative GI", "Positive GI"),
                     values = c("blue", "#F5E700")) +
  xlab("Kd (M)") +
  ylab("|Mean GI Score|") +
  ggtitle("Yeast")

#Yeast_GI_Kd_lin_plot

## Negative GI
Yeast_Kd_negGI_Deciles_DF <- Yeast_Kd_GI_Deciles_DF[Yeast_Kd_GI_Deciles_DF$GI_Score_Type == "Negative",]
Yeast_negGI_Kd_lin_plot <-
  ggplot(data = Yeast_Kd_negGI_Deciles_DF,
         mapping = aes(x = Kd_mean,
                       y = GI_mean)) +
  geom_point(size = geompoint_size,
             color = "blue") +
  ggpubr::stat_cor(method="spearman",
                   show.legend = FALSE,
                   cor.coef.name = "rho",
                   alternative = "less",
                   label.x = -6.5,
                   size = 2) +
  geom_function(data = Yeast_Kd_negGI_Deciles_DF,
                fun = function(x) (Y_NegGI_a/(1 + exp(Y_NegGI_b*(log10(x)-(Y_NegGI_d))))) + Y_NegGI_c,
                linewidth = geomline_linewidth,
                color = "blue") +
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
  scale_y_continuous(breaks = seq(0.05, 0.25, by = 0.05),
                     limits = c(0.05, 0.25)) +
  scale_x_log10(limits = c(10^-8, 10^-3),
                breaks = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^"-8"),
                           expression(10^"-7"), expression(10^"-6"),
                           expression(10^"-5"), expression(10^"-4"),
                           expression(10^"-3"))) +
  scale_color_manual(labels = c("Negative GI"),
                     values = c("blue")) +
  xlab("Kd (M)") +
  ylab("|Mean GI Score|") +
  ggtitle("Yeast")

#Yeast_negGI_Kd_lin_plot

## Positive GI
Yeast_Kd_posGI_Deciles_DF <- Yeast_Kd_GI_Deciles_DF[Yeast_Kd_GI_Deciles_DF$GI_Score_Type == "Positive",]
Yeast_posGI_Kd_lin_plot <-
  ggplot(data = Yeast_Kd_posGI_Deciles_DF,
         mapping = aes(x = Kd_mean,
                       y = GI_mean)) +
  geom_point(size = geompoint_size,
             color = "#F5E700") +
  ggpubr::stat_cor(method="spearman",
                   show.legend = FALSE,
                   cor.coef.name = "rho",
                   alternative = "less",
                   label.x = -6.5,
                   size = 2) +
  geom_function(data = Yeast_Kd_posGI_Deciles_DF,
                fun = function(x) (Y_PosGI_a/(1 + exp(Y_PosGI_b*(log10(x)-(Y_PosGI_d))))) + Y_PosGI_c,
                linewidth = geomline_linewidth,
                color = "#F5E700") +
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
  scale_y_continuous(breaks = seq(0.04, 0.08, by = 0.01),
                     limits = c(0.04, 0.08)) +
  scale_x_log10(limits = c(10^-8, 10^-3),
                breaks = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^"-8"),
                           expression(10^"-7"), expression(10^"-6"),
                           expression(10^"-5"), expression(10^"-4"),
                           expression(10^"-3"))) +
  scale_color_manual(labels = c("Negative GI"),
                     values = c("blue")) +
  xlab("Kd (M)") +
  ylab("|Mean GI Score|") +
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


#### Human GI Kd Lineplot ####

## Negative is a sigmoidal function of general form (1/(1 + e^x))
## Specific form is (a/(1 + e^b*(x-d)) + c)
Human_NegGI_Model <- nls(MeanGI ~ (a/(1 + exp(b*(log10_Kd_mean-d)))+c), data = Human_Kd_GI_Quintiles_DF[Human_Kd_GI_Quintiles_DF$GI_Type == "Negative",], start = list(b = 5.2, a = 0.8, c = 1.1, d = -4.8))

## Positive is a linear equation a*x+b
Human_PosGI_Model <- nls(MeanGI ~ a*log10_Kd_mean+b, data = Human_Kd_GI_Quintiles_DF[Human_Kd_GI_Quintiles_DF$GI_Type == "Positive",], start = list(b = 1, a = 1))

## Get parameters for negative GI in Human Cell Line
Human_NegGI_Model_sum <- summary(Human_NegGI_Model)
H_NegGI_a <- Human_NegGI_Model_sum$parameters['a','Estimate']
H_NegGI_b <- Human_NegGI_Model_sum$parameters['b','Estimate']
H_NegGI_c <- Human_NegGI_Model_sum$parameters['c','Estimate']
H_NegGI_d <- Human_NegGI_Model_sum$parameters['d','Estimate']

## Get parameters for positive GI in Human Cell Line
Human_PosGI_Model_sum <- summary(Human_PosGI_Model)
H_PosGI_a <- Human_PosGI_Model_sum$parameters['a','Estimate']
H_PosGI_b <- Human_PosGI_Model_sum$parameters['b','Estimate']

Human_GI_Kd_lin_plot <- 
  ggplot(data = Human_Kd_GI_Quintiles_DF,
         mapping = aes(x = Kd_mean,
                       y = MeanGI,
                       color = GI_Type)) +
  geom_point(size = geompoint_size) +
  ggpubr::stat_cor(method="spearman",
                   cor.coef.name = "rho",
                   alternative = "less",
                   label.x = -5.25,
                   size = 2,
                   show.legend = FALSE) +
  geom_function(data = Human_Kd_GI_Quintiles_DF[Human_Kd_GI_Quintiles_DF$GI_Type == "Negative",],
                fun = function(x) (H_NegGI_a/(1 + exp(H_NegGI_b*(log10(x)-(H_NegGI_d))) )) + H_NegGI_c,
                linewidth = geomline_linewidth) +
  geom_function(data = Human_Kd_GI_Quintiles_DF[Human_Kd_GI_Quintiles_DF$GI_Type == "Positive",],
                fun = function(x) H_PosGI_a*log10(x) + H_PosGI_b,
                linewidth = geomline_linewidth) +
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
  scale_x_log10(breaks = c(10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^"-7"), 
                           expression(10^"-6"), expression(10^"-5"), 
                           expression(10^"-4"), expression(10^"-3"))) +
  #scale_y_continuous(limits = c(0.8, 1.85),
  #                   breaks = c(0, 0.05, 0.1, 0.15, 0.2),
  #                   labels = c(0, 0.05, 0.1, 0.15, 0.2)) +
  xlab("Kd (M)") +
  ylab("|Mean GI Score|") +
  scale_color_manual(labels = c("Negative GI", "Positive GI"),
                     values = c("blue", "#F5E700")) +
  ggtitle("Human")

#Human_GI_Kd_lin_plot


## Negative GI
Human_Kd_negGI_Quintiles_DF <- Human_Kd_GI_Quintiles_DF[Human_Kd_GI_Quintiles_DF$GI_Type == "Negative",]
Human_negGI_Kd_lin_plot <- 
  ggplot(data = Human_Kd_negGI_Quintiles_DF,
         mapping = aes(x = Kd_mean,
                       y = MeanGI)) +
  geom_point(size = geompoint_size,
             color = "blue") +
  ggpubr::stat_cor(method="spearman",
                   cor.coef.name = "rho",
                   alternative = "less",
                   label.x = -5.25,
                   size = 2,
                   show.legend = FALSE) +
  geom_function(data = Human_Kd_negGI_Quintiles_DF,
                fun = function(x) (H_NegGI_a/(1 + exp(H_NegGI_b*(log10(x)-(H_NegGI_d))) )) + H_NegGI_c,
                linewidth = geomline_linewidth,
                color = "blue") +
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
  scale_x_log10(breaks = c(10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^"-7"), 
                           expression(10^"-6"), expression(10^"-5"), 
                           expression(10^"-4"), expression(10^"-3"))) +
  scale_y_continuous(limits = c(1, 1.85),
                     breaks = c(1, 1.25, 1.5, 1.75, 2),
                     labels = c(1, 1.25, 1.5, 1.75, 2)) +
  xlab("Kd (M)") +
  ylab("|Mean GI Score|") +
  ggtitle("Human")

#Human_negGI_Kd_lin_plot


## Positive GI
Human_Kd_posGI_Quintiles_DF <- Human_Kd_GI_Quintiles_DF[Human_Kd_GI_Quintiles_DF$GI_Type == "Positive",]
Human_posGI_Kd_lin_plot <- 
  ggplot(data = Human_Kd_posGI_Quintiles_DF,
         mapping = aes(x = Kd_mean,
                       y = MeanGI)) +
  geom_point(size = geompoint_size,
             color = "#F5E700") +
  ggpubr::stat_cor(method="spearman",
                   cor.coef.name = "rho",
                   alternative = "less",
                   label.x = -5.25,
                   size = 2,
                   show.legend = FALSE) +
  geom_function(data = Human_Kd_GI_Quintiles_DF[Human_Kd_GI_Quintiles_DF$GI_Type == "Positive",],
                fun = function(x) H_PosGI_a*log10(x) + H_PosGI_b,
                linewidth = geomline_linewidth,
                color = "#F5E700") +
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
  scale_x_log10(breaks = c(10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^"-7"), 
                           expression(10^"-6"), expression(10^"-5"), 
                           expression(10^"-4"), expression(10^"-3"))) +
  scale_y_continuous(limits = c(0.8, 1.2),
                     breaks = c(0.8, 0.9, 1, 1.1, 1.2),
                     labels = c(0.8, 0.9, 1, 1.1, 1.2)) +
  xlab("Kd (M)") +
  ylab("|Mean GI Score|") +
  ggtitle("Human")

#Human_posGI_Kd_lin_plot


#### Combine plots ####

# Get legends
yeastGILegend <- ggpubr::get_legend(Yeast_GI_Kd_lin_plot)

# Remove legends from yeast plots
Yeast_GI_Kd_lin_plot_noLeg <- Yeast_GI_Kd_lin_plot + theme(legend.position = "none")


TopPlot2 <-
  plot_grid(Yeast_negGI_Kd_lin_plot,
            Human_negGI_Kd_lin_plot,
            ncol = 2,
            labels = c("a", "b"),
            label_size = 8,
            label_fontface = "bold",
            rel_widths = c(1,1))

MidPlot2 <-
  plot_grid(Yeast_posGI_Kd_lin_plot,
            Human_posGI_Kd_lin_plot,
            ncol = 2,
            labels = c("c", "d"),
            label_size = 8,
            label_fontface = "bold",
            rel_widths = c(1,1))

legend <-
  plot_grid(yeastGILegend,
            ncol = 1,
            labels = NULL)

BotPlot2 <-
  plot_grid(NULL,
            NULL,
            ncol = 2,
            labels = c("e", "f"),
            label_size = 8,
            label_fontface = "bold")

CombinedPlot <-
  plot_grid(TopPlot2,
            MidPlot2,
            legend,
            BotPlot2,
            ncol = 1,
            rel_heights = c(1,1,0.1,1)
            #rel_heights = c(0.03, 1, 0.03, 1),
            #labels = c("Yeast", "Human"),
            #label_size = 8,
            #label_fontface = "bold",
            #label_x = 0.4,
            #label_y = 1.05)
  )


png("../results/Figures/Main_Figures/Figure_2_properdim.png", units="mm", width=89, height=120, res=300)
CombinedPlot
dev.off()

mm_to_inch = 0.0393701
svg("../results/Figures/Main_Figures/Figure_2_properdim.svg", width=89*mm_to_inch, height=120*mm_to_inch)
CombinedPlot
dev.off()
