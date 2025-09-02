# Xavier Castellanos-Girouard
#
# Statistical Controls for Network Features
#
# Date First Created: June 5 2024
# Date Last Modified: August 7th 2025

#### Import Libraries ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(readxl)

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

#### Import Data ####

### Modules
ModuleOverlap_Control_DF <- read.csv("../results/Topological_Analyses/GI_PPI_randomized_ModuleOverlap.csv", row.names = 1)
rownames(ModuleOverlap_Control_DF) <- NULL
colnames(ModuleOverlap_Control_DF)[1] <- "Cluster_ID"

ModuleOverlap_DF <- read.csv("../results/Topological_Analyses/GI_PPI_optimal_module_overlap.csv", row.names = 1)


### Connectors
ConnectorOverlap_Control_DF <- read.csv("../data/Precomputed_data/GI_PPI_ConnectorOverlap_Controls.csv")

ConnectorOverlap_DF <- read.csv("../results/Topological_Analyses/Shortest_Path_Connector_Overlap.csv", row.names = 1)
row.names(ConnectorOverlap_DF) <- NULL

#### Overall Module Overlap Analysis ####

mean_Jaccard_scores <- c()
for (i in unique(ModuleOverlap_Control_DF$iter_ID)){
  mean_Jaccard <- mean(ModuleOverlap_Control_DF$Jaccard_index[ModuleOverlap_Control_DF$iter_ID==i])
  mean_Jaccard_scores <- c(mean_Jaccard_scores, mean_Jaccard)
}


## Stats
zscore = (mean(ModuleOverlap_DF$Jaccard_index) - mean(mean_Jaccard_scores))/sd(mean_Jaccard_scores)
pval = pnorm(zscore, mean = 0, sd = 1, lower.tail = FALSE)

module_exp_vs_ctrl_p <- 
  ggplot(data=NULL,
         mapping = aes(x = mean_Jaccard_scores)) +
  geom_histogram(fill = "black",
                 color = "white",
                 bins = 50) +
  geom_vline(xintercept = mean(ModuleOverlap_DF$Jaccard_index),
             color = "red",
             linewidth = stoichplot_dashed_linewidth) +
  annotate("text", 
           x=0.20, 
           y=1600, 
           label = paste0("p-value = ",
                          signif(pval, 
                                 digits=3)),
           size = 2) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = Axis_text_size), 
    axis.text.y = element_text(colour="black", size = Axis_text_size),
    axis.ticks = element_line(color="black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.line.x.bottom=element_line(color="black", linewidth = xy_plotline_width),
    axis.line.y.left=element_line(color="black", linewidth = xy_plotline_width),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position = "none") +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 0.375),
                     breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("Mean Jaccard Score") +
  ylab("Count")

## Stats
module_exp_vs_ctrl_p

#### Overall Connector Overlap Analysis ####

# Find instances where a new Control starts (i.e. index is 0)
DF_starts <- rev(which(ConnectorOverlap_Control_DF$X == 0))

number_of_overlapping_connectors <- c()
init = nrow(ConnectorOverlap_Control_DF)+1
for (i in DF_starts){
  iter_length <- init - i # length of this iteration of controls for connector overlaps
  number_of_overlapping_connectors <- c(number_of_overlapping_connectors, iter_length)
  init = init-iter_length
}

## Visualize 

connector_zscore = (nrow(ConnectorOverlap_DF) - mean(number_of_overlapping_connectors))/sd(number_of_overlapping_connectors)

connector_pval = pnorm(connector_zscore, mean = 0, sd = 1, lower.tail = FALSE)


connector_exp_vs_ctrl_p <- 
  ggplot(data=NULL,
         mapping = aes(x = number_of_overlapping_connectors)) +
  geom_histogram(fill = "black",
                 color = "white",
                 bins = 50) +
  geom_vline(xintercept = nrow(ConnectorOverlap_DF),
             color = "red",
             linewidth = stoichplot_dashed_linewidth) +
  annotate("text", 
           x=400, 
           y=115, 
           label = paste0("p-value = ",
                          signif(connector_pval, 
                                 digits=3)),
           size = 2) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = Axis_text_size), 
    axis.text.y = element_text(colour="black", size = Axis_text_size),
    axis.ticks = element_line(color="black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.line.x.bottom=element_line(color="black", linewidth = xy_plotline_width),
    axis.line.y.left=element_line(color="black", linewidth = xy_plotline_width),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position = "none") +
  scale_x_continuous(expand = c(0,0),
                     limits = c(50, 750),
                     breaks = c(100, 200, 300, 400, 500, 600, 700)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("Number of Overlapping Shortest Paths") +
  ylab("Count")


#### Make plots ####


Combined_plot <-
  plot_grid(module_exp_vs_ctrl_p,
            connector_exp_vs_ctrl_p,
            labels = c("a", "b"),
            rel_widths = c(1, 1),
            ncol = 2,
            label_size = 8,
            label_fontface = "bold")

png("../results/Figures/Extended_Data_Figures/SupplementaryFigure12.png", units = "mm", res = 300, width = 120, height = 60)
Combined_plot
dev.off()

#mm_to_inch = 0.0393701
#svg("../../../graphs/ModuleConnector_Overlap_Ctrls_properdim.svg", width=120*mm_to_inch, height=60*mm_to_inch)
#Combined_plot
#dev.off()