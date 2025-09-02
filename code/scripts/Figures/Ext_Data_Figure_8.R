# Xavier Castellanos-Girouard

# Date First Created: May 22 2024
# Date Last Modified: October 16 2024


#### Import Libraries ####

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(cowplot)

#### Import data ####

quantile_Fitting_DF <- read.csv("../results/Kd_and_GI/Quantile_Dist_GI_Table.csv",
                                     row.names = 1)
row.names(quantile_Fitting_DF) <- NULL

GI_Kd_distance_distribution_DF <- read.csv("../results/Kd_and_GI/GI_Kd_distance_fit.csv",
                                  row.names = 1)
row.names(GI_Kd_distance_distribution_DF) <- NULL

# Import epistasis network
GI_Network_DF <- read.csv("../results/Stoichiometry_and_GI/costanzo_2016_longer_withoutReps.csv")

GI_Network_DT <- 
  setDT(x = GI_Network_DF,
        key = c("ORF_query", "ORF_array"))

# Import Recontructed network
Recontructed_GI_Network_DF <- read.csv("../results/Kd_and_GI/Clustered_Kd_inferred_GI.csv")

colnames(Recontructed_GI_Network_DF)[1] <- "Source"

Recontructed_GI_Network_longer_DF <- 
  Recontructed_GI_Network_DF %>%
  tidyr::pivot_longer(cols = 2:ncol(Recontructed_GI_Network_DF),
                      names_to = "Target")

Recontructed_GI_Network_longer_DF <- Recontructed_GI_Network_longer_DF[!is.na(Recontructed_GI_Network_longer_DF$value),]

colnames(Recontructed_GI_Network_longer_DF)[3] <- "Inferred_GI"

Recontructed_GI_Network_longer_DT <- 
  setDT(Recontructed_GI_Network_longer_DF, 
        key = c("Source", "Target"))

#### Distribution of Quantile Means and Fitting ####

GI_Kd_distance_distribution_p <- 
  ggplot(data = NULL,
         mapping = aes(x = 5:(nrow(GI_Kd_distance_distribution_DF)+4), # Quantile search starts at 5
                       y = GI_Kd_distance_distribution_DF$exponentialFit_RMSE)) +
  geom_point() +
  #geom_vline(xintercept = 25) +
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
  scale_x_continuous(expand = c(0.02, 0)) +
  xlab("Number of Quantiles") +
  ylab("Root Mean Square Error")

#png("../graphs/Estimated_vs_Lit.png", units="mm", width=120, height=80, res=300)
#GI_Kd_distance_distribution_p

#### Distribution of Quantile Means and Fitting (zoomed) ####

GI_Kd_distance_distribution_zoom_p <- 
  ggplot(data = NULL,
         mapping = aes(x = 5:(nrow(GI_Kd_distance_distribution_DF[1:200,])+4), # Quantile search starts at 5
                       y = GI_Kd_distance_distribution_DF$exponentialFit_RMSE[1:200])) +
  geom_point() +
  geom_vline(xintercept = 25,
             color = "red",
             linewidth = 1) +
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
  scale_x_continuous(expand = c(0.02, 0)) +
  xlab("Number of Quantiles") +
  ylab("Root Mean Square Error")

#png("../graphs/Estimated_vs_Lit.png", units="mm", width=120, height=80, res=300)
#GI_Kd_distance_distribution_zoom_p

#### Relation (Fit) between weight and negative GI at 25 quantiles ####

quantile_Fitting_p <- 
  ggplot(data = quantile_Fitting_DF,
         mapping = aes(x = Association,
                       y = abs(scores))) +
  geom_point() +
  stat_function(fun = function(x) 0.026271*exp(0.26976*x),
                color = "red",
                linewidth = 1) +
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
  xlab("Weight") +
  ylab("Negative GI Score\n(Absolute Value)")

#png("../graphs/quantile_Fitting.png", units="mm", width=120, height=80, res=300)
#quantile_Fitting_p
#dev.off()

#### GI explainability from PPI ####

## Merge Reconstructed Network with experimental GI

# Combo 1
Comparison_GI_Network_DT_1 <-
  merge.data.table(x = Recontructed_GI_Network_longer_DT,
                   y = GI_Network_DT,
                   by.x = c("Source", "Target"),
                   by.y = c("ORF_query", "ORF_array"),
                   all.x = FALSE,
                   all.y = FALSE)

# Combo 2
Comparison_GI_Network_DT_2 <-
  merge.data.table(x = Recontructed_GI_Network_longer_DT,
                   y = GI_Network_DT,
                   by.x = c("Target", "Source"),
                   by.y = c("ORF_query", "ORF_array"),
                   all.x = FALSE,
                   all.y = FALSE)

# Merge combinations
Comparison_GI_Network_DT <-
  rbind(Comparison_GI_Network_DT_1,
        Comparison_GI_Network_DT_2)

# Only keep negative GI score
Comparison_GI_Network_DT <- Comparison_GI_Network_DT[Comparison_GI_Network_DT$scores<0,]

# Extract pvalue
pval <- cor.test(Comparison_GI_Network_DT$Inferred_GI, 
                 Comparison_GI_Network_DT$scores)

# Calculate percentage of variance explained
(pval$estimate**2)*100

# Calculate proportion of PPIs in network
Kd_GI_Network_DF <- read.csv("../results/Kd_and_GI/Yeast_Kd_GI.csv",
                              row.names = 1)

number_of_GIs <- ((ncol(Recontructed_GI_Network_DF)-1)*nrow(Recontructed_GI_Network_DF)-1)/2

nrow(Kd_GI_Network_DF)/number_of_GIs*100


#### Combine Plots ####

combined_plot <- 
  plot_grid(quantile_Fitting_p,
            NULL,
            GI_Kd_distance_distribution_p,
            GI_Kd_distance_distribution_zoom_p,
            labels = c("A", "", "B", "C"),
            rel_widths = c(1, 1, 1, 1),
            ncol = 2)

png("../results/Figures/Extended_Data_Figures/SupplementaryFigure8.png", units="mm", width=240, height=160, res=300)
combined_plot
dev.off()
