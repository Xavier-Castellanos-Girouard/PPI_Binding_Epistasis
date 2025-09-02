# Xavier Castellanos-Girouard

# Date First Created: March 14 2024 
# Date Last Modified: May 29 2024


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

#### All Ka Distribution plot ####


ka_distribution_p <-
  ggplot(data = all_paralogs_kaks_DF,
       mapping = aes(x = ka)) +
  geom_histogram(bins = 20,
                 color = "white",
                 fill = "black") +
  geom_vline(xintercept = 0.2, 
             color = "red",
             linewidth = 1) +
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
    axis.title.y.left = element_text(color = "black", size = 14)) +
  scale_x_continuous() +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Sequence Divergence (dN)") +
  ylab("Count")

#ka_distribution_p

#### P-value robustness (One Paralog) ####

# For One paralog, only one ka value is needed
OnePara_Ka_DF <- 
  Kd_GI_AllPara_DF %>%
  dplyr::filter((Type == "OneParalog") & (!is.na(source_Ka) | !is.na(target_Ka)))

# Necessary to avoid bugs
OnePara_Ka_DF$ka_merged <- pmax(OnePara_Ka_DF$source_Ka, OnePara_Ka_DF$target_Ka, na.rm = TRUE)

onePara_pval_vec <- c()
for (i in seq(0.001, 1, 0.001)) {
  OnePara_Ka_DF$Ka_div <- NA
  
  Ka_thresh <- i
  
  # Assign High and low
  OnePara_Ka_DF$Ka_div[(OnePara_Ka_DF$ka_merged > Ka_thresh)] <- "High"
  OnePara_Ka_DF$Ka_div[(OnePara_Ka_DF$ka_merged < Ka_thresh)] <- "Low"
  
  wilcox_test_pval <- wilcox.test(OnePara_Ka_DF$Kd[OnePara_Ka_DF$Ka_div == "High"], OnePara_Ka_DF$Kd[OnePara_Ka_DF$Ka_div == "Low"])$p.value
  
  onePara_pval_vec <- c(onePara_pval_vec, wilcox_test_pval)
}

## Make plot
OnePara_SeqDiv_p <- 
  ggplot(data = NULL,
         mapping = aes(x = seq(0.001, 1, 0.001),
                       y = -log10(onePara_pval_vec))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), color = 'red') +
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
    plot.title = element_text(hjust = 0.5)) +
  #legend.position = "none") +
  xlab("Sequence Divergence (dN) Threshold") +
  ylab("-log10 p-value") +
  ggtitle("One Paralog")

#OnePara_SeqDiv_p

#### P-value robustness (Two Paralog) ####

# For Two paralog, both ka values are needed. Filter NA values
TwoPara_Ka_DF <- 
  Kd_GI_AllPara_DF %>%
  dplyr::filter((Type == "TwoParalog") & (!is.na(source_Ka) & !is.na(target_Ka)))

twoPara_pval_vec <- c()
for (i in seq(0.001, 0.75, 0.001)) {
  TwoPara_Ka_DF$Ka_div <- NA
  
  Ka_thresh <- i
  
  # Assign High and low
  TwoPara_Ka_DF$Ka_div[((TwoPara_Ka_DF$source_Ka > Ka_thresh) & (TwoPara_Ka_DF$target_Ka > Ka_thresh))] <- "High"
  TwoPara_Ka_DF$Ka_div[((TwoPara_Ka_DF$source_Ka < Ka_thresh) | (TwoPara_Ka_DF$target_Ka < Ka_thresh))] <- "Low"
  
  wilcox_test_pval <- wilcox.test(TwoPara_Ka_DF$Kd[TwoPara_Ka_DF$Ka_div == "High"], TwoPara_Ka_DF$Kd[TwoPara_Ka_DF$Ka_div == "Low"])$p.value
  
  twoPara_pval_vec <- c(twoPara_pval_vec, wilcox_test_pval)
}


TwoPara_SeqDiv_p <- 
  ggplot(data = NULL,
         mapping = aes(x = seq(0.001, 0.75, 0.001),
                       y = -log10(twoPara_pval_vec))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), color = "red") +
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
    plot.title = element_text(hjust = 0.5)) +
  #legend.position = "none") +
  xlab("Sequence Divergence (dN) Threshold") +
  ylab("-log10 p-value") +
  ggtitle("Two Paralogs")

#TwoPara_SeqDiv_p

#### Combine Plot ####

combined_plot <-
  plot_grid(
    ka_distribution_p,
    OnePara_SeqDiv_p,
    TwoPara_SeqDiv_p,
    ncol = 3,
    rel_widths = c(1,1,1),
    labels = c("A", "B", "C"),
    vjust = 1.5
  )

#combined_plot

png("../results/Figures/Extended_Data_Figures/SupplementaryFigure10.png", units="mm", width=300, height = 60, res = 300)
combined_plot
dev.off()
