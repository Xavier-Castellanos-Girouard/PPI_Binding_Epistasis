# Xavier Castellanos-Girouard
# 
#
# Date First Created: Adapted from a script created sometime Summer 2022
# Date Last Modified: May 30 2024

#### Import libraries ####

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggsignif)


#### Import and format data from Horlbeck ####

# Importing GI data from Horlbeck et al. (Cell 2018) supplementary data S12
# The "gene GI scores and correlations" is the main interest for this analysis
GI_Score_Corr_Network_DF <- 
  read_excel("../data/GI_data/GI_Horlbeck.xlsx", 
             sheet = "gene GI scores and correlations",
             col_types = c("text", "text", rep('numeric', 12)))

# Retreiving only GI scores (replicate means). Cell lines separated into distinct datasets.
K562_GI_Score <- GI_Score_Corr_Network_DF[4:nrow(GI_Score_Corr_Network_DF), c(1,2,7)]
Jurkat_GI_Score <- GI_Score_Corr_Network_DF[4:nrow(GI_Score_Corr_Network_DF), c(1,2,13)]


#Renaming columns for clarity
colnames(K562_GI_Score) <- c("Gene_1", "Gene_2", "GI_Score")
colnames(Jurkat_GI_Score) <- c("Gene_1", "Gene_2", "GI_Score")


#### Import and format Stoichiometry data from Hein and Opencell ####


### Import data from opencell
OpenCell_Stoich_PPI_DF <- read.csv("../data/Mass_Spectrometry_Interactomes/opencell-protein-interactions.csv")

# Remove infinite, zero, and NA values
OpenCell_Stoich_PPI_DF <- 
  OpenCell_Stoich_PPI_DF %>%
  dplyr::filter(!is.na(interaction_stoichiometry) & !is.na(abundance_stoichiometry)) %>%
  dplyr::filter((interaction_stoichiometry!= Inf) & (interaction_stoichiometry!= 0)) %>%
  dplyr::filter(abundance_stoichiometry != 0)

# Remove columns with no use
OpenCell_Stoich_PPI_DF <- 
  OpenCell_Stoich_PPI_DF %>%
  dplyr::select("target_gene_name", "interactor_gene_name",
                "interaction_stoichiometry", "abundance_stoichiometry")

# Rename columns for compatibility with Hein data
colnames(OpenCell_Stoich_PPI_DF) <- c("Bait", "Prey", "Interaction_Stoichiometry", "Abundance_Stoichiometry")

# Certain Preys will have multiple possible genes. Separate for merge with
#    the Horlbeck Dataframe
OpenCell_Stoich_PPI_DF <- 
  OpenCell_Stoich_PPI_DF %>%
  tidyr::separate_longer_delim(cols = "Bait",
                               delim = ";") %>%
  tidyr::separate_longer_delim(cols = "Prey",
                               delim = ";")

## Reset indices and add column to indicate origin of data
row.names(OpenCell_Stoich_PPI_DF) <- NULL
OpenCell_Stoich_PPI_DF$Dataset <- "OpenCell"

### Import PPI stoichiometries from Hein et al. data.
Hein_Stoich_PPI_DF <- read_excel("../data/Mass_Spectrometry_Interactomes/Stoichiometries_Hein.xlsx", sheet = "interactions")

# Stoichiometry is currently presented as Log(Prey/Bait).
#Log10.bait.prey.ratio is from DNA pull-down, (Interaction Stoichiometry)
#log10.bait.prey.expression.ratio is from HeLa proteome (co-IP w/ LC-MS) (Abundance Stoichiometry)

Hein_Stoich_PPI_DF <- Hein_Stoich_PPI_DF[,c(5,6,12,13)]

## Certain Preys will have multiple possible genes. Separate for merge with
##    the Horlbeck Dataframe
Hein_Stoich_PPI_DF <- 
  Hein_Stoich_PPI_DF %>%
  tidyr::separate_longer_delim(cols = "prey.Gene.name",
                               delim = ";") %>%
  tidyr::separate_longer_delim(cols = "bait.Gene.name",
                               delim = ";")

# Rename columns so Stoichiometry type is clear
colnames(Hein_Stoich_PPI_DF) <- c("Bait", "Prey", 
                                  "Interaction_Stoichiometry",
                                  "Abundance_Stoichiometry")

# Store Stoichimetries as numeric type
Hein_Stoich_PPI_DF$Interaction_Stoichiometry <- as.numeric(Hein_Stoich_PPI_DF$Interaction_Stoichiometry)
Hein_Stoich_PPI_DF$Abundance_Stoichiometry <- as.numeric(Hein_Stoich_PPI_DF$Abundance_Stoichiometry)

# Keep Stoichiometries in linear space (non-log)
Hein_Stoich_PPI_DF$Interaction_Stoichiometry <- 10^Hein_Stoich_PPI_DF$Interaction_Stoichiometry
Hein_Stoich_PPI_DF$Abundance_Stoichiometry <- 10^Hein_Stoich_PPI_DF$Abundance_Stoichiometry

## Filter Data from self interacting proteins creating peaks at 1
Hein_Stoich_PPI_DF <- Hein_Stoich_PPI_DF[Hein_Stoich_PPI_DF$Interaction_Stoichiometry != 1,]

# Remove null values
Hein_Stoich_PPI_DF <- 
  Hein_Stoich_PPI_DF %>%
  dplyr::filter(!is.na(Interaction_Stoichiometry))

## Deal with duplicates
Hein_Stoich_PPI_DF <-
  Hein_Stoich_PPI_DF %>%
  dplyr::group_by(Bait, Prey) %>%
  dplyr::summarise(Interaction_Stoichiometry = mean(Interaction_Stoichiometry),
                   Abundance_Stoichiometry = mean(Abundance_Stoichiometry, na.rm = TRUE))

# Reset row names and indicate source of data
row.names(Hein_Stoich_PPI_DF) <- NULL
Hein_Stoich_PPI_DF$Dataset <- "Hein"

### Merge Hein data and Opencell data
Human_Stoich_PPI_DF <- rbind(Hein_Stoich_PPI_DF, OpenCell_Stoich_PPI_DF)

#### Merge Stoichiometry data with Jurkat Horlbeck data ####

#Merge of the two dataframes by gene pair:
merge_combo1 <- 
  merge(x = Human_Stoich_PPI_DF,
        y = Jurkat_GI_Score,
        by.x = c("Bait", "Prey"),
        by.y = c("Gene_1", "Gene_2"),
        all.x = FALSE,
        all.y = FALSE)

length((merge_combo1$GI_Score[!(is.na(merge_combo1$GI_Score))]))
# 322 common interactions.


# Second merge combination:
merge_combo2 <- 
  merge(x = Human_Stoich_PPI_DF,
        y = Jurkat_GI_Score, 
        by.x = c("Bait", "Prey"),
        by.y = c("Gene_2", "Gene_1"),
        all.x = FALSE,
        all.y = FALSE)

length((merge_combo2$GI_Score[!(is.na(merge_combo2$GI_Score))]))
# 301 common interactions for this combination.

# Total of 623 Common interactions.
Human_Stoich_GI_DF <- rbind(merge_combo1, merge_combo2)

Human_Stoich_GI_DF <- Human_Stoich_GI_DF[!is.na(Human_Stoich_GI_DF$GI_Score),]


#### Human Stoichiometry StatPlot with Mean GIs ####

## Positive GI stat Plot
Human_PosGI_Stoichiometry_plot <-
  ggplot(data = Human_Stoich_GI_DF[(Human_Stoich_GI_DF$GI_Score > 0),],
         mapping = aes(x = Interaction_Stoichiometry, 
                       y = Abundance_Stoichiometry,
                       z = GI_Score)) +
  stat_summary_2d(bins = 16, fun = "mean") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = 0.5),
        axis.line.y = element_line(color ="black", linewidth = 0.5),
        axis.ticks = element_line(color = "black"),
        axis.text.x.bottom = element_text(color = "black", size = 12),
        axis.text.y.left = element_text(color = "black", size = 12),
        axis.title.x.bottom = element_text(color = "black", size = 14),
        axis.title.y.left = element_text(color = "black", size = 14),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    linewidth=0.5),
        #legend.text.align = 0,
        #legend.key = element_rect(fill = "white"),
        #legend.position = "none",
        legend.key = element_rect(fill = "black", linewidth = 4)) +
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
                                               ticks = TRUE,
                                               ticks.colour='black')) +
  xlab("Interaction Stoichiometry") +
  ylab("Abundance\nStoichiometry") +
  labs(fill = '   Mean\nGI Score\n') +
  ggtitle("Positive GI")

#Human_PosGI_Stoichiometry_plot

## Negative GI Stat Plot
Human_NegGI_Stoichiometry_plot <-
  ggplot(data = Human_Stoich_GI_DF[(Human_Stoich_GI_DF$GI_Score < 0),],
         mapping = aes(x = Interaction_Stoichiometry, 
                       y = Abundance_Stoichiometry,
                       z = abs(GI_Score))) +
  stat_summary_2d(bins = 16, fun = "mean") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = 0.5),
        axis.line.y = element_line(color ="black", linewidth = 0.5),
        axis.ticks = element_line(color = "black"),
        axis.text.x.bottom = element_text(color = "black", size = 12),
        axis.text.y.left = element_text(color = "black", size = 12),
        axis.title.x.bottom = element_text(color = "black", size = 14),
        axis.title.y.left = element_text(color = "black", size = 14),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    linewidth=0.5),
        #legend.text.align = 0,
        #legend.key = element_rect(fill = "white"),
        #legend.position = "none",
        legend.key = element_rect(fill = "black", linewidth = 4)) +
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
                                               ticks = TRUE,
                                               ticks.colour='black')) +
  xlab("Interaction Stoichiometry") +
  ylab("Abundance\nStoichiometry") +
  labs(fill = '   Mean\nGI Score\n') +
  ggtitle("Negative GI")

#Human_NegGI_Stoichiometry_plot

### Divide into Hein regions (circle)

# Make new factor column for regions
Human_Stoich_GI_DF$Region <- NA


## Region 1 Subset. Radius 0.75, Center (-0.5, 0) ; (IntStoich, AbStoich)

# For every interaction calculate distance from circle center of region 1
Reg1_distances <- 
  sqrt((log10(Human_Stoich_GI_DF$Interaction_Stoichiometry) - (-0.5))^2 +
         (log10(Human_Stoich_GI_DF$Abundance_Stoichiometry) - 0)^2)

# Assign region 1 to any point within the circle (distance lower than 1.0)
Human_Stoich_GI_DF$Region[Reg1_distances < 0.75] <- 1 # Assign region 1 to any 


# Region 2 Subset. Center (-3.5, 0) 
Reg2_distances <- 
  sqrt((log10(Human_Stoich_GI_DF$Interaction_Stoichiometry) - (-3.5))^2 +
         (log10(Human_Stoich_GI_DF$Abundance_Stoichiometry) - 0)^2)

Human_Stoich_GI_DF$Region[Reg2_distances < 0.75] <- 2

# Region 3 Subset. Center (-0.5, 1.5) 
Reg3_distances <- 
  sqrt((log10(Human_Stoich_GI_DF$Interaction_Stoichiometry) - (-0.5))^2 +
         (log10(Human_Stoich_GI_DF$Abundance_Stoichiometry) - 1.5)^2)

Human_Stoich_GI_DF$Region[Reg3_distances < 0.75] <- 3

# Region 4 Subset. Center (-3.5, -1.5) 
Reg4_distances <- 
  sqrt((log10(Human_Stoich_GI_DF$Interaction_Stoichiometry) - (-3.5))^2 +
         (log10(Human_Stoich_GI_DF$Abundance_Stoichiometry) - (-1.5))^2)

Human_Stoich_GI_DF$Region[Reg4_distances < 0.75] <- 4

Human_Stoich_GI_DF <- Human_Stoich_GI_DF[!is.na(Human_Stoich_GI_DF$Region),]


#### Human Stoichiometry Region vs GI Violin pot ####

## Set to factor to conserve order in graph
Human_Stoich_GI_DF$Region <- 
  factor(Human_Stoich_GI_DF$Region, 
         levels = 1:4)

Human_PosGI_StoichReg_violin_plot <-
  ggplot(data = Human_Stoich_GI_DF[(Human_Stoich_GI_DF$GI_Score > 0),],
         mapping = aes(x = Region,
                       y = abs(GI_Score))) +
  geom_violin(trim = TRUE,
              draw_quantiles = c(0.25, 0.5, 0.75),
              fill = "yellow2",
              color = "black") +
  geom_signif(comparisons = list(c(2,3)),
              tip_length = 0.01,
              step_increase = 0.06,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color ="black", linewidth = 0.5),
    axis.line.y = element_line(color ="black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text.x.bottom = element_text(color = "black", size = 12),
    axis.text.y.left = element_text(color = "black", size = 12),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_y_log10(breaks = c(10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-3"),
                           expression(10^"-2"), expression(10^"-1"), 
                           expression(10^0), expression(10^1),
                           expression(10^2)),
                limits = c(10^-3, 10^2)) +
  scale_x_discrete(labels = c("Reg. 1", "Reg. 2", "Reg. 3", "Reg. 4")) +
  xlab("Stoichiometry Region") +
  ylab("Positive GI")

#Human_PosGI_StoichReg_violin_plot



Human_NegGI_StoichReg_violin_plot <-
  ggplot(data = Human_Stoich_GI_DF[(Human_Stoich_GI_DF$GI_Score < 0),],
         mapping = aes(x = Region,
                       y = abs(GI_Score))) +
  geom_violin(trim = TRUE,
              draw_quantiles = c(0.25, 0.5, 0.75),
              fill = "royalblue2",
              color = "black") +
  geom_signif(comparisons = list(c(1,2), c(1,3)),
              tip_length = 0.01,
              step_increase = 0.06,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color ="black", linewidth = 0.5),
    axis.line.y = element_line(color ="black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text.x.bottom = element_text(color = "black", size = 12),
    axis.text.y.left = element_text(color = "black", size = 12),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_y_log10(breaks = c(10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression("-"*10^"-3"),
                           expression("-"*10^"-2"), expression("-"*10^"-1"), 
                           expression("-"*10^0), expression("-"*10^1),
                           expression("-"*10^2)),
                limits = c(10^-3, 10^2)) +
  scale_x_discrete(labels = c("Reg. 1", "Reg. 2", "Reg. 3", "Reg. 4")) +
  xlab("Stoichiometry Region") +
  ylab("Negative GI")

#Human_NegGI_StoichReg_violin_plot



#### Combined plot Display ####


# Yeast and Human Combined
CombinedPlot <-
  plot_grid(Human_NegGI_Stoichiometry_plot,
            Human_PosGI_Stoichiometry_plot,
            Human_NegGI_StoichReg_violin_plot,
            Human_PosGI_StoichReg_violin_plot,
            labels = c("A", "B", "C", "D"),
            rel_heights = c(1, 1, 1, 1),
            ncol = 2)

png("../results/Figures/Extended_Data_Figures/SupplementaryFigure3.png", units="mm", width=180, height=120, res=300)
CombinedPlot
dev.off()
