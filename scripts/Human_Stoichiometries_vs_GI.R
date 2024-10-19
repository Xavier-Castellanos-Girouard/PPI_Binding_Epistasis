# Xavier Castellanos-Girouard

# Date First Created: Sometime Summer 2022
# Date Last Modified: April 22 2024


#### Import Libraries ####

library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

#### Import and format data from Horlbeck ####

dir <- "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/"

# Importing GI data from Horlbeck et al. (Cell 2018) supplementary data S12
# The "gene GI scores and correlations" is the main interest for this analysis
GI_Score_Corr_Network_DF <- 
  read_excel(paste0(dir, "R_Stoichiometries_vs_GI/data/GI_Horlbeck.xlsx"), 
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
OpenCell_Stoich_PPI_DF <- read.csv(paste0("R_Stoichiometries_vs_GI/data/opencell-protein-interactions.csv"))

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
Hein_Stoich_PPI_DF <- read_excel(paste0(dir, "R_Stoichiometries_vs_GI/data/Stoichiometries_Hein.xlsx"), sheet = "interactions")

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

# export
#write.csv(Human_Stoich_PPI_DF, "../results/Human_Interactome_Stoich.csv")

#### Merge Stoichiometry data with K562 Horlbeck data ####

#Merge of the two dataframes by gene pair:
merge_combo1 <- 
  merge(x = Human_Stoich_PPI_DF,
        y = K562_GI_Score,
        by.x = c("Bait", "Prey"),
        by.y = c("Gene_1", "Gene_2"),
        all.x = FALSE,
        all.y = FALSE)

length((merge_combo1$GI_Score[!(is.na(merge_combo1$GI_Score))]))
# 234 common interactions.


# Second merge combination:
merge_combo2 <- 
  merge(x = Human_Stoich_PPI_DF,
        y = K562_GI_Score, 
        by.x = c("Bait", "Prey"),
        by.y = c("Gene_2", "Gene_1"),
        all.x = FALSE,
        all.y = FALSE)

length((merge_combo2$GI_Score[!(is.na(merge_combo2$GI_Score))]))
# 246 common interactions for this combination.

# Total of 507 Common interactions.
Human_Stoich_GI_DF <- rbind(merge_combo1, merge_combo2)

Human_Stoich_GI_DF <- Human_Stoich_GI_DF[!is.na(Human_Stoich_GI_DF$GI_Score),]

#write.csv(Human_Stoich_GI_DF, "../results/Human_Stoich_GI.csv")

#### Separate Stoichiometry plot into regions described by Hein et al. ####

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

#write.csv(Human_Stoich_GI_DF, paste0(dir, "R_Stoichiometries_vs_GI/results/Human_Stoich_GI_wRegions.csv"))

#### Statistical Tests (Positive GI) ####

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 1) & (Human_Stoich_GI_DF$GI_Score > 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 2) & (Human_Stoich_GI_DF$GI_Score > 0)])
# p-value = 0.01644

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 1) & (Human_Stoich_GI_DF$GI_Score > 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 3) & (Human_Stoich_GI_DF$GI_Score > 0)])
# p-value = 0.009135

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 1) & (Human_Stoich_GI_DF$GI_Score > 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 4) & (Human_Stoich_GI_DF$GI_Score > 0)])
# p-value = 0.2327

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 2) & (Human_Stoich_GI_DF$GI_Score > 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 3) & (Human_Stoich_GI_DF$GI_Score > 0)])
# p-value = 0.6106

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 2) & (Human_Stoich_GI_DF$GI_Score > 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 4) & (Human_Stoich_GI_DF$GI_Score > 0)])
# p-value = 0.08472

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 3) & (Human_Stoich_GI_DF$GI_Score > 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 4) & (Human_Stoich_GI_DF$GI_Score > 0)])
# p-value = 0.09073

#### Statistical Tests (Negative GI) ####

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 1) & (Human_Stoich_GI_DF$GI_Score < 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 2) & (Human_Stoich_GI_DF$GI_Score < 0)])
# p-value = 0.002793

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 1) & (Human_Stoich_GI_DF$GI_Score < 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 3) & (Human_Stoich_GI_DF$GI_Score < 0)])
# p-value = 0.02846

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 1) & (Human_Stoich_GI_DF$GI_Score < 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 4) & (Human_Stoich_GI_DF$GI_Score < 0)])
# p-value = 0.08335

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 2) & (Human_Stoich_GI_DF$GI_Score < 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 3) & (Human_Stoich_GI_DF$GI_Score < 0)])
# p-value = 0.9784

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 2) & (Human_Stoich_GI_DF$GI_Score < 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 4) & (Human_Stoich_GI_DF$GI_Score < 0)])
# p-value = 0.7762

wilcox.test(Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 3) & (Human_Stoich_GI_DF$GI_Score < 0)],
            Human_Stoich_GI_DF$GI_Score[(Human_Stoich_GI_DF$Region == 4) & (Human_Stoich_GI_DF$GI_Score < 0)])
# p-value = 0.9165