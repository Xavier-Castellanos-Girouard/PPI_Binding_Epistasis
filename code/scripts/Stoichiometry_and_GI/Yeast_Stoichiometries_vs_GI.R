# Xavier Castellanos-Girouard
# 
# Explore the Relationship between Stoichiometries and GIs
#
# Date First Created: November 9 2023
# Date Last Modified: April 22 2024

#### Import libraries ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

#### Import Data ####

# Import epistasis network
GI_Network_DF <- read.csv("../results/Stoichiometry_and_GI/costanzo_2016_longer_withoutReps.csv")

# Import Stoichiometry Network
Stoich_PPI_Network <- read.csv("../results/Stoichiometry_and_GI/Yeast_Interactome_Stoich.csv", row.names = 1)
row.names(Stoich_PPI_Network) <- NULL

#### Merge Epistasis values and Stoichiometries ####

Stoich_PPI_Network_DT <- setDT(Stoich_PPI_Network, 
                               key = c("source",
                                       "target"))

GI_Network_DT <- setDT(GI_Network_DF,
                           key = c("ORF_query",
                                   "ORF_array"))

# Get epistasis values from interactions in Stoich network (match A:B with A:B)
Stoich_GI_Network_DT_1 <-
  merge.data.table(x = Stoich_PPI_Network_DT,
                   y = GI_Network_DT,
                   by.x = c("source",
                            "target"),
                   by.y = c("ORF_query",
                            "ORF_array"),
                   all.x = FALSE,
                   all.y = FALSE)

# Get epistasis values from interactions in Stoich network (match A:B with B:A)
Stoich_GI_Network_DT_2 <-
  merge.data.table(x = Stoich_PPI_Network_DT,
                   y = GI_Network_DT,
                   by.x = c("source",
                            "target"),
                   by.y = c("ORF_array",
                            "ORF_query"),
                   all.x = FALSE,
                   all.y = FALSE)


## Remove NA values (no match)
Stoich_GI_Network_DT_1 <- na.omit(Stoich_GI_Network_DT_1, cols = "scores")
Stoich_GI_Network_DT_2 <- na.omit(Stoich_GI_Network_DT_2, cols = "scores")

## Combine data tables
Stoich_GI_Network_DT <- 
  rbind(Stoich_GI_Network_DT_1, 
        Stoich_GI_Network_DT_2)

# Remove duplicates
Stoich_GI_Network_DT <- 
  unique(Stoich_GI_Network_DT, 
         by = c("source",
                "target"))

### remove unneeded columns and rename others

Stoich_GI_Network_DT <-
  Stoich_GI_Network_DT %>%
  dplyr::select(-c(strain_query, strain_array))


colnames(Stoich_GI_Network_DT)[7:8] <- c("GI_pval", "GI_scores")

# Export merged stoichiometry and GI datasets
write.csv(Stoich_GI_Network_DT, "../results/Stoichiometry_and_GI/Yeast_Stoich_GI.csv")


#### Separate Stoichiometry plot into regions described by Hein et al. ####

# Data reset
Stoich_GI_Network_DF <- as.data.frame(Stoich_GI_Network_DT)

### Divide into Hein regions (circle)
Stoich_GI_Network_DF$Region <- NA


## Region 1 Subset. Radius 0.75, Center (-0.5, 0) ; (IntStoich, AbStoich)

# For every interaction calculate distance from circle center of region 1
Reg1_distances <- 
  sqrt((log10(Stoich_GI_Network_DF$Interaction_Stoichiometry) - (-0.5))^2 +
         (log10(Stoich_GI_Network_DF$Abundance_Stoichiometry) - 0)^2)

# Assign region 1 to any point within the circle (distance lower than 1.0)
Stoich_GI_Network_DF$Region[Reg1_distances < 0.75] <- 1 # Assign regio 1 to any 


# Region 2 Subset. Center (-3, 0) 
Reg2_distances <- 
  sqrt((log10(Stoich_GI_Network_DF$Interaction_Stoichiometry) - (-2.5))^2 +
         (log10(Stoich_GI_Network_DF$Abundance_Stoichiometry) - 0)^2)

Stoich_GI_Network_DF$Region[Reg2_distances < 0.75] <- 2

# Region 3 Subset
Reg3_distances <- 
  sqrt((log10(Stoich_GI_Network_DF$Interaction_Stoichiometry) - (-0.5))^2 +
         (log10(Stoich_GI_Network_DF$Abundance_Stoichiometry) - 1.5)^2)

Stoich_GI_Network_DF$Region[Reg3_distances < 0.75] <- 3

# Region 4 Subset
Reg4_distances <- 
  sqrt((log10(Stoich_GI_Network_DF$Interaction_Stoichiometry) - (-2.5))^2 +
         (log10(Stoich_GI_Network_DF$Abundance_Stoichiometry) - (-1.5))^2)

Stoich_GI_Network_DF$Region[Reg4_distances < 0.75] <- 4

Stoich_GI_Network_DF <- Stoich_GI_Network_DF[!is.na(Stoich_GI_Network_DF$Region),]

write.csv(Stoich_GI_Network_DF, "../results/Stoichiometry_and_GI/Yeast_Stoich_GI_wRegions.csv")

#### Statistical Tests (Positive GI) ####

print("Statistical Tests (Positive GI)")

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 1) & (Stoich_GI_Network_DF$GI_scores > 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 2) & (Stoich_GI_Network_DF$GI_scores > 0)])
# p-value = 0.0004896

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 1) & (Stoich_GI_Network_DF$GI_scores > 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 3) & (Stoich_GI_Network_DF$GI_scores > 0)])
# p-value = 0.01316

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 1) & (Stoich_GI_Network_DF$GI_scores > 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 4) & (Stoich_GI_Network_DF$GI_scores > 0)])
# p-value = 0.1625

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 2) & (Stoich_GI_Network_DF$GI_scores > 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 3) & (Stoich_GI_Network_DF$GI_scores > 0)])
# p-value = 0.7413

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 2) & (Stoich_GI_Network_DF$GI_scores > 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 4) & (Stoich_GI_Network_DF$GI_scores > 0)])
# p-value = 0.1655

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 3) & (Stoich_GI_Network_DF$GI_scores > 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 4) & (Stoich_GI_Network_DF$GI_scores > 0)])
# p-value = 0.2505

#### Statistical Tests (Negative GI) ####

print("Statistical Tests (Negative GI)")

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 1) & (Stoich_GI_Network_DF$GI_scores < 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 2) & (Stoich_GI_Network_DF$GI_scores < 0)])
# p-value = 9.21e-10

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 1) & (Stoich_GI_Network_DF$GI_scores < 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 3) & (Stoich_GI_Network_DF$GI_scores < 0)])
# p-value = 6.596e-07

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 1) & (Stoich_GI_Network_DF$GI_scores < 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 4) & (Stoich_GI_Network_DF$GI_scores < 0)])
# p-value = 1.511e-15

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 2) & (Stoich_GI_Network_DF$GI_scores < 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 3) & (Stoich_GI_Network_DF$GI_scores < 0)])
# p-value = 0.006832

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 2) & (Stoich_GI_Network_DF$GI_scores < 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 4) & (Stoich_GI_Network_DF$GI_scores < 0)])
# p-value = 0.0002067

wilcox.test(Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 3) & (Stoich_GI_Network_DF$GI_scores < 0)],
            Stoich_GI_Network_DF$GI_scores[(Stoich_GI_Network_DF$Region == 4) & (Stoich_GI_Network_DF$GI_scores < 0)])
# p-value = 0.6186
