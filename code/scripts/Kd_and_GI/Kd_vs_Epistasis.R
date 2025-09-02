# Xavier Castellanos-Girouard
# 
# Determining the Relationship between Kd and GI
#
# Date First Created: November 9 2023
# Date Last Modified: May 22 2024

#### Import libraries ####

library(dplyr)
library(tidyr)
library(data.table)
library(readxl)

#### Import Data ####

# Import Yeast GI network
GI_Network_DF <- read.csv("../results/Stoichiometry_and_GI/costanzo_2016_longer_withoutReps.csv")
row.names(GI_Network_DF) <- NULL

# Import Yeast Kd values
Kd_PPI_Network_DF <- read.csv("../results/Kd_and_GI/Yeast_Estimated_Kd.csv",
                                  row.names = 1)
row.names(Kd_PPI_Network_DF) <- NULL

# Import Human Kd values
Human_Kd_PPI_Network_DF <- read.csv("../results/Kd_and_GI/Human_Estimated_Kd.csv")

#### Merge GI values and Kds (Yeast) ####

Kd_PPI_Network_DT <- setDT(Kd_PPI_Network_DF, 
                           key = c("source",
                                   "target"))

GI_Network_DT <- setDT(GI_Network_DF,
                       key = c("ORF_query",
                               "ORF_array"))

# Get GI values from interactions in Kd network (match A:B with A:B)
Kd_GI_Network_DT_1 <-
  merge.data.table(x = Kd_PPI_Network_DT,
                   y = GI_Network_DT,
                   by.x = c("source",
                            "target"),
                   by.y = c("ORF_query",
                            "ORF_array"),
                   all.x = FALSE,
                   all.y = FALSE)

# Get GI values from interactions in Kd network (match A:B with B:A)
Kd_GI_Network_DT_2 <-
  merge.data.table(x = Kd_PPI_Network_DT,
                   y = GI_Network_DT,
                   by.x = c("source",
                            "target"),
                   by.y = c("ORF_array",
                            "ORF_query"),
                   all.x = FALSE,
                   all.y = FALSE)

## Combine data tables
Kd_GI_Network_DT <- rbind(Kd_GI_Network_DT_1, Kd_GI_Network_DT_2)


# Export
write.csv(Kd_GI_Network_DT, "../results/Kd_and_GI/Yeast_Kd_GI.csv")


#### Separate data into Quantiles according to Kd (Yeast) ####

Kd_GI_Network_DF <- as.data.frame(Kd_GI_Network_DT)

### Scatterplots

## Negative GI
Kd_NegGI_Network_DF <- 
  as.data.frame(Kd_GI_Network_DF) %>%
  dplyr::filter(scores < 0)

# Divide into deciles according to Interaction Stoichiometry
Kd_NegGI_Network_DF$Division <- 
  dplyr::ntile(x = Kd_NegGI_Network_DF$Kd, 
               n = 10)

# Set Division to factor type
Kd_NegGI_Network_DF$Division <- 
  factor(Kd_NegGI_Network_DF$Division, 
         levels = as.character(sort(unique(Kd_NegGI_Network_DF$Division))))

Mean_NegGI_KdDeciles <-
  Kd_NegGI_Network_DF %>%
  group_by(Division) %>%
  summarize(GI_mean = mean(abs(scores)),
            log10_Kd_mean = mean(log10(Kd)))

Mean_NegGI_KdDeciles$GI_Score_Type <- "Negative"


## Positive GI
# Divide into deciles according to Kd
Kd_PosGI_Network_DF <- 
  as.data.frame(Kd_GI_Network_DF) %>%
  dplyr::filter(scores > 0)

# Divide into deciles according to Kd
Kd_PosGI_Network_DF$Division <- 
  dplyr::ntile(x = Kd_PosGI_Network_DF$Kd, 
               n = 10)

# Set Division to factor type
Kd_PosGI_Network_DF$Division <- 
  factor(Kd_PosGI_Network_DF$Division, 
         levels = as.character(sort(unique(Kd_PosGI_Network_DF$Division))))

Mean_PosGI_KdDeciles <-
  Kd_PosGI_Network_DF %>%
  group_by(Division) %>%
  summarize(GI_mean = mean(abs(scores)),
            log10_Kd_mean = mean(log10(Kd)))

Mean_PosGI_KdDeciles$GI_Score_Type <- "Positive"

### Model relationship between GI and Kd

Mean_GI_KdDeciles <- rbind(Mean_PosGI_KdDeciles, Mean_NegGI_KdDeciles)

Mean_GI_KdDeciles$Kd_mean <- 10^Mean_GI_KdDeciles$log10_Kd_mean

#NegGI_Model <- nls(GI_mean ~ (a/(1 + exp(b*(log10_Kd_mean-d)))+c), data = Mean_GI_KdDeciles[Mean_GI_KdDeciles$GI_Score_Type == "Negative",], start = list(b = 3, a = 0.1, c = 0.1, d = -6))
#NegGI_Model <- nls(GI_mean ~ (1/(1 + exp(b*(log10_Kd_mean-d)))), data = Mean_GI_KdDeciles[Mean_GI_KdDeciles$GI_Score_Type == "Negative",], start = list(b = 0.26, d = -12))
#NegGI_Model <- nls(GI_mean ~ (1/(1 + exp(log10_Kd_mean-d))), data = Mean_GI_KdDeciles[Mean_GI_KdDeciles$GI_Score_Type == "Negative",], start = list(d = -12))
PosGI_Model <- nls(GI_mean ~ a*log10_Kd_mean+b, data = Mean_GI_KdDeciles[Mean_GI_KdDeciles$GI_Score_Type == "Positive",], start = list(b = 1, a = 1))

write.csv(Mean_GI_KdDeciles, "../results/Kd_and_GI/Yeast_Mean_GI_KdDeciles.csv")


#### Import and format data from Horlbeck ####

# * Maybe make script in GI network formatting to not have to do this again

# Importing GI data from Horlbeck et al. (Cell 2018) supplementary data S12
# The "gene GI scores and correlations" is the main interest for this analysis
Gene_GI_Score_Correlation <- 
  read_excel("../data/GI_data/GI_Horlbeck.xlsx", 
             sheet = "gene GI scores and correlations",
             col_types = c("text", "text", rep('numeric', 12)))

# Organizing GI scores and GI correlations (replicate averages) as distinct datasets.
K562_GI_Score <- Gene_GI_Score_Correlation[4:nrow(Gene_GI_Score_Correlation), c(1,2,7)]
K652_GI_Correlation <- Gene_GI_Score_Correlation[4:nrow(Gene_GI_Score_Correlation), c(1,2,8)]
Jurkat_GI_Score <- Gene_GI_Score_Correlation[4:nrow(Gene_GI_Score_Correlation), c(1,2,13)]
Jurkat_GI_Correlation <- Gene_GI_Score_Correlation[4:nrow(Gene_GI_Score_Correlation), c(1,2,14)]
# Note: Main focus will be on the GI scores. It is the metric used in the paper


#Renaming columns for clarity
colnames(K562_GI_Score) <- c("Gene_1", "Gene_2", "GI_Score")
colnames(K652_GI_Correlation) <- c("Gene_1", "Gene_2", "GI_Correlation")
colnames(Jurkat_GI_Score) <- c("Gene_1", "Gene_2", "GI_Score")
colnames(Jurkat_GI_Correlation) <- c("Gene_1", "Gene_2", "GI_Correlation")

#### Merge GI values and Kds (Human) ####

#Merge of the two dataframes by gene pair:
Kd_GI_merge_combo1 <- 
  merge(x = Human_Kd_PPI_Network_DF,
        y = K562_GI_Score,
        by.x = c("Bait", "Prey"),
        by.y = c("Gene_1", "Gene_2"),
        all.x = FALSE,
        all.y = FALSE)

length((Kd_GI_merge_combo1$GI_Score[!(is.na(Kd_GI_merge_combo1$GI_Score))]))
# 302


# Second merge combination:
Kd_GI_merge_combo2 <- 
  merge(x = Human_Kd_PPI_Network_DF,
        y = K562_GI_Score, 
        by.x = c("Bait", "Prey"),
        by.y = c("Gene_2", "Gene_1"),
        all.x = FALSE,
        all.y = FALSE)

length((Kd_GI_merge_combo2$GI_Score[!(is.na(Kd_GI_merge_combo2$GI_Score))]))
# 326 common interactions for this combination.

Human_Kd_GI_Network_DF <- rbind(Kd_GI_merge_combo1, Kd_GI_merge_combo2)

Human_Kd_GI_Network_DF <- Human_Kd_GI_Network_DF[!duplicated(Human_Kd_GI_Network_DF),]

Human_Kd_GI_Network_DF <- Human_Kd_GI_Network_DF[!is.na(Human_Kd_GI_Network_DF$GI_Score),]
# Total of 620 Common interactions.

write.csv(Human_Kd_GI_Network_DF, "../results/Kd_and_GI/Human_Kd_GI.csv")

#### Separate data into Quantiles according to Kd (Human) ####

### Negative GI
Human_Kd_GI_Network_DF$Quintile[Human_Kd_GI_Network_DF$GI_Score < 0] <- 
  Human_Kd_GI_Network_DF$Kd[Human_Kd_GI_Network_DF$GI_Score < 0] %>%
  ntile(n = 5)

Human_Kd_GI_Network_DF$Quintile[Human_Kd_GI_Network_DF$GI_Score < 0] <-
  factor(Human_Kd_GI_Network_DF$Quintile[Human_Kd_GI_Network_DF$GI_Score < 0],
         levels = 1:5)

negquintile_mean <-
  Human_Kd_GI_Network_DF[Human_Kd_GI_Network_DF$GI_Score < 0,] %>%
  dplyr::group_by(Quintile) %>%
  summarize(MeanGI = mean(abs(GI_Score)),
            log10_Kd_mean = mean(log10(Kd)))

### Positive GI
Human_Kd_GI_Network_DF$Quintile[Human_Kd_GI_Network_DF$GI_Score > 0] <- 
  Human_Kd_GI_Network_DF$Kd[Human_Kd_GI_Network_DF$GI_Score > 0] %>%
  ntile(n = 5)

Human_Kd_GI_Network_DF$Quintile[Human_Kd_GI_Network_DF$GI_Score > 0] <-
  factor(Human_Kd_GI_Network_DF$Quintile[Human_Kd_GI_Network_DF$GI_Score > 0],
         levels = 1:5)


posquintile_mean <-
  Human_Kd_GI_Network_DF[Human_Kd_GI_Network_DF$GI_Score > 0,] %>%
  group_by(Quintile) %>%
  summarize(MeanGI = mean(abs(GI_Score)),
            log10_Kd_mean = mean(log10(Kd)))


## Combine both
negquintile_mean$GI_Type <- "Negative"
posquintile_mean$GI_Type <- "Positive" 

GIQuintile_mean <- rbind(negquintile_mean, posquintile_mean)

GIQuintile_mean$Kd_mean <- 10^GIQuintile_mean$log10_Kd_mean


Human_NegGI_Model <- nls(MeanGI ~ (a/(1 + exp(b*(log10_Kd_mean-d)))+c), data = GIQuintile_mean[GIQuintile_mean$GI_Type == "Negative",], start = list(b = 5.2, a = 0.8, c = 1.1, d = -4.8))
Human_NegGI_Model <- nls(MeanGI ~ a*log10_Kd_mean+b, data = GIQuintile_mean[GIQuintile_mean$GI_Type == "Negative",], start = list(b = 1, a = 1))
Human_PosGI_Model <- nls(MeanGI ~ a*log10_Kd_mean+b, data = GIQuintile_mean[GIQuintile_mean$GI_Type == "Positive",], start = list(b = 1, a = 1))

write.csv(GIQuintile_mean, "../results/Kd_and_GI/Human_Mean_GI_KdQuintiles.csv")
