# Xavier Castellanos-Girouard

# Date First Created: June 5 2024
# Date Last Modified: June 5 2024


#### Import Libraries ####

library(dplyr)
library(tidyr)
library(ggplot2)


#### Import Data ####

dir <- "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/"

ModuleOverlap_Control_DF <- read.csv(paste0(dir, "Python_Module_Overlap/results/GI_PPI_ModuleOverlap_Controls.csv"))
colnames(ModuleOverlap_Control_DF)[1] <- "Cluster_ID"

ModuleOverlap_DF <- read.csv(paste0(dir, "Python_Module_Overlap/results/GI_PPI_optimal_module_overlap.csv"))

#### Mean Jaccard Score Analysis ####

mean_Jaccard_scores <- c()
curr_row = 112
for (i in seq(1, nrow(ModuleOverlap_Control_DF), 112)){
  mean_Jaccard <- mean(ModuleOverlap_Control_DF$Jaccard_index[i:curr_row])
  mean_Jaccard_scores <- c(mean_Jaccard_scores, mean_Jaccard)
  curr_row = curr_row + 112
  }

ggplot(data=NULL,
       mapping = aes(x = mean_Jaccard_scores)) +
  geom_histogram(fill = "black",
                 color = "white",
                 bins = 100) +
  geom_vline(xintercept = mean(ModuleOverlap_DF$Jaccard_index),
             color = "red") +
  theme_classic() +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0, 0.275)) + 
  scale_y_continuous(expand = c(0,0))

hist(mean_Jaccard_scores, nclass = 50)

View(ModuleOverlap_Control_DF)
