# Xavier Castellanos-Girouard

# Date First Created: June 25 2024
# Date Last Modified: June 25 2024


#### Import Libraries ####

library(dplyr)
library(tidyr)
library(ggplot2)

#### Import Data ####

dir <- "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/"

ConnectorOverlap_Control_DF <- read.csv(paste0(dir, "Python_Connector_Overlap/results/GI_PPI_ConnectorOverlap_Controls.csv"))

ConnectorOverlap_DF <- read.csv(paste0(dir, "Python_Connector_Overlap/results/Shortest_Path_Connector_Overlap.csv"))

#### Overlap Analysis ####

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
ggplot(data=NULL,
       mapping = aes(x = number_of_overlapping_connectors)) +
  geom_histogram(fill = "black",
                 color = "white",
                 bins = 100) +
  geom_vline(xintercept = nrow(ConnectorOverlap_DF),
             color = "red") +
  theme_classic() +
  scale_x_continuous(expand = c(0,0),
                     limits = c(50, 750)) + 
  scale_y_continuous(expand = c(0,0))
