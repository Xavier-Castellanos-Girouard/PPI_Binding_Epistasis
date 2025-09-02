# Xavier Castellanos-Girouard

# Date First Created: October 16 2024
# Date Last Modified: October 16 2024

#### Import libraries ####

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)

#### Import data ####

Costanzo_ExE_DF <- read.csv("../data/GI_data/SGA_ExE.txt", sep = "\t")
Costanzo_NxN_DF <- read.csv("../data/GI_data/SGA_NxN.txt", sep = "\t")

Human_Kd_GI_DF <- read.csv("../results/Kd_and_GI/Human_Kd_GI.csv", row.names = 1)

# Import Yeast Kd values
Yeast_Kd_GI_Network_DF <- 
  read.csv("../results/Kd_and_GI/Yeast_Kd_GI.csv",
           row.names = 1)
row.names(Yeast_Kd_GI_Network_DF) <- NULL

#### Format Data ####

## Extract Essential Genes
Essential_Query <- Costanzo_ExE_DF$Query.Strain.ID
Essential_Array <- Costanzo_ExE_DF$Array.Strain.ID

Essential_Query <- 
  gsub(pattern = "_.*",
       replacement = "",
       x = Essential_Query)

Essential_Array <- 
  gsub(pattern = "_.*",
       replacement = "",
       x = Essential_Array)

Essential_Genes <- unique(c(Essential_Query, Essential_Array))

## Extract NonEssential Genes
NonEssential_Query <- Costanzo_NxN_DF$Query.Strain.ID
NonEssential_Array <- Costanzo_NxN_DF$Array.Strain.ID

NonEssential_Query <- 
  gsub(pattern = "_.*",
       replacement = "",
       x = NonEssential_Query)

NonEssential_Array <- 
  gsub(pattern = "_.*",
       replacement = "",
       x = NonEssential_Array)

NonEssential_Genes <- unique(c(NonEssential_Query, NonEssential_Array))

#### Separate GIs according to Essentiality ####

Yeast_Kd_GI_E_Network_DF <-
  Yeast_Kd_GI_Network_DF %>%
  dplyr::filter((source %in% Essential_Genes) & (target %in% Essential_Genes))

Yeast_Kd_GI_NE_Network_DF <-
  Yeast_Kd_GI_Network_DF %>%
  dplyr::filter((source %in% NonEssential_Genes) & (target %in% NonEssential_Genes))

Yeast_Kd_GI_X_Network_DF <-
  Yeast_Kd_GI_Network_DF %>%
  dplyr::filter(xor((source %in% NonEssential_Genes), (target %in% NonEssential_Genes)))

#### Separate Essential data into Quantiles according to Kd (Yeast) ####

### Scatterplots

## Negative GI
Yeast_Kd_NegGI_E_Network_DF <- 
  as.data.frame(Yeast_Kd_GI_E_Network_DF) %>%
  dplyr::filter(scores < 0)

# Divide into deciles according to Interaction Stoichiometry
Yeast_Kd_NegGI_E_Network_DF$Division <- 
  dplyr::ntile(x = Yeast_Kd_NegGI_E_Network_DF$Kd, 
               n = 10)

# Set Division to factor type
Yeast_Kd_NegGI_E_Network_DF$Division <- 
  factor(Yeast_Kd_NegGI_E_Network_DF$Division, 
         levels = as.character(sort(unique(Yeast_Kd_NegGI_E_Network_DF$Division))))

Mean_NegGI_E_KdDeciles <-
  Yeast_Kd_NegGI_E_Network_DF %>%
  group_by(Division) %>%
  summarize(GI_mean = mean(abs(scores)),
            log10_Kd_mean = mean(log10(Kd)))

Mean_NegGI_E_KdDeciles$GI_Score_Type <- "Negative"


## Positive GI
# Divide into deciles according to Kd
Yeast_Kd_PosGI_E_Network_DF <- 
  as.data.frame(Yeast_Kd_GI_E_Network_DF) %>%
  dplyr::filter(scores > 0)

# Divide into deciles according to Kd
Yeast_Kd_PosGI_E_Network_DF$Division <- 
  dplyr::ntile(x = Yeast_Kd_PosGI_E_Network_DF$Kd, 
               n = 10)

# Set Division to factor type
Yeast_Kd_PosGI_E_Network_DF$Division <- 
  factor(Yeast_Kd_PosGI_E_Network_DF$Division, 
         levels = as.character(sort(unique(Yeast_Kd_PosGI_E_Network_DF$Division))))

Mean_PosGI_E_KdDeciles <-
  Yeast_Kd_PosGI_E_Network_DF %>%
  group_by(Division) %>%
  summarize(GI_mean = mean(abs(scores)),
            log10_Kd_mean = mean(log10(Kd)))

Mean_PosGI_E_KdDeciles$GI_Score_Type <- "Positive"

### Model relationship between GI and Kd

Mean_GI_E_KdDeciles <- rbind(Mean_PosGI_E_KdDeciles, Mean_NegGI_E_KdDeciles)

Mean_GI_E_KdDeciles$Kd_mean <- 10^Mean_GI_E_KdDeciles$log10_Kd_mean


#### Yeast Essential GI Kd Lineplot ####

## Linear equation a*x+b
Yeast_NegGI_E_Model <- nls(GI_mean ~ a*log10_Kd_mean+b, data = Mean_GI_E_KdDeciles[Mean_GI_E_KdDeciles$GI_Score_Type == "Negative",], start = list(b = 1, a = 1))
Yeast_PosGI_E_Model <- nls(GI_mean ~ a*log10_Kd_mean+b, data = Mean_GI_E_KdDeciles[Mean_GI_E_KdDeciles$GI_Score_Type == "Positive",], start = list(b = 1, a = 1))

## Get parameters for negative GI in Yeast
Yeast_NegGI_E_Model_sum <- summary(Yeast_NegGI_E_Model)
Y_NegGI_E_a <- Yeast_NegGI_E_Model_sum$parameters['a','Estimate']
Y_NegGI_E_b <- Yeast_NegGI_E_Model_sum$parameters['b','Estimate']

## Get parameters for positive GI in Yeast
Yeast_PosGI_E_Model_sum <- summary(Yeast_PosGI_E_Model)
Y_PosGI_E_a <- Yeast_PosGI_E_Model_sum$parameters['a','Estimate']
Y_PosGI_E_b <- Yeast_PosGI_E_Model_sum$parameters['b','Estimate']


Yeast_E_GI_Kd_lin_plot <-
  ggplot(data = Mean_GI_E_KdDeciles,
         mapping = aes(x = Kd_mean,
                       y = GI_mean,
                       color = GI_Score_Type)) +
  geom_point(size = 2) +
  ggpubr::stat_cor(method="pearson",
                   show.legend = TRUE,
                   cor.coef.name = "R",
                   label.x = -6.5) +
  geom_function(data = Mean_GI_E_KdDeciles[Mean_GI_E_KdDeciles$GI_Score_Type == "Negative",],
                fun = function(x) Y_NegGI_E_a*log10(x) + Y_NegGI_E_b,
                linewidth = 1) +
  geom_function(data = Mean_GI_E_KdDeciles[Mean_GI_E_KdDeciles$GI_Score_Type == "Positive",],
                fun = function(x) Y_PosGI_E_a*log10(x) + Y_PosGI_E_b,
                linewidth = 1) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text.x.bottom = element_text(color = "black", size = 12),
    axis.text.y.left = element_text(color = "black", size = 12),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, 0.3, by = 0.05),
                     limits = c(0, 0.32)) +
  scale_x_log10(limits = c(10^-8.25, 10^-3.7),
                breaks = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^"-8"),
                           expression(10^"-7"), expression(10^"-6"),
                           expression(10^"-5"), expression(10^"-4"),
                           expression(10^"-3"))) +
  scale_color_manual(values = c("blue", "#F5E700")) +
  xlab("Kd (M)") +
  ylab("|Mean GI Score|") +
  ggtitle("Essential x Essential")

#Yeast_E_GI_Kd_lin_plot

#### Separate NonEssential data into Quantiles according to Kd (Yeast) ####

### Scatterplots

## Negative GI
Yeast_Kd_NegGI_NE_Network_DF <- 
  as.data.frame(Yeast_Kd_GI_NE_Network_DF) %>%
  dplyr::filter(scores < 0)

# Divide into deciles according to Interaction Stoichiometry
Yeast_Kd_NegGI_NE_Network_DF$Division <- 
  dplyr::ntile(x = Yeast_Kd_NegGI_NE_Network_DF$Kd, 
               n = 10)

# Set Division to factor type
Yeast_Kd_NegGI_NE_Network_DF$Division <- 
  factor(Yeast_Kd_NegGI_NE_Network_DF$Division, 
         levels = as.character(sort(unique(Yeast_Kd_NegGI_NE_Network_DF$Division))))

Mean_NegGI_NE_KdDeciles <-
  Yeast_Kd_NegGI_NE_Network_DF %>%
  group_by(Division) %>%
  summarize(GI_mean = mean(abs(scores)),
            log10_Kd_mean = mean(log10(Kd)))

Mean_NegGI_NE_KdDeciles$GI_Score_Type <- "Negative"


## Positive GI
# Divide into deciles according to Kd
Yeast_Kd_PosGI_NE_Network_DF <- 
  as.data.frame(Yeast_Kd_GI_NE_Network_DF) %>%
  dplyr::filter(scores > 0)

# Divide into deciles according to Kd
Yeast_Kd_PosGI_NE_Network_DF$Division <- 
  dplyr::ntile(x = Yeast_Kd_PosGI_NE_Network_DF$Kd, 
               n = 10)

# Set Division to factor type
Yeast_Kd_PosGI_NE_Network_DF$Division <- 
  factor(Yeast_Kd_PosGI_NE_Network_DF$Division, 
         levels = as.character(sort(unique(Yeast_Kd_PosGI_NE_Network_DF$Division))))

Mean_PosGI_NE_KdDeciles <-
  Yeast_Kd_PosGI_NE_Network_DF %>%
  group_by(Division) %>%
  summarize(GI_mean = mean(abs(scores)),
            log10_Kd_mean = mean(log10(Kd)))

Mean_PosGI_NE_KdDeciles$GI_Score_Type <- "Positive"

### Model relationship between GI and Kd

Mean_GI_NE_KdDeciles <- rbind(Mean_PosGI_NE_KdDeciles, Mean_NegGI_NE_KdDeciles)

Mean_GI_NE_KdDeciles$Kd_mean <- 10^Mean_GI_NE_KdDeciles$log10_Kd_mean

#### Yeast NonEssential GI Kd Lineplot ####

## Linear equation a*x+b
Yeast_NegGI_NE_Model <- nls(GI_mean ~ a*log10_Kd_mean+b, data = Mean_GI_NE_KdDeciles[Mean_GI_NE_KdDeciles$GI_Score_Type == "Negative",], start = list(b = 1, a = 1))
Yeast_PosGI_NE_Model <- nls(GI_mean ~ a*log10_Kd_mean+b, data = Mean_GI_NE_KdDeciles[Mean_GI_NE_KdDeciles$GI_Score_Type == "Positive",], start = list(b = 1, a = 1))

## Get parameters for negative GI in Yeast
Yeast_NegGI_NE_Model_sum <- summary(Yeast_NegGI_NE_Model)
Y_NegGI_NE_a <- Yeast_NegGI_NE_Model_sum$parameters['a','Estimate']
Y_NegGI_NE_b <- Yeast_NegGI_NE_Model_sum$parameters['b','Estimate']

## Get parameters for positive GI in Yeast
Yeast_PosGI_NE_Model_sum <- summary(Yeast_PosGI_NE_Model)
Y_PosGI_NE_a <- Yeast_PosGI_NE_Model_sum$parameters['a','Estimate']
Y_PosGI_NE_b <- Yeast_PosGI_NE_Model_sum$parameters['b','Estimate']


Yeast_NE_GI_Kd_lin_plot <-
  ggplot(data = Mean_GI_NE_KdDeciles,
         mapping = aes(x = Kd_mean,
                       y = GI_mean,
                       color = GI_Score_Type)) +
  geom_point(size = 2) +
  ggpubr::stat_cor(method="pearson",
                   show.legend = TRUE,
                   cor.coef.name = "R",
                   #alternative = "less",
                   label.x = -6.5) +
  geom_function(data = Mean_GI_NE_KdDeciles[Mean_GI_NE_KdDeciles$GI_Score_Type == "Negative",],
                fun = function(x) Y_NegGI_NE_a*log10(x) + Y_NegGI_NE_b,
                linewidth = 1) +
  geom_function(data = Mean_GI_NE_KdDeciles[Mean_GI_NE_KdDeciles$GI_Score_Type == "Positive",],
                fun = function(x) Y_PosGI_NE_a*log10(x) + Y_PosGI_NE_b,
                linewidth = 1) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text.x.bottom = element_text(color = "black", size = 12),
    axis.text.y.left = element_text(color = "black", size = 12),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, 0.3, by = 0.05),
                     limits = c(0, 0.32)) +
  scale_x_log10(limits = c(10^-8.25, 10^-3),
                breaks = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^"-8"),
                           expression(10^"-7"), expression(10^"-6"),
                           expression(10^"-5"), expression(10^"-4"),
                           expression(10^"-3"))) +
  scale_color_manual(values = c("blue", "#F5E700")) +
  xlab("Kd (M)") +
  ylab("|Mean GI Score|") +
  ggtitle("Non-Essential x Non-Essential")

#Yeast_NE_GI_Kd_lin_plot


#### Separate Essential X NonEssential data into Quantiles according to Kd (Yeast) ####

### Scatterplots

## Negative GI
Yeast_Kd_NegGI_X_Network_DF <- 
  as.data.frame(Yeast_Kd_GI_X_Network_DF) %>%
  dplyr::filter(scores < 0)

# Divide into deciles according to Interaction Stoichiometry
Yeast_Kd_NegGI_X_Network_DF$Division <- 
  dplyr::ntile(x = Yeast_Kd_NegGI_X_Network_DF$Kd, 
               n = 10)

# Set Division to factor type
Yeast_Kd_NegGI_X_Network_DF$Division <- 
  factor(Yeast_Kd_NegGI_X_Network_DF$Division, 
         levels = as.character(sort(unique(Yeast_Kd_NegGI_X_Network_DF$Division))))

Mean_NegGI_X_KdDeciles <-
  Yeast_Kd_NegGI_X_Network_DF %>%
  group_by(Division) %>%
  summarize(GI_mean = mean(abs(scores)),
            log10_Kd_mean = mean(log10(Kd)))

Mean_NegGI_X_KdDeciles$GI_Score_Type <- "Negative"


## Positive GI
# Divide into deciles according to Kd
Yeast_Kd_PosGI_X_Network_DF <- 
  as.data.frame(Yeast_Kd_GI_X_Network_DF) %>%
  dplyr::filter(scores > 0)

# Divide into deciles according to Kd
Yeast_Kd_PosGI_X_Network_DF$Division <- 
  dplyr::ntile(x = Yeast_Kd_PosGI_X_Network_DF$Kd, 
               n = 10)

# Set Division to factor type
Yeast_Kd_PosGI_X_Network_DF$Division <- 
  factor(Yeast_Kd_PosGI_X_Network_DF$Division, 
         levels = as.character(sort(unique(Yeast_Kd_PosGI_X_Network_DF$Division))))

Mean_PosGI_X_KdDeciles <-
  Yeast_Kd_PosGI_X_Network_DF %>%
  group_by(Division) %>%
  summarize(GI_mean = mean(abs(scores)),
            log10_Kd_mean = mean(log10(Kd)))

Mean_PosGI_X_KdDeciles$GI_Score_Type <- "Positive"

### Model relationship between GI and Kd

Mean_GI_X_KdDeciles <- rbind(Mean_PosGI_X_KdDeciles, Mean_NegGI_X_KdDeciles)

Mean_GI_X_KdDeciles$Kd_mean <- 10^Mean_GI_X_KdDeciles$log10_Kd_mean


#### Yeast Essential X NonEssential GI Kd Lineplot ####

## Linear equation a*x+b
Yeast_NegGI_X_Model <- nls(GI_mean ~ a*log10_Kd_mean+b, data = Mean_GI_X_KdDeciles[Mean_GI_X_KdDeciles$GI_Score_Type == "Negative",], start = list(b = 1, a = 1))
Yeast_PosGI_X_Model <- nls(GI_mean ~ a*log10_Kd_mean+b, data = Mean_GI_X_KdDeciles[Mean_GI_X_KdDeciles$GI_Score_Type == "Positive",], start = list(b = 1, a = 1))

## Get parameters for negative GI in Yeast
Yeast_NegGI_X_Model_sum <- summary(Yeast_NegGI_X_Model)
Y_NegGI_X_a <- Yeast_NegGI_X_Model_sum$parameters['a','Estimate']
Y_NegGI_X_b <- Yeast_NegGI_X_Model_sum$parameters['b','Estimate']

## Get parameters for positive GI in Yeast
Yeast_PosGI_X_Model_sum <- summary(Yeast_PosGI_X_Model)
Y_PosGI_X_a <- Yeast_PosGI_X_Model_sum$parameters['a','Estimate']
Y_PosGI_X_b <- Yeast_PosGI_X_Model_sum$parameters['b','Estimate']


Yeast_X_GI_Kd_lin_plot <-
  ggplot(data = Mean_GI_X_KdDeciles,
         mapping = aes(x = Kd_mean,
                       y = GI_mean,
                       color = GI_Score_Type)) +
  geom_point(size = 2) +
  ggpubr::stat_cor(method="pearson",
                   show.legend = TRUE,
                   cor.coef.name = "R",
                   label.x = -6.5) +
  geom_function(data = Mean_GI_X_KdDeciles[Mean_GI_X_KdDeciles$GI_Score_Type == "Negative",],
                fun = function(x) Y_NegGI_X_a*log10(x) + Y_NegGI_X_b,
                linewidth = 1) +
  geom_function(data = Mean_GI_X_KdDeciles[Mean_GI_X_KdDeciles$GI_Score_Type == "Positive",],
                fun = function(x) Y_PosGI_X_a*log10(x) + Y_PosGI_X_b,
                linewidth = 1) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text.x.bottom = element_text(color = "black", size = 12),
    axis.text.y.left = element_text(color = "black", size = 12),
    axis.title.x.bottom = element_text(color = "black", size = 14),
    axis.title.y.left = element_text(color = "black", size = 14),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, 0.3, by = 0.05),
                     limits = c(0, 0.32)) +
  scale_x_log10(limits = c(10^-8.25, 10^-2.5),
                breaks = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(expression(10^"-8"),
                           expression(10^"-7"), expression(10^"-6"),
                           expression(10^"-5"), expression(10^"-4"),
                           expression(10^"-3"))) +
  scale_color_manual(values = c("blue", "#F5E700")) +
  xlab("Kd (M)") +
  ylab("|Mean GI Score|") +
  ggtitle("Non-Essential x Essential")

#Yeast_X_GI_Kd_lin_plot


#### Combine Figures ####

# Yeast and Human Combined
CombinedPlot <-
  plot_grid(Yeast_E_GI_Kd_lin_plot,
            Yeast_X_GI_Kd_lin_plot,
            Yeast_NE_GI_Kd_lin_plot,
            labels = c("a", "b", "c"),
            rel_heights = c(1, 1, 1),
            ncol = 3)

png("../results/Figures/Extended_Data_Figures/SupplementaryFigure7.png", units="mm", width=280, height=80, res=300)
CombinedPlot
dev.off()
