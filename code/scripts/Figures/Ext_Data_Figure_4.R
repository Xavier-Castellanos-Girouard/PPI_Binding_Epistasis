# Xavier Castellanos-Girouard
#
#
# Check functional enrichment of specific pathway types in stoichiometry map
#
# Date First Created: February 15 2024
# Date Last Modified: August 12 2025

#### Import Libraries ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(readxl)
library(stringr)
library(org.Sc.sgd.db)

#### Graph parameters ####

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

#### Import data ####

#ComplexPortal <-
#  read.csv(paste0("../data/ComplexPortal_559292.tsv"),
#           sep = "\t")

Yeast_PPI_Stoich_DF <-
  read.csv("../results/Stoichiometry_and_GI/Yeast_Interactome_Stoich.csv", row.names = 1)

Yeast_GI_Stoich_DF <- read.csv("../results/Stoichiometry_and_GI/Yeast_Stoich_GI.csv", row.names = 1)

Yeast_PPI_Kd_DF <- read.csv("../results/Kd_and_GI/Yeast_Estimated_Kd.csv", row.names = 1)
Yeast_Kd_GI_DF <- read.csv("../results/Kd_and_GI/results/Yeast_Kd_GI.csv", row.names = 1)

## Import kinase dataset from BioGRID (yeast kinome project)

# Import yeast kinase list from BioGRID (yeast kinome project)
Yeast_Kinase_list_DF <- read.csv("../data/Miscellaneous/BIOGRID-PROJECT-kinome_project_sc-GENES-4.4.247.projectindex.txt", sep = "\t")

## Import Yeast Metabolism from BioCyc
Metabolic_Network <- read.csv("../data/Miscellaneous/Copy-of-All-pathways-of-S.-cerevisiae-S288c.txt", sep = "\t")
Metabolic_Network <- 
  Metabolic_Network %>%
  dplyr::select(c("Pathways", "Common.Name", "Accession.1")) %>%
  tidyr::separate_longer_delim(col = "Accession.1",
                               delim = " // ")

#Metabolic_Network <- Metabolic_Network[grepl(pattern = "glycolysis", Metabolic_Network$Common.Name),]

ribosome_large_subunit <- c("YPL220W", "YGL135W", "YFR031C-A", "YIL018W", 
                            "YOR063W", "YBR031W", "YDR012W", "YPL131W", 
                            "YML073C", "YLR448W", "YGL076C", "YPL198W", 
                            "YHL033C", "YLL045C", "YGL147C", "YNL067W", 
                            "YLR075W", "YPR102C", "YGR085C", "YEL054C", 
                            "YDR418W", "YDL082W", "YMR142C", "YKL006W",
                            "YHL001W", "YLR029C", "YMR121C", "YIL133C",
                            "YNL069C", "YKL180W", "YJL177W", "YOL120C",
                            "YNL301C", "YBR084C-A", "YBL027W", "YMR242C", 
                            "YOR312C", "YBR191W", "YPL079W", "YLR061W",
                            "YFL034C-A", "YBL087C", "RPL23B", "YGL031C",
                            "YGR148C", "YOL127W", "YLR344W", "YGR034W",
                            "YHR010W", "YDR471W", "YGL103W", "YFR032C-A",
                            "YGL030W", "YDL075W", "YLR406C", "YBL092W",
                            "YPL143W", "YOR234C", "YER056C-A", "YIL052C",
                            "YDL191W", "YDL136W", "YMR194W", "YPL249C-A",
                            "YLR185W", "YDR500C", "YLR325C", "YJL189W",
                            "YIL148W", "YKR094C", "YDL184C", "YDL133C-A",
                            "YNL162W", "YHR141C", "YPR043W", "YJR094W-A")

proteasome_26S <- c("YFR004W", "YDL007W", "YDR394W", "YKL145W", "YGL048C",
                    "YOR117W", "YMR314W", "YML092C", "YOR362C", "YGR135W",
                    "YGR253C", "YGL011C", "YOL038W", "YBL041W", "YER012W",
                    "YER094C", "YFR050C", "YPR103W", "YJL001W", "YOR157C",
                    "YHR027C", "YIL075C", "YER021W", "YDL147W", "YPR108W",
                    "YDR427W", "YHR200W", "YFR052W", "YOR259C", "YDL097C",
                    "YDR363W-A", "YLR421C", "YOR261C")

#### Format and merge Kinome Data (OUTDATED) ####

## Yeast Kinome
Yeast_Kinome_DF <-
  Yeast_Kinome_DF %>%
  filter(Category %in% c("Kinase Catalytic", "Phosphatase Catalytic"))

Yeast_Kinome_DF$Interactors <- str_remove_all(Yeast_Kinome_DF$Interactors, pattern = "(\\[|\\])")

Yeast_Kinome_DF <-
  Yeast_Kinome_DF %>%
  separate_longer_delim(cols = "Interactors",
                        delim = ", ")

Yeast_Kinome_DF$Interactor_ORF <- NA

map_Common2ORF <- as.list(org.Sc.sgdCOMMON2ORF)
# Remove probes that do not map in COMMON2ORF
map_Common2ORF <- map_Common2ORF[!is.na(map_Common2ORF)]

i = 0
for (gene in Yeast_Kinome_DF$Interactors){
  print(gene)
  i = i + 1
  ORF <- map_Common2ORF[gene]
  
  if (is.null(unname(unlist(ORF)))){ # If ORF is null, then interactor gene name was the ORF ID
    Yeast_Kinome_DF$Interactor_ORF[i] <- Yeast_Kinome_DF$Interactors[i]
  } else{
    Yeast_Kinome_DF$Interactor_ORF[i] <- map_Common2ORF[gene]
  }
}

## Remove instances where multiple ORFs can be mapped
multiple_match <- 
  lapply(Yeast_Kinome_DF$Interactor_ORF, 
         FUN = function(x){
           if (length(x) > 1){
             return(FALSE)
           } else {return(TRUE)}
         })

Yeast_Kinome_DF <- Yeast_Kinome_DF[unname(unlist(multiple_match)),]

Yeast_Kinome_Kd_comb1_DF <- 
  merge(x = Yeast_Kinome_DF,
        y = Yeast_PPI_Stoich_DF,
        by.x = c("Orf Name", "Interactor_ORF"),
        by.y = c("source", "target"),
        all.x = FALSE,
        all.y = FALSE)

Yeast_Kinome_Kd_comb2_DF <- 
  merge(x = Yeast_Kinome_DF,
        y = Yeast_PPI_Stoich_DF,
        by.x = c("Orf Name", "Interactor_ORF"),
        by.y = c("target", "source"),
        all.x = FALSE,
        all.y = FALSE)

Yeast_Kinome_Kd_DF <- rbind(Yeast_Kinome_Kd_comb1_DF, Yeast_Kinome_Kd_comb2_DF)

wilcox.test(x = Yeast_Kinome_Kd_DF$Interaction_Stoichiometry, 
            y = Yeast_PPI_Stoich_DF$Interaction_Stoichiometry)

#### Yeast Kinase interaction stoichiometry and GI ####
Yeast_Kinase_list_DF <- 
  Yeast_Kinase_list_DF %>%
  filter(CATEGORY.VALUES=="Kinase")

## Stoichiometry
Yeast_PPI_Stoich_DF$source_isKinase <- Yeast_PPI_Stoich_DF$source %in% Yeast_Kinase_list_DF$SYSTEMATIC.NAME
Yeast_PPI_Stoich_DF$target_isKinase <- Yeast_PPI_Stoich_DF$target %in% Yeast_Kinase_list_DF$SYSTEMATIC.NAME

Yeast_PPI_Stoich_DF$Kinase_Interaction <- Yeast_PPI_Stoich_DF$source_isKinase | Yeast_PPI_Stoich_DF$target_isKinase

Yeast_PPI_Stoich_DF$Kinase_Interaction <- 
  factor(Yeast_PPI_Stoich_DF$Kinase_Interaction,
         levels = c(TRUE, FALSE))

Yeast_PPI_Stoich_DF <- Yeast_PPI_Stoich_DF[order(Yeast_PPI_Stoich_DF$Kinase_Interaction, decreasing = TRUE),]

## GI
Yeast_GI_Stoich_DF$source_isKinase <- Yeast_GI_Stoich_DF$source %in% Yeast_Kinase_list_DF$SYSTEMATIC.NAME
Yeast_GI_Stoich_DF$target_isKinase <- Yeast_GI_Stoich_DF$target %in% Yeast_Kinase_list_DF$SYSTEMATIC.NAME

Yeast_GI_Stoich_DF$Kinase_Interaction <- Yeast_GI_Stoich_DF$source_isKinase | Yeast_GI_Stoich_DF$target_isKinase

Yeast_GI_Stoich_DF$Kinase_Interaction <- 
  factor(Yeast_GI_Stoich_DF$Kinase_Interaction,
         levels = c(TRUE, FALSE))

Yeast_GI_Stoich_DF <- Yeast_GI_Stoich_DF[order(Yeast_GI_Stoich_DF$Kinase_Interaction, decreasing = TRUE),]

#### Kinase Graphs ####

## Stoichiometry
#kinase_stoich_map_p <-
#  ggplot(data = Yeast_PPI_Stoich_DF,
#         mapping = aes(x = Interaction_Stoichiometry,
#                       y = Abundance_Stoichiometry,
#                       colour = Kinase_Interaction)) +
#  geom_point() +
#  theme(panel.background = element_rect(fill = "white"),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        axis.line.x = element_line(color = "black", linewidth = 0.5),
#        axis.line.y = element_line(color = "black", linewidth = 0.5),
#        axis.ticks = element_line(color = "black"),
#        axis.text.x.bottom = element_text(color = "black", size = 12),
#        axis.text.y.left = element_text(color = "black", size = 12),
#        axis.title.x.bottom = element_text(color = "black", size = 14),
#        axis.title.y.left = element_text(color = "black", size = 14),
#        plot.title = element_text(hjust = 0.5,
#                                  color = "black"),
#        panel.border = element_rect(colour = "black",
#                                    fill = NA,
#                                    linewidth = 1),
#        legend.text=element_text(size=12)) +
#  scale_x_log10() +
#  scale_y_log10() +
#  scale_colour_manual(values = c("red", "black")) +
#  xlab("Abundance Stoichiometry") +
#  ylab("Interaction Stoichiometry") +
#  labs(colour = "Protein-protein interaction\n        involves kinase?")

#kinase_stoich_map_p


kinase_stoich_map_p <-
  ggplot() +
  geom_density_2d(data = Yeast_PPI_Stoich_DF,
                  mapping = aes(x = Interaction_Stoichiometry,
                                y = Abundance_Stoichiometry,
                                colour = Kinase_Interaction),
                  color = "black") +
  geom_density_2d(data = Yeast_PPI_Stoich_DF[Yeast_PPI_Stoich_DF$Kinase_Interaction==TRUE,],
                  mapping = aes(x = Interaction_Stoichiometry,
                                y = Abundance_Stoichiometry,
                                colour = Kinase_Interaction),
                  color = "red",
                  bins = 5) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color = "black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text = element_text(colour = "black", 
                                   size = Axis_text_size),
        legend.title = element_text(colour = "black",
                                    size = Axis_title_size),
        legend.position = "none") +
  scale_x_log10(limits = c(10^-5.3, 10^2),
                breaks = c(10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-5"),
                           expression(10^"-4"), expression(10^"-3"),
                           expression(10^"-2"), expression(10^"-1"),
                           expression(10^0), expression(10^1),
                           expression(10^2))) +
  scale_y_log10(limits = c(10^-2, 10^1.1),
                breaks = c(10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-2"), expression(10^"-1"), 
                           expression(10^"-0"), expression(10^1),
                           expression(10^2))) +
  xlab("Abundance Stoichiometry") +
  ylab("Interaction Stoichiometry") +
  labs(colour = "Protein-protein interaction\n        involves kinase?")

kinase_stoich_map_p

kinase_IntStoich_violin_p <- 
  ggplot(data = Yeast_PPI_Stoich_DF,
         mapping = aes(x = Kinase_Interaction,
                       y = Interaction_Stoichiometry)) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "grey") +
  geom_signif(comparisons = list(c(1,2)), # Only keep significant
              #position="dodge",
              #tip_length = 0.01,
              #step_increase = 0.05,
              #vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text=element_text(size=12)) +
  scale_y_log10(limits = c(10^-6, 10^6),
                breaks = c(10^-6, 10^-4, 10^-2, 10^0, 10^2, 10^4),
                labels = c(expression(10^"-6"), expression(10^"-4"), 
                           expression(10^"-2"), expression(10^0),
                           expression(10^2), expression(10^4))) +
  xlab("PPI Involves Kinase?") +
  ylab("Interaction Stoichiometry")

kinase_IntStoich_violin_p

wilcox.test(x = Yeast_PPI_Stoich_DF$Interaction_Stoichiometry[Yeast_PPI_Stoich_DF$Kinase_Interaction==FALSE], 
            y = Yeast_PPI_Stoich_DF$Interaction_Stoichiometry[Yeast_PPI_Stoich_DF$Kinase_Interaction==TRUE])

## GI
Yeast_negGI_Stoich_DF <- Yeast_GI_Stoich_DF[Yeast_GI_Stoich_DF$GI_scores<0,]
kinase_negGI_violin_p <- 
  ggplot(data = Yeast_negGI_Stoich_DF,
         mapping = aes(x = Kinase_Interaction,
                       y = abs(GI_scores))) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "blue") +
  geom_signif(comparisons = list(c(1,2)), # Only keep significant
              #position="dodge",
              #tip_length = 0.01,
              #step_increase = 0.05,
              #vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text=element_text(size=12)) +
  scale_y_log10(limits = c(10^-4, 10^1),
                breaks = c(10^-4, 10^-2, 10^0),
                labels = c(expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0))) +
  xlab("PPI Involves Kinase?") +
  ylab("Negative GI Score")

kinase_negGI_violin_p

Yeast_posGI_Stoich_DF <- Yeast_GI_Stoich_DF[Yeast_GI_Stoich_DF$GI_scores>0,]
kinase_posGI_violin_p <- 
  ggplot(data = Yeast_posGI_Stoich_DF,
         mapping = aes(x = Kinase_Interaction,
                       y = abs(GI_scores))) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "yellow2") +
  geom_signif(comparisons = list(c(1,2)), # Only keep significant
              #position="dodge",
              #tip_length = 0.01,
              #step_increase = 0.05,
              #vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text=element_text(size=12)) +
  scale_y_log10(limits = c(10^-4, 10^1),
                breaks = c(10^-4, 10^-2, 10^0),
                labels = c(expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0))) +
  xlab("PPI Involves Kinase?") +
  ylab("Positive GI Score")

kinase_posGI_violin_p

#### Yeast Metabolism interaction stoichiometry and GI ####

## Stoichiometry
Yeast_PPI_Stoich_DF$source_isMetabolism <- Yeast_PPI_Stoich_DF$source %in% Metabolic_Network$Accession.1
Yeast_PPI_Stoich_DF$target_isMetabolism <- Yeast_PPI_Stoich_DF$target %in% Metabolic_Network$Accession.1

Yeast_PPI_Stoich_DF$Metabolism_Interaction <- Yeast_PPI_Stoich_DF$source_isMetabolism | Yeast_PPI_Stoich_DF$target_isMetabolism

Yeast_PPI_Stoich_DF$Metabolism_Interaction <- 
  factor(Yeast_PPI_Stoich_DF$Metabolism_Interaction,
         levels = c(TRUE, FALSE))

Yeast_PPI_Stoich_DF <- Yeast_PPI_Stoich_DF[order(Yeast_PPI_Stoich_DF$Metabolism_Interaction, decreasing = TRUE),]

## GI
Yeast_GI_Stoich_DF$source_isMetabolism <- Yeast_GI_Stoich_DF$source %in% Metabolic_Network$Accession.1
Yeast_GI_Stoich_DF$target_isMetabolism <- Yeast_GI_Stoich_DF$target %in% Metabolic_Network$Accession.1

Yeast_GI_Stoich_DF$Metabolism_Interaction <- Yeast_GI_Stoich_DF$source_isMetabolism | Yeast_GI_Stoich_DF$target_isMetabolism

Yeast_GI_Stoich_DF$Metabolism_Interaction <- 
  factor(Yeast_GI_Stoich_DF$Metabolism_Interaction,
         levels = c(TRUE, FALSE))

Yeast_GI_Stoich_DF <- Yeast_GI_Stoich_DF[order(Yeast_GI_Stoich_DF$Metabolism_Interaction, decreasing = TRUE),]

#### Metabolism Graphs ####

enzyme_stoich_map_p <-
  ggplot() +
  geom_density_2d(data = Yeast_PPI_Stoich_DF,
                  mapping = aes(x = Interaction_Stoichiometry,
                                y = Abundance_Stoichiometry,
                                colour = Metabolism_Interaction),
                  color = "black") +
  geom_density_2d(data = Yeast_PPI_Stoich_DF[Yeast_PPI_Stoich_DF$Metabolism_Interaction==TRUE,],
                  mapping = aes(x = Interaction_Stoichiometry,
                                y = Abundance_Stoichiometry,
                                colour = Metabolism_Interaction),
                  color = "red",
                  bins = 5) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color = "black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text = element_text(colour = "black", 
                                   size = Axis_text_size),
        legend.title = element_text(colour = "black",
                                    size = Axis_title_size),
        legend.position = "none") +
  scale_x_log10(limits = c(10^-5.3, 10^2),
                breaks = c(10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-5"),
                           expression(10^"-4"), expression(10^"-3"),
                           expression(10^"-2"), expression(10^"-1"),
                           expression(10^0), expression(10^1),
                           expression(10^2))) +
  scale_y_log10(limits = c(10^-2, 10^1.1),
                breaks = c(10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-2"), expression(10^"-1"), 
                           expression(10^"-0"), expression(10^1),
                           expression(10^2))) +
  scale_colour_manual(values = c("red", "black")) +
  xlab("Abundance Stoichiometry") +
  ylab("Interaction Stoichiometry") +
  labs(colour = "Interaction between\n     two enzymes?")

enzyme_stoich_map_p

enzyme_IntStoich_violin_p <- 
  ggplot(data = Yeast_PPI_Stoich_DF,
         mapping = aes(x = Metabolism_Interaction,
                       y = Interaction_Stoichiometry)) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "grey") +
  geom_signif(comparisons = list(c(1,2)), # Only keep significant
              #position="dodge",
              #tip_length = 0.01,
              #step_increase = 0.05,
              #vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text=element_text(size=12)) +
  scale_y_log10(limits = c(10^-6, 10^6),
                breaks = c(10^-6, 10^-4, 10^-2, 10^0, 10^2, 10^4),
                labels = c(expression(10^"-6"), expression(10^"-4"), 
                           expression(10^"-2"), expression(10^0),
                           expression(10^2), expression(10^4))) +
  xlab("PPI Involves Enzyme?") +
  ylab("Interaction Stoichiometry")

enzyme_IntStoich_violin_p

wilcox.test(x = Yeast_PPI_Stoich_DF$Interaction_Stoichiometry[Yeast_PPI_Stoich_DF$Metabolism_Interaction==FALSE], 
            y = Yeast_PPI_Stoich_DF$Interaction_Stoichiometry[Yeast_PPI_Stoich_DF$Metabolism_Interaction==TRUE])

## GI
Yeast_negGI_Stoich_DF <- Yeast_GI_Stoich_DF[Yeast_GI_Stoich_DF$GI_scores<0,]
enzyme_negGI_violin_p <- 
  ggplot(data = Yeast_negGI_Stoich_DF,
         mapping = aes(x = Metabolism_Interaction,
                       y = abs(GI_scores))) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "blue") +
  geom_signif(comparisons = list(c(1,2)), # Only keep significant
              #position="dodge",
              #tip_length = 0.01,
              #step_increase = 0.05,
              #vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text=element_text(size=12)) +
  scale_y_log10(limits = c(10^-4, 10^1),
                breaks = c(10^-4, 10^-2, 10^0),
                labels = c(expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0))) +
  xlab("PPI Involves Enzyme?") +
  ylab("Negative GI Score")

enzyme_negGI_violin_p

Yeast_posGI_Stoich_DF <- Yeast_GI_Stoich_DF[Yeast_GI_Stoich_DF$GI_scores>0,]
enzyme_posGI_violin_p <- 
  ggplot(data = Yeast_posGI_Stoich_DF,
         mapping = aes(x = Metabolism_Interaction,
                       y = abs(GI_scores))) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "yellow2") +
  geom_signif(comparisons = list(c(1,2)), # Only keep significant
              #position="dodge",
              #tip_length = 0.01,
              #step_increase = 0.05,
              #vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text=element_text(size=12)) +
  scale_y_log10(limits = c(10^-4, 10^1),
                breaks = c(10^-4, 10^-2, 10^0),
                labels = c(expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0))) +
  xlab("PPI Involves Enzyme?") +
  ylab("Positive GI Score")

enzyme_posGI_violin_p

#### Yeast Large Ribosomal Subunit interaction stoichiometry and GI ####

## Stoichiometry
Yeast_PPI_Stoich_DF$source_isRibosome <- Yeast_PPI_Stoich_DF$source %in% ribosome_large_subunit
Yeast_PPI_Stoich_DF$target_isRibosome <- Yeast_PPI_Stoich_DF$target %in% ribosome_large_subunit

Yeast_PPI_Stoich_DF$Ribosome_Interaction <- Yeast_PPI_Stoich_DF$source_isRibosome & Yeast_PPI_Stoich_DF$target_isRibosome

Yeast_PPI_Stoich_DF$Ribosome_Interaction <- 
  factor(Yeast_PPI_Stoich_DF$Ribosome_Interaction,
         levels = c(TRUE, FALSE))

Yeast_PPI_Stoich_DF <- Yeast_PPI_Stoich_DF[order(Yeast_PPI_Stoich_DF$Ribosome_Interaction, decreasing = TRUE),]

## GI
Yeast_GI_Stoich_DF$source_isRibosome <- Yeast_GI_Stoich_DF$source %in% ribosome_large_subunit
Yeast_GI_Stoich_DF$target_isRibosome <- Yeast_GI_Stoich_DF$target %in% ribosome_large_subunit

Yeast_GI_Stoich_DF$Ribosome_Interaction <- Yeast_GI_Stoich_DF$source_isRibosome & Yeast_GI_Stoich_DF$target_isRibosome

Yeast_GI_Stoich_DF$Ribosome_Interaction <- 
  factor(Yeast_GI_Stoich_DF$Ribosome_Interaction,
         levels = c(TRUE, FALSE))

Yeast_GI_Stoich_DF <- Yeast_GI_Stoich_DF[order(Yeast_GI_Stoich_DF$Ribosome_Interaction, decreasing = TRUE),]

#### Ribosome Graphs ####

Ribosome_stoich_map_p <-
  ggplot() +
  geom_density_2d(data = Yeast_PPI_Stoich_DF,
                  mapping = aes(x = Interaction_Stoichiometry,
                                y = Abundance_Stoichiometry,
                                colour = Ribosome_Interaction),
                  color = "black") +
  geom_density_2d(data = Yeast_PPI_Stoich_DF[Yeast_PPI_Stoich_DF$Ribosome_Interaction==TRUE,],
                  mapping = aes(x = Interaction_Stoichiometry,
                                y = Abundance_Stoichiometry,
                                colour = Ribosome_Interaction),
                  color = "red",
                  bins = 5) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color = "black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text = element_text(colour = "black", 
                                   size = Axis_text_size),
        legend.title = element_text(colour = "black",
                                    size = Axis_title_size),
        legend.position = "none") +
  scale_x_log10(limits = c(10^-5.3, 10^2),
                breaks = c(10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-5"),
                           expression(10^"-4"), expression(10^"-3"),
                           expression(10^"-2"), expression(10^"-1"),
                           expression(10^0), expression(10^1),
                           expression(10^2))) +
  scale_y_log10(limits = c(10^-2, 10^1.1),
                breaks = c(10^-2, 10^-1, 10^0, 10^1, 10^2),
                labels = c(expression(10^"-2"), expression(10^"-1"), 
                           expression(10^"-0"), expression(10^1),
                           expression(10^2))) +
  xlab("Abundance Stoichiometry") +
  ylab("Interaction Stoichiometry") +
  labs(colour = "Interaction within Large\n  Ribosomal Subunit?")

Ribosome_stoich_map_p

Ribosome_IntStoich_violin_p <- 
  ggplot(data = Yeast_PPI_Stoich_DF,
         mapping = aes(x = Ribosome_Interaction,
                       y = Interaction_Stoichiometry)) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "grey") +
  geom_signif(comparisons = list(c(1,2)), # Only keep significant
              #position="dodge",
              #tip_length = 0.01,
              #step_increase = 0.05,
              #vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text=element_text(size=12)) +
  scale_y_log10(limits = c(10^-6, 10^6),
                breaks = c(10^-6, 10^-4, 10^-2, 10^0, 10^2, 10^4),
                labels = c(expression(10^"-6"), expression(10^"-4"), 
                           expression(10^"-2"), expression(10^0),
                           expression(10^2), expression(10^4))) +
  xlab("Interaction within Large\n  Ribosomal Subunit?") +
  ylab("Interaction Stoichiometry")

Ribosome_IntStoich_violin_p

wilcox.test(x = Yeast_PPI_Stoich_DF$Interaction_Stoichiometry[Yeast_PPI_Stoich_DF$Metabolism_Interaction==FALSE], 
            y = Yeast_PPI_Stoich_DF$Interaction_Stoichiometry[Yeast_PPI_Stoich_DF$Metabolism_Interaction==TRUE])

## GI
Yeast_negGI_Stoich_DF <- Yeast_GI_Stoich_DF[Yeast_GI_Stoich_DF$GI_scores<0,]
Ribosome_negGI_violin_p <- 
  ggplot(data = Yeast_negGI_Stoich_DF,
         mapping = aes(x = Ribosome_Interaction,
                       y = abs(GI_scores))) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "blue") +
  geom_signif(comparisons = list(c(1,2)), # Only keep significant
              #position="dodge",
              #tip_length = 0.01,
              #step_increase = 0.05,
              #vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text=element_text(size=12)) +
  scale_y_log10(limits = c(10^-4, 10^1),
                breaks = c(10^-4, 10^-2, 10^0),
                labels = c(expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0))) +
  xlab("Interaction within Large\n  Ribosomal Subunit?") +
  ylab("Negative GI Score")

Ribosome_negGI_violin_p

Yeast_posGI_Stoich_DF <- Yeast_GI_Stoich_DF[Yeast_GI_Stoich_DF$GI_scores>0,]
Ribosome_posGI_violin_p <- 
  ggplot(data = Yeast_posGI_Stoich_DF,
         mapping = aes(x = Ribosome_Interaction,
                       y = abs(GI_scores))) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "yellow2") +
  geom_signif(comparisons = list(c(1,2)), # Only keep significant
              #position="dodge",
              #tip_length = 0.01,
              #step_increase = 0.05,
              #vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color ="black", linewidth = xy_plotline_width),
        axis.line.y = element_line(color ="black", linewidth = xy_plotline_width),
        axis.ticks = element_line(color = "black", linewidth = Axis_tick_width),
        axis.ticks.length = unit(Axis_tick_length, "mm"),
        axis.text.x.bottom = element_text(color = "black", size = Axis_text_size),
        axis.text.y.left = element_text(color = "black", size = Axis_text_size),
        axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
        axis.title.y.left = element_text(color = "black", size = Axis_title_size),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = xy_plotline_width),
        legend.text=element_text(size=12)) +
  scale_y_log10(limits = c(10^-4, 10^1),
                breaks = c(10^-4, 10^-2, 10^0),
                labels = c(expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0))) +
  xlab("Interaction within Large\n  Ribosomal Subunit?") +
  ylab("Positive GI Score")

Ribosome_posGI_violin_p


#### Graphs ####


KinasePlots <-
  plot_grid(
    kinase_stoich_map_p,
    kinase_IntStoich_violin_p,
    kinase_negGI_violin_p,
    kinase_posGI_violin_p,
    labels = c("h", "i", "j", "k"),
    label_size = 8,
    label_fontface = "bold",
    ncol = 4
  )

KinasePlots

MetabolismPlots <-
  plot_grid(
    enzyme_stoich_map_p,
    enzyme_IntStoich_violin_p,
    enzyme_negGI_violin_p,
    enzyme_posGI_violin_p,
    labels = c("e", "f", "g", "h"),
    label_size = 8,
    label_fontface = "bold",
    ncol = 4
  )

MetabolismPlots

RibosomePlots <-
  plot_grid(
    Ribosome_stoich_map_p,
    Ribosome_IntStoich_violin_p,
    Ribosome_negGI_violin_p,
    Ribosome_posGI_violin_p,
    labels = c("a", "b", "c", "d"),
    label_size = 8,
    label_fontface = "bold",
    ncol = 4
  )

RibosomePlots

# Combined
CombinedPlot <-
  plot_grid(NULL,
            RibosomePlots,
            NULL,
            MetabolismPlots,
            NULL,
            KinasePlots,
            labels = c("Large Ribosomal Subunit", "",
                       "Metabolic enzymes", "",
                       "Kinases", ""),
            label_size = 8,
            hjust = c(-0.1, -0.1, -0.1, -0.1, -0.2, -0.2),
            rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1),
            ncol = 1)

png("../results/Figures/Extended_Data_Figures/SupplementaryFigure4.png", units="mm", width=180, height=120, res=300)
CombinedPlot
dev.off()

#mm_to_inch = 0.0393701
#svg("../graphs/Supplementary_Figure_properdim.svg", width=180*mm_to_inch, height=120*mm_to_inch)
#CombinedPlot
#dev.off()
