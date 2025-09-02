# Xavier Castellanos-Girouard

# Date First Created: March 14 2024 
# Date Last Modified: September 4 2024

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

## Stoichiometry/GI dataset with Paralog Type
Stoich_GI_AllPara_DF <- 
  read.csv("../results/Paralog_Analyses/Stoich_GI_ParalogType.csv",
           row.names = 1)


#### Set Text Sizes ####

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


#### Kd Graphs ####

Kd_GI_AllPara_DF$Type <-
  factor(Kd_GI_AllPara_DF$Type,
         levels = c("NoParalog", "OneParalog", "TwoParalog", "IntParalog"))

## Check distribution of Kds in paralog types
Kd_dist_para_plot <- 
  ggplot(data = Kd_GI_AllPara_DF,
       mapping = aes(x = Type,
                     y = Kd,
                     fill = Type)) +
  geom_violin(trim = TRUE,
              color = "black",
              draw_quantiles = c(0.25, 0.5, 0.75),
              linewidth = geom_violin_linewidth) +
  geom_signif(comparisons = list(c(1,2), c(2,3), c(1,3), c(2,4), c(3,4)), # Only keep significant
              #position="dodge",
              tip_length = 0.02,
              y_position = c(-1, -1, -0.4, 0.3, -1),
              #step_increase = 0.05,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = Axis_text_size), 
    axis.text.y = element_text(colour="black", size = Axis_text_size),
    axis.ticks = element_line(color="black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.line.x.bottom=element_line(color="black", linewidth = xy_plotline_width),
    axis.line.y.left=element_line(color="black", linewidth = xy_plotline_width),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position = "none") +
  scale_fill_manual(values = c("grey50", 
                                "royalblue2", 
                                "green2",
                                "violet")) +
  ylab("Kd (M)") +
  xlab("Paralog Type") +
  scale_y_log10(limits = c(10^-12, 10^1),
                breaks = c(10^-12, 10^-10, 10^-8, 10^-6, 10^-4, 10^-2, 10^0),
                labels = c(expression(10^"-12"), expression(10^"-10"), 
                           expression(10^"-8"), expression(10^"-6"),
                           expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0))) +
  scale_x_discrete(labels = c("None", "One", "Two", "Paralogous\nHeteromer"))
Kd_dist_para_plot


#### GI distribution in paralog types ####

## Check distribution of GIs in paralog types
GI_dist_para_plot <- 
  ggplot(data = Kd_GI_AllPara_DF,
       mapping = aes(x = Type,
                     y = abs(scores),
                     fill = Type)) +
  geom_violin(trim = TRUE,
              color = "black",
              draw_quantiles = c(0.25, 0.5, 0.75),
              linewidth = geom_violin_linewidth) +
  geom_signif(comparisons = list(c(1,2), c(2,3), c(1,3), c(2,4), c(3,4), c(1,4)), # Only keep significant
              #position="dodge",
              tip_length = 0.02,
              y_position = c(0, 0, 0.25, 0.5, 0, 0.75),
              #step_increase = 0.05,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = Axis_text_size), 
    axis.text.y = element_text(colour="black", size = Axis_text_size),
    axis.ticks = element_line(color="black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.line.x.bottom=element_line(color="black", linewidth = xy_plotline_width),
    axis.line.y.left=element_line(color="black", linewidth = xy_plotline_width),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),
    legend.position = "none") +
  scale_fill_manual(values = c("grey50", 
                               "royalblue2", 
                               "green2",
                               "violet")) +
  ylab("|GI Score|") +
  xlab("Paralog Type") +
  scale_x_discrete(labels = c("None", "One", "Two", "Paralogous\nHeteromer")) +
  scale_y_log10(limits = c(10^-4, 10^1),
                breaks = c(10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1),
                labels = c(expression(10^"-4"), expression(10^"-3"), 
                           expression(10^"-2"), expression(10^"-1"),
                           expression(10^0), expression(10^1)))

#GI_dist_para_plot


#### Check how Kd separates according to sequence divergence (Two paralog) ####

# For Two paralog, both ka values are needed. Filter NA values
TwoPara_Ka_DF <- 
  Kd_GI_AllPara_DF %>%
  dplyr::filter((Type == "TwoParalog") & (!is.na(source_Ka) & !is.na(target_Ka)))

# Make a Dataframe that will contain both High and low values (important for Violin plot)
TwoParaAll_ka_DF <- TwoPara_Ka_DF
TwoParaAll_ka_DF$Ka_div <- "All\n"

TwoPara_Ka_DF$Ka_div <- NA

Ka_thresh <- 0.2

## For two paralogs
TwoPara_Ka_DF$Ka_div[((TwoPara_Ka_DF$source_Ka > Ka_thresh) & (TwoPara_Ka_DF$target_Ka > Ka_thresh))] <- "High"
TwoPara_Ka_DF$Ka_div[((TwoPara_Ka_DF$source_Ka < Ka_thresh) | (TwoPara_Ka_DF$target_Ka < Ka_thresh))] <- "Low"

TwoPara_Ka_DF <- rbind(TwoPara_Ka_DF, TwoParaAll_ka_DF)

TwoPara_SeqDiv_p <- 
  ggplot(data = TwoPara_Ka_DF,
         mapping = aes(x = Ka_div,
                       y = Kd)) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "green2",
              linewidth = geom_violin_linewidth) +
  geom_signif(comparisons = list(c(2,3)), # Only keep significant
              #position="dodge",
              tip_length = 0.01,
              step_increase = 0.05,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = Axis_text_size), 
    axis.text.y = element_text(colour="black", size = Axis_text_size),
    axis.ticks = element_line(color="black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.line.x.bottom=element_line(color="black", linewidth = xy_plotline_width),
    axis.line.y.left=element_line(color="black", linewidth = xy_plotline_width),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),) +
  #legend.position = "none") +
  scale_y_log10(limits = c(10^-12, 10^1),
                breaks = c(10^-12, 10^-10, 10^-8, 10^-6, 10^-4, 10^-2, 10^0),
                labels = c(expression(10^"-12"), expression(10^"-10"), 
                           expression(10^"-8"), expression(10^"-6"),
                           expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0))) +
  xlab("Sequence Divergence (dN)") +
  ylab("Kd (M)")


#tiff("../graphs/Paralogs_Kd_vs_dN.tiff", units = "mm", height = 120, width = 180, res = 300)
#TwoPara_SeqDiv_p
#dev.off()

#### Check how Kd separates according to sequence divergence (One paralog) ####

# For One paralog, only one ka values is needed
OnePara_Ka_DF <- 
  Kd_GI_AllPara_DF %>%
  dplyr::filter((Type == "OneParalog") & (!is.na(source_Ka) | !is.na(target_Ka)))

# Necessary to avoid bugs
OnePara_Ka_DF$ka_merged <- pmax(OnePara_Ka_DF$source_Ka, OnePara_Ka_DF$target_Ka, na.rm = TRUE)

OnePara_All_Ka_DF <- OnePara_Ka_DF
OnePara_All_Ka_DF$ka_merged <- NA
OnePara_All_Ka_DF$Ka_div <- "All\n"

Ka_thresh <- 0.2

# Assign High and low
OnePara_Ka_DF$Ka_div[(OnePara_Ka_DF$ka_merged > Ka_thresh)] <- "High"
OnePara_Ka_DF$Ka_div[(OnePara_Ka_DF$ka_merged < Ka_thresh)] <- "Low"

OnePara_Ka_DF <- rbind(OnePara_Ka_DF, OnePara_All_Ka_DF)

plot_SeqDiv_OnePara <- 
  ggplot(data = OnePara_Ka_DF,
         mapping = aes(x = Ka_div,
                       y = Kd)) +
  geom_violin(draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black",
              fill = "royalblue2",
              linewidth = geom_violin_linewidth) +
  geom_signif(comparisons = list(c(2,3)), # Only keep significant
              #position="dodge",
              tip_length = 0.01,
              step_increase = 0.05,
              vjust = 0.6,
              map_signif_level=TRUE,
              test = "wilcox.test",
              size = geom_signif_linewidth,
              textsize = geom_signif_starsize) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black", size = Axis_text_size), 
    axis.text.y = element_text(colour="black", size = Axis_text_size),
    axis.ticks = element_line(color="black", linewidth = Axis_tick_width),
    axis.ticks.length = unit(Axis_tick_length, "mm"),
    axis.line.x.bottom=element_line(color="black", linewidth = xy_plotline_width),
    axis.line.y.left=element_line(color="black", linewidth = xy_plotline_width),
    axis.title.x.bottom = element_text(color = "black", size = Axis_title_size),
    axis.title.y.left = element_text(color = "black", size = Axis_title_size),) +
  #legend.position = "none") +
  scale_y_log10(limits = c(10^-12, 10^1),
                breaks = c(10^-12, 10^-10, 10^-8, 10^-6, 10^-4, 10^-2, 10^0),
                labels = c(expression(10^"-12"), expression(10^"-10"), 
                           expression(10^"-8"), expression(10^"-6"),
                           expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0))) +
  xlab("Sequence Divergence (dN)") +
  ylab("Kd (M)")


#tiff("../graphs/Paralogs_Kd_vs_dN.tiff", units = "mm", height = 120, width = 180, res = 300)
#plot_SeqDiv_OnePara
#dev.off()



#### Make plots ####


Combined_plot <-
  plot_grid(GI_dist_para_plot,
            Kd_dist_para_plot,
            plot_SeqDiv_OnePara,
            TwoPara_SeqDiv_p,
            labels = c("a", "b", "c", "d"),
            label_size = 8,
            rel_widths = c(1, 1, 0.75, 0.75),
            ncol = 4)

png("../results/Figures/Main_Figures/Figure_3.png", units="mm", width=180, height=40, res=300)
Combined_plot
dev.off()

mm_to_inch = 0.0393701
svg("../results/Figures/Main_Figures/Figure_3.svg", width=180*mm_to_inch, height=40*mm_to_inch)
Combined_plot
dev.off()
