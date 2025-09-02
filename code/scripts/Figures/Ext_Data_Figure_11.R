# Xavier Castellanos-Girouard
# 
# Analysis of localization for paralog proteins
#
# Date First Created: August 4th 2025 
# Date Last Modified: August 5th 2025

#### Import Libraries ####

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(stringr)
library(ggstatsplot) # Make sure you have: "sudo apt-get install libgmp-dev" and "sudo apt-get install libmpfr-dev"
library(cowplot)

#### Set Text Sizes ####

Axis_text_size = 7
Axis_title_size = 8
plot_title_size = 7

Axis_tick_width = 0.5
Axis_tick_length = 1

geom_point_size = 0.5
geom_violin_linewidth = 0.2

xy_plotline_width = 0.5

geom_signif_linewidth = 0.3
geom_signif_starsize = 2.5


#### Import Data ####
Yeast_GFP_DF <- read.csv("../data/Miscellaneous/YeastGFP.tsv", sep = "\t")
Yeast_GFP_DF <- Yeast_GFP_DF %>% select(-c("X", "X.1"))

Yeast_Kd_GI_Prg_DF <- read.csv("../results/Paralog_Analyses/Kd_GI_Divergence.csv", row.names = 1)
Yeast_paralogs_DF <- read.csv("../results/Paralog_Analyses/all_paralog_kaks.csv", row.names = 1)


#### Format and merge data ####
Yeast_GFP_DF <- 
  Yeast_GFP_DF %>%
  select("yORF", "localization.summary")

Yeast_GFP_DF$localization.summary <- 
  str_split(Yeast_GFP_DF$localization.summary, pattern = ",")

Yeast_Prg_GFP_DF <-
  merge(x = Yeast_Kd_GI_Prg_DF,
        y = Yeast_GFP_DF,
        by.x = "source",
        by.y = "yORF",
        all.x = FALSE,
        all.y = FALSE)

Yeast_Prg_GFP_DF <-
  merge(x = Yeast_Prg_GFP_DF,
        y = Yeast_GFP_DF,
        by.x = "target",
        by.y = "yORF",
        all.x = FALSE,
        all.y = FALSE)

colnames(Yeast_Prg_GFP_DF)[23:24] <- c("source_localization", "target_localization")

## format paralog list 
Yeast_paralogs_DF <-
  Yeast_paralogs_DF %>%
  separate_wider_delim(cols = "paralogs",
                       delim = "_",
                       names = c("paralog1", "paralog2")) %>%
  select(c("paralog1", "paralog2"))

## Assign paralog to source and target
Yeast_Prg_GFP_DF$source_paralog <-
  sapply(Yeast_Prg_GFP_DF$source,
         FUN = function(x){
           comb1 <- Yeast_paralogs_DF$paralog2[which(Yeast_paralogs_DF$paralog1 == x)]
           comb2 <- Yeast_paralogs_DF$paralog1[which(Yeast_paralogs_DF$paralog2 == x)]
           if (length(comb1)==1){return(comb1)}
           else if (length(comb2)==1){return(comb2)}
           else {return(NA)}
         })
Yeast_Prg_GFP_DF$target_paralog <-
  sapply(Yeast_Prg_GFP_DF$target,
         FUN = function(x){
           comb1 <- Yeast_paralogs_DF$paralog2[which(Yeast_paralogs_DF$paralog1 == x)]
           comb2 <- Yeast_paralogs_DF$paralog1[which(Yeast_paralogs_DF$paralog2 == x)]
           if (length(comb1)==1){return(comb1)}
           else if (length(comb2)==1){return(comb2)}
           else {return(NA)}
         })

## Add localization information to paralogs
Yeast_Prg_GFP_DF <-
  merge(x = Yeast_Prg_GFP_DF,
        y = Yeast_GFP_DF,
        by.x = "source_paralog",
        by.y = "yORF",
        all.x = TRUE,
        all.y = FALSE)

Yeast_Prg_GFP_DF <-
  merge(x = Yeast_Prg_GFP_DF,
        y = Yeast_GFP_DF,
        by.x = "target_paralog",
        by.y = "yORF",
        all.x = TRUE,
        all.y = FALSE)

colnames(Yeast_Prg_GFP_DF)[27:28] <- c("source_paralog_localization", "target_paralog_localization")

## Determine if source and target are colocalized with their paralogs
Yeast_Prg_GFP_DF$source_colocalized <- 
  sapply(seq_along(Yeast_Prg_GFP_DF$source),
         FUN = function(i){
           colocalization <- 
             intersect(Yeast_Prg_GFP_DF$source_localization[[i]],
                       Yeast_Prg_GFP_DF$source_paralog_localization[[i]])
           if (length(colocalization)>0){TRUE}else{return(FALSE)}
         })

Yeast_Prg_GFP_DF$target_colocalized <- 
  sapply(seq_along(Yeast_Prg_GFP_DF$target),
         FUN = function(i){
           colocalization <- 
             intersect(Yeast_Prg_GFP_DF$target_localization[[i]],
                       Yeast_Prg_GFP_DF$target_paralog_localization[[i]])
           if (length(colocalization)>0){TRUE}else{return(FALSE)}
         })

Yeast_Prg_GFP_DF$Colocalized <- Yeast_Prg_GFP_DF$source_colocalized | Yeast_Prg_GFP_DF$target_colocalized

## Assign divergence group (high/low)
Yeast_Prg_GFP_DF$Divergence <- NA
Yeast_Prg_GFP_DF$Divergence[(Yeast_Prg_GFP_DF$source_Ka > 0.2) & (Yeast_Prg_GFP_DF$target_Ka > 0.2)] = "High"
Yeast_Prg_GFP_DF$Divergence[(Yeast_Prg_GFP_DF$source_Ka > 0.2) & is.na(Yeast_Prg_GFP_DF$target_Ka)] = "High"
Yeast_Prg_GFP_DF$Divergence[is.na(Yeast_Prg_GFP_DF$source_Ka) & (Yeast_Prg_GFP_DF$target_Ka > 0.2)] = "High"
Yeast_Prg_GFP_DF$Divergence[(Yeast_Prg_GFP_DF$source_Ka < 0.2)] = "Low"
Yeast_Prg_GFP_DF$Divergence[(Yeast_Prg_GFP_DF$target_Ka < 0.2)] = "Low"

Yeast_Prg_GFP_DF <- Yeast_Prg_GFP_DF[!is.na(Yeast_Prg_GFP_DF$Divergence),]

Yeast_Prg_GFP_DF <- 
  Yeast_Prg_GFP_DF %>%
  mutate(
    Divergence = factor(Divergence, levels = c("High", "Low")),
    Colocalized = factor(Colocalized, levels = c("TRUE", "FALSE")))

#### Graphs ####

comparisons <- list(c("FALSE", "TRUE"))
paralogLoc_p <- 
  ggplot() +
  geom_violin(data = Yeast_Prg_GFP_DF,
              mapping = aes(x = Colocalized,
                            y = Kd,
                            fill = Divergence),
              draw_quantiles = c(0.25, 0.50, 0.75),
              trim = TRUE,
              color = "black") +
  facet_wrap(~ Divergence) +
  geom_signif(
    data = Yeast_Prg_GFP_DF %>% filter(Divergence == "High"),
    aes(x = Colocalized, y = Kd),
    comparisons = list(c("FALSE", "TRUE")),
    map_signif_level = TRUE,
    tip_length = 0.01,
    textsize = 3,
    inherit.aes = FALSE
  ) +
  geom_signif(
    data = Yeast_Prg_GFP_DF %>% filter(Divergence == "Low"),
    aes(x = Colocalized, y = Kd),
    comparisons = list(c("FALSE", "TRUE")),
    map_signif_level = TRUE,
    tip_length = 0.01,
    textsize = 3,
    inherit.aes = FALSE
  ) +
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
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.text = element_text(colour = "black", 
                               size = Axis_text_size),
    legend.title = element_text(colour = "black",
                                size = Axis_title_size),
    #legend.text.align = 0,
    #legend.key = element_rect(fill = "white"),
    #legend.position = "none",
    #legend.key = element_rect(fill = "black", 
    #                          linewidth = 0),
    #legend.key.size = unit(3, "mm"),
    ) +
  #legend.position = "none") +
  scale_y_log10(limits = c(10^-12, 10^1),
                breaks = c(10^-12, 10^-10, 10^-8, 10^-6, 10^-4, 10^-2, 10^0),
                labels = c(expression(10^"-12"), expression(10^"-10"), 
                           expression(10^"-8"), expression(10^"-6"),
                           expression(10^"-4"), expression(10^"-2"), 
                           expression(10^0))) +
  xlab("Are Paralogs Colocalized?")


#### Fishers exact test ####

##                High Div.  Low Div.
## Colocal          379       514
## Not Colocal      209       128

## Masiac plot
LocTesting <- data.frame(
  "High Divergence" = c(Yeast_Prg_GFP_DF %>% filter((Divergence == "High") & (Colocalized == "TRUE")) %>% nrow(), 
                        Yeast_Prg_GFP_DF %>% filter((Divergence == "High") & (Colocalized == "FALSE")) %>% nrow()),
  "Low Divergence" = c(Yeast_Prg_GFP_DF %>% filter((Divergence == "Low") & (Colocalized == "TRUE")) %>% nrow(), 
                       Yeast_Prg_GFP_DF %>% filter((Divergence == "Low") & (Colocalized == "FALSE")) %>% nrow()),
  row.names = c("Colocalized", "Not Colocalized"),
  stringsAsFactors = FALSE
)
colnames(LocTesting) <- c("High Divergence", "Low Divergence")

mosaicplot(LocTesting,
           main = "Mosaic plot",
           color = TRUE
)

fisher.test(LocTesting, alternative = "less")

## Mosaic plot with values
x <- c()
for (row in rownames(LocTesting)) {
  for (col in colnames(LocTesting)) {
    x <- rbind(x, matrix(rep(c(row, col), LocTesting[row, col]), ncol = 2, byrow = TRUE))
  }
}
df <- as.data.frame(x)
colnames(df) <- c("Colocalization", "Sequence Divergence")
df

test <- fisher.test(table(df))

# combine plot and statistical test with ggbarstats
p <- ggbarstats(
  df, Colocalization, `Sequence Divergence`,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    formatC(test$p.value, format = "e", digits = 2)
  )
)  +
  theme(
    axis.text.x = element_text(color = "black", size = Axis_text_size),  # X-axis: High/Low Divergence
    axis.text.y = element_text(color = "black", size = Axis_text_size),  # Y-axis tick labels
    axis.title.x = element_text(color = "black", size = Axis_title_size, face = "plain"),
    axis.title.y = element_text(color = "black", size = Axis_title_size),
    legend.title = element_text(size = Axis_title_size, face = "plain"),
    legend.text = element_text(size = Axis_title_size),
    plot.subtitle = element_text(size = Axis_title_size)
  )

p$layers[[2]]$aes_params$size <- 3  # percentages
p$layers[[3]]$aes_params$size <- 3  # counts (n=)
p

# Interpretation: The odds of being colocalized are lower in High Divergence than in Low Divergence

### combine plots 
Combined_plot <-
  plot_grid(p,
            paralogLoc_p,
            labels = c("a", "b"),
            label_size = 8,
            ncol = 2)

png("../results/Figures/Extended_Data_Figures/SupplementaryFigure11.png", units="mm", width=180, height=90, res=300)
Combined_plot
dev.off()

#mm_to_inch = 0.0393701
#svg("../graphs/ParalogLoc_properdim.svg", width=180*mm_to_inch, height=90*mm_to_inch)
#Combined_plot
#dev.off()
