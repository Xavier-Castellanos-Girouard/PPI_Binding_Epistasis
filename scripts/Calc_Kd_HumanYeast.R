# Xavier Castellanos-Girouard

# Date First Created: Sometime Summer 2022
# Date Last Modified: October 15 2024


#### Import Libraries ####

library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)

#### Import data ####

dir <- "/home/xavier/Desktop/Cell_interactome_stoichiometry/Yeast_Epistasis_Stoichiometries/Third_Run/"

# Import human stoichiometries
Human_Stoich_DF <- read.csv(paste0(dir, "R_Stoichiometries_vs_GI/results/Human_Interactome_Stoich.csv"), row.names = 1)

# Import abundance data for Hein et al. 
Hein_Proteome <- readxl::read_xlsx(paste0(dir, "R_Calc_Kd_Benchmark/data/Hein_TableS3.xlsx"), sheet = 2)

# Import abundance data for Opencell
# From https://opencell.czbiohub.org/download
OpenCell_Proteome <- read.csv(paste0(dir, "R_Calc_Kd_Benchmark/data/opencell-protein-abundance.csv"))

# Import human litterature kd values
Human_Litt_Kd_DF <- read.csv(paste0(dir, "R_Calc_Kd_Benchmark/results/Human_Kd_literature.csv"), row.names = 1)

# Computational dG
file_list = list.files(path = paste0(dir, "R_Calc_Kd_Benchmark/results/PRODIGY_outputs/"), pattern="\\.out$")
tsv_files = lapply(file_list, 
                   function(x){
                     new_dir <- paste0(dir, "R_Calc_Kd_Benchmark/results/PRODIGY_outputs/", x)
                     return(read.csv(new_dir, sep ="\t"))
                   })

# Import yeast stoichiometries
Yeast_Stoich_DF <- read.csv(paste0(dir, "R_Calc_Stoichiometries/results/Yeast_Interactome_Stoich.csv"), row.names = 1)

# Import yeast litterature kd values
Yeast_Litt_Kd_DF <- read.csv(paste0(dir, "R_Calc_Kd_Benchmark/results/Yeast_Kd_literature.csv"), row.names = 1)

# Import yeast PCA data
Yeast_PCA_Kd_DF <- read.csv(paste0(dir, "R_Calc_Kd_Benchmark/results/Estimated_Kd2_PCA_intensitiesNorm_LinFit_MeanCopyNumberUsed.csv"), row.names = 1)
colnames(Yeast_PCA_Kd_DF)[14] <- "PCA_Kd"

#### Calculating Kd for Hein data ####

# Filter stoichiometry data for Hein
Hein_Stoich_DF <- 
  Human_Stoich_DF %>%
  dplyr::filter(Dataset == "Hein")

# Separate genes names when there are several in a same row
Hein_Proteome <-
  Hein_Proteome %>%
  dplyr::select(c(3,17,18)) %>%
  separate_longer_delim(cols = "Gene.name",
                        delim = ";")

# Sum the copy number different proteins forms for each gene
Hein_Proteome <-
  Hein_Proteome %>%
  dplyr::group_by(Gene.name) %>%
  dplyr::summarise(Copy.number = sum(Copy.number)) %>%
  as.data.frame()

## Merge Abundance Values to hein PPI dataframe

# Merge bait abundances
Hein_Stoich_DF <-
  merge(x = Hein_Stoich_DF,
        y = Hein_Proteome[,c(1,2)],
        by.x = "Bait",
        by.y = "Gene.name",
        all.x = FALSE,
        all.y = FALSE)

colnames(Hein_Stoich_DF)[6] <- "Bait_Abundance"

# merge prey abundances
Hein_Stoich_DF <-
  merge(x = Hein_Stoich_DF,
        y = Hein_Proteome[,c(1,2)],
        by.x = "Prey",
        by.y = "Gene.name",
        all.x = FALSE,
        all.y = FALSE)

colnames(Hein_Stoich_DF)[7] <- "Prey_Abundance"


## Volume of HeLa Cell (doi: 10.1002/nbm.1173) 2.6 × 10^3 μm^3
HeLa_vol <- (2.6*10^3)*10^-15 # Convert micrometer cubed to liters
Avogadro <- 6.023*10^23

# Estimate Complex Copy Numbers
Hein_Stoich_DF$Complex_CopyNumber <-
  Hein_Stoich_DF$Interaction_Stoichiometry*Hein_Stoich_DF$Bait_Abundance

# Estimate Free Abundances (in Copy Numbers)
Hein_Stoich_DF$Bait_FreeAbundance <- Hein_Stoich_DF$Bait_Abundance - Hein_Stoich_DF$Complex_CopyNumber
Hein_Stoich_DF$Prey_FreeAbundance <- Hein_Stoich_DF$Prey_Abundance - Hein_Stoich_DF$Complex_CopyNumber

## Convert Copy Numbers to concentrations
Hein_Stoich_DF$Complex_Conc <- (Hein_Stoich_DF$Complex_CopyNumber/Avogadro)/HeLa_vol
Hein_Stoich_DF$Bait_FreeAbundance_Conc <- (Hein_Stoich_DF$Bait_FreeAbundance/Avogadro)/HeLa_vol
Hein_Stoich_DF$Prey_FreeAbundance_Conc <- (Hein_Stoich_DF$Prey_FreeAbundance/Avogadro)/HeLa_vol

# Calculate estimated Kd (Unit is mol/L)
Hein_Stoich_DF$Kd <- Hein_Stoich_DF$Bait_FreeAbundance_Conc*Hein_Stoich_DF$Prey_FreeAbundance_Conc/Hein_Stoich_DF$Complex_Conc

# Remove negative Kd's, or positive Kd's that are a result of multiplying two
#   negative free abundance values
Hein_Stoich_DF <- 
  Hein_Stoich_DF %>%
  dplyr::filter(Kd > 0) %>%
  dplyr::filter((Prey_FreeAbundance > 0) & (Bait_FreeAbundance > 0))


#### Calculate Kd for OpenCell Interactions ####


# Filter stoichiometry data for Opencell
Opencell_Stoich_DF <- 
  Human_Stoich_DF %>%
  dplyr::filter(Dataset == "OpenCell")

## Merge Abundance Values to dataframe
Opencell_Stoich_DF <-
  merge(x = Opencell_Stoich_DF,
        y = OpenCell_Proteome[,c(1,8)],
        by.x = "Bait",
        by.y = "gene_name",
        all.x = FALSE,
        all.y = FALSE)

colnames(Opencell_Stoich_DF)[6] <- "Bait_Abundance"

Opencell_Stoich_DF <-
  merge(x = Opencell_Stoich_DF,
        y = OpenCell_Proteome[,c(1,8)],
        by.x = "Prey",
        by.y = "gene_name",
        all.x = FALSE,
        all.y = FALSE)

colnames(Opencell_Stoich_DF)[7] <- "Prey_Abundance"


### Estimate Kd

# Estimate number of complexes (copy number)
Opencell_Stoich_DF$Complex_CopyNumber <-
  Opencell_Stoich_DF$Interaction_Stoichiometry*Opencell_Stoich_DF$Bait_Abundance


## Estimate Free Abundances of binders (Copy Numbers)
Opencell_Stoich_DF$Bait_FreeAbundance <- Opencell_Stoich_DF$Bait_Abundance-Opencell_Stoich_DF$Complex_CopyNumber
Opencell_Stoich_DF$Prey_FreeAbundance <- Opencell_Stoich_DF$Prey_Abundance-Opencell_Stoich_DF$Complex_CopyNumber


## Calculate Concentrations

# Diameter of a HEK293 cell: 13 micrometers, Bionumbers ID: 108893
# Estimated volume (perfect sphere): 1.15*10^3 um^3

HEK_vol <- (1.15*10^3)*10^-15 # Convert micrometer cubed to liters
Avogadro <- 6.023*10^23

Opencell_Stoich_DF$Complex_Conc <- (Opencell_Stoich_DF$Complex_CopyNumber/Avogadro)/HEK_vol
Opencell_Stoich_DF$Bait_FreeAbundance_Conc <- (Opencell_Stoich_DF$Bait_FreeAbundance/Avogadro)/HEK_vol
Opencell_Stoich_DF$Prey_FreeAbundance_Conc <- (Opencell_Stoich_DF$Prey_FreeAbundance/Avogadro)/HEK_vol

# Estimate Kd
Opencell_Stoich_DF$Kd <- 
  Opencell_Stoich_DF$Bait_FreeAbundance_Conc*Opencell_Stoich_DF$Prey_FreeAbundance_Conc/Opencell_Stoich_DF$Complex_Conc


# Remove negative Kd's or Kd's resulting from negative concentrations
Opencell_Stoich_DF <- 
  Opencell_Stoich_DF %>%
  dplyr::filter(Kd > 0) %>%
  dplyr::filter((Prey_FreeAbundance > 0) & (Bait_FreeAbundance > 0))


# Remove NA Kd values
Opencell_Stoich_DF <- 
  Opencell_Stoich_DF %>%
  dplyr::filter(!is.na(Kd))

#### Verify Estimated Kd values with that of literature (Human) ####

## Combine Hein and Opencell Kds
Human_Kd_DF <- rbind(Hein_Stoich_DF, Opencell_Stoich_DF)

## Merge literature and estimated Kd

# Merge by interactions. First possible combination
Human_Est_vs_Lit_DF1 <-
  merge(x = Human_Kd_DF,
        y = Human_Litt_Kd_DF,
        by.x = c("Bait", "Prey"),
        by.y = c("protein1_GeneName", "protein2_GeneName"),
        all.x = FALSE,
        all.y = FALSE)

# Merge by interactions. Second possible combination
Human_Est_vs_Lit_DF2 <-
  merge(x = Human_Kd_DF,
        y = Human_Litt_Kd_DF,
        by.x = c("Bait", "Prey"),
        by.y = c("protein2_GeneName", "protein1_GeneName"),
        all.x = FALSE,
        all.y = FALSE)

# Combine interaction merges
Human_Est_vs_Lit_DF <- rbind(Human_Est_vs_Lit_DF1, Human_Est_vs_Lit_DF2)


## Remove duplicated structures with identical estimated and litterature Kd's
# Note: This makes the correlation a little worse but is more rigorous
Human_Est_vs_Lit_DF <- Human_Est_vs_Lit_DF[!duplicated(Human_Est_vs_Lit_DF[,c(1:4,14,16)]),]


#plot(x = log10(Human_Est_vs_Lit_DF$Value),
#     y = log10(Human_Est_vs_Lit_DF$Kd),
#     xlab = "Literature Kd (M)",
#     ylab = "Estimated Kd")

cor.test(x = log10(Human_Est_vs_Lit_DF$Value),
         y = log10(Human_Est_vs_Lit_DF$Kd),
         method = "pearson")

plot(x = log10(Human_Est_vs_Lit_DF$Value),
     y = log10(Human_Est_vs_Lit_DF$Interaction_Stoichiometry),
     xlab = "Literature Kd (M)",
     ylab = "Interaction Stoichiometry")

cor.test(x = log10(Human_Est_vs_Lit_DF$Value),
         y = log10(Human_Est_vs_Lit_DF$Interaction_Stoichiometry),
         method = "pearson")


Estimated_Kd_vs_Literature_Kd_p <- 
  ggplot(data = Human_Est_vs_Lit_DF,
         mapping = aes(x = Kd,
                       y = Value)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "red") +
  ggpubr::stat_cor(method="pearson",
                   show.legend = FALSE,
                   digits = 2) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black"), 
    axis.text.y = element_text(colour="black"),
    axis.ticks = element_line(color="black"),
    axis.line.x.bottom=element_line(color="black"),
    axis.line.y.left=element_line(color="black"),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_x_log10(breaks = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2),
                labels = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2)) +
  scale_y_log10(breaks = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3),
                labels = c(10^-9, 10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3)) +
  xlab("Estimated Kd (M)") +
  ylab("Literature Kd (M)")

#png("../graphs/Estimated_vs_Lit.png", units="mm", width=120, height=80, res=300)
#Estimated_Kd_vs_Literature_Kd_p
#dev.off()


## Export estimated Human Kd values and Est. vs Lit kd table
#write.csv(Human_Est_vs_Lit_DF, paste0(dir, "R_Calc_Kd_Benchmark/results/Human_Kd_Est_vs_Lit.csv"))

## Export estimated human Kd values
#write.csv(Human_Kd_DF, paste0(dir, "R_Calc_Kd_Benchmark/results/Human_Estimated_Kd.csv"))


#### Benchmarking with PRODIGY (computed dG) ####

PDB_ID_list <- 
  gsub(x = file_list,
       pattern = ".out",
       fixed = TRUE,
       replacement = "")

PDB_ID_list <- tolower(PDB_ID_list)

Kd_value_list <- c()
for (table in tsv_files){
  Kd_row <- nrow(table)
  Kd_info <- table[Kd_row,1]
  Kd_value <- 
    gsub(x = Kd_info,
         pattern = "[++] Predicted dissociation constant (M) at 25.0˚C:  ", 
         replacement = "",
         fixed = TRUE)
  Kd_value <- as.numeric(Kd_value)
  Kd_value_list <- c(Kd_value_list, Kd_value)
}

PRODIGY_Kd_DF <- 
  data.frame(PDB_ID = PDB_ID_list,
             Computed_Kd = Kd_value_list)

Human_Est_vs_Lit_vs_Comp_DF <-
  merge(x = Human_Est_vs_Lit_DF,
        y = PRODIGY_Kd_DF,
        by.x = "PDB_ID",
        by.y = "PDB_ID",
        all.x = TRUE,
        all.y = FALSE)

cor.test(-log10(Human_Est_vs_Lit_vs_Comp_DF$Computed_Kd),
         -log10(Human_Est_vs_Lit_vs_Comp_DF$Value))

plot(x = -log10(Human_Est_vs_Lit_vs_Comp_DF$Computed_Kd),
     y = -log10(Human_Est_vs_Lit_vs_Comp_DF$Value))

#write.csv(Human_Est_vs_Lit_vs_Comp_DF, paste0(dir, "R_Calc_Kd_Benchmark/results/Human_Kd_Est_vs_Lit_vs_Comp.csv"))

#### Calculating Kd for Yeast data ####

### Estimate number of free proteins
Yeast_Stoich_DF$Complex_CopyNumber <- Yeast_Stoich_DF$Interaction_Stoichiometry*Yeast_Stoich_DF$source_copy_number

Yeast_Stoich_DF$Source_FreeAbundance <- Yeast_Stoich_DF$source_copy_number-Yeast_Stoich_DF$Complex_CopyNumber
Yeast_Stoich_DF$Target_FreeAbundance <- Yeast_Stoich_DF$target_copy_number-Yeast_Stoich_DF$Complex_CopyNumber

### Calculate concentrations and Kd

# Haploid yeast cell volume (DOI: 10.1126/science.1070850)
volume = 4.2*10^-14 # Litres
Avogadro = 6.022*10^23 # Molecules per mole

# Calculate concentration in mol/L (M)
Yeast_Stoich_DF$Complex_Conc <- (Yeast_Stoich_DF$Complex_CopyNumber/(6.022*10^23))/8.2e-14 
Yeast_Stoich_DF$Source_FreeConc <- (Yeast_Stoich_DF$Source_FreeAbundance/(6.022*10^23))/8.2e-14
Yeast_Stoich_DF$Target_FreeConc <- (Yeast_Stoich_DF$Target_FreeAbundance/(6.022*10^23))/8.2e-14

# Calculate dissociation constant
Yeast_Stoich_DF$Kd <- (Yeast_Stoich_DF$Source_FreeConc*Yeast_Stoich_DF$Target_FreeConc)/Yeast_Stoich_DF$Complex_Conc

# Remove any negative concentrations and Kd's
Yeast_Stoich_DF <-
  Yeast_Stoich_DF %>%
  dplyr::filter((Source_FreeConc > 0) & (Target_FreeConc > 0)) %>%
  dplyr::filter(Kd > 0)

#### Verify Estimated Kd values with that of literature (Yeast) ####

## Merge literature and estimated Kd

# Merge by interactions. First possible combination
Yeast_Est_vs_Lit_DF1 <-
  merge(x = Yeast_Stoich_DF,
        y = Yeast_Litt_Kd_DF,
        by.x = c("source", "target"),
        by.y = c("protein1_ORF", "protein2_ORF"),
        all.x = FALSE,
        all.y = FALSE)

# Merge by interactions. Second possible combination
Yeast_Est_vs_Lit_DF2 <-
  merge(x = Yeast_Stoich_DF,
        y = Yeast_Litt_Kd_DF,
        by.x = c("source", "target"),
        by.y = c("protein2_ORF", "protein1_ORF"),
        all.x = FALSE,
        all.y = FALSE)

# Combine interaction merges
Yeast_Est_vs_Lit_DF <- rbind(Yeast_Est_vs_Lit_DF1, Yeast_Est_vs_Lit_DF2)


## Note: Only four values


#### Corroborate Kd values with that of DHFR PCA (Yeast) ####


Combined_DF1 <- 
  merge(x = Yeast_Stoich_DF[, colnames(Yeast_Stoich_DF) %in% c("source", "target", "Layer", "Kd")],
        y = Yeast_PCA_Kd_DF[, colnames(Yeast_PCA_Kd_DF) %in% c("MATa.orf_name", "MATalpha.orf_name", "PCA_Kd")],
        by.x = c("source", "target"),
        by.y = c("MATa.orf_name", "MATalpha.orf_name"),
        all.x = FALSE,
        all.y = FALSE)

Combined_DF2 <- 
  merge(x = Yeast_Stoich_DF[, colnames(Yeast_Stoich_DF) %in% c("source", "target", "Layer", "Kd")],
        y = Yeast_PCA_Kd_DF[, colnames(Yeast_PCA_Kd_DF) %in% c("MATa.orf_name", "MATalpha.orf_name", "PCA_Kd")],
        by.x = c("target", "source"),
        by.y = c("MATa.orf_name", "MATalpha.orf_name"),
        all.x = FALSE,
        all.y = FALSE)

Combined_DF <- rbind(Combined_DF1, Combined_DF2)


cor.test(x = log10(Combined_DF$Kd),
         y = log10(Combined_DF$PCA_Kd))

ggplot(data = Combined_DF,
       mapping = aes(x = Kd,
                     y = PCA_Kd)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "red") +
  ggpubr::stat_cor(method="pearson",
                   show.legend = FALSE,
                   digits = 2) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour="black"), 
    axis.text.y = element_text(colour="black"),
    axis.ticks = element_line(color="black"),
    axis.line.x.bottom=element_line(color="black"),
    axis.line.y.left=element_line(color="black"),
    legend.position="none",
    plot.title = element_text(hjust = 0.5)) +
  scale_x_log10(breaks = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2),
                labels = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4, 10^-3, 10^-2)) +
  scale_y_log10(breaks = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4),
                labels = c(10^-8, 10^-7, 10^-6, 10^-5, 10^-4)) +
  xlab("Mass Spec. Estimated Kd (M)") +
  ylab("PCA Estimated Kd (M)")
  

## Export Yeast
#write.csv(Yeast_Stoich_DF, paste0(dir, "R_Calc_Kd_Benchmark/results/Yeast_Estimated_Kd.csv"))
#write.csv(Combined_DF, paste0(dir, "R_Calc_Kd_Benchmark/results/Yeast_Kd_Est_vs_PCA_Kd.csv"))
