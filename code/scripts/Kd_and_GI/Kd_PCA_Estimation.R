## Xavier Castellanos-Girouard
## Date First Created: November 6 2023
## Date Last Modified: April 27 2024

# Using Baryshnikova, Tarassov and Levy data to estimate Kds


#### Import Libraries ####

library(tidyr)
library(dplyr)
library(readxl)
library(ggplot2)
#library(stringr)
#library(modelr)

#### Import data ####

# Protein abundance data Ho et al. 2018
Protein_Abundances_comp_DF <- read_xlsx("../data/ProteinAbundances/Abundances_copies_per_cell.xlsx")
Protein_Abundances_comp_colnames <- Protein_Abundances_comp_DF[2,]

# Protein abundance data from Levy et al 2014
Protein_Abundance_PCA_DF <- read_xls("../data/PCA_data/Abundace_Levy_2014.xls")
Protein_Abundance_PCA_colnames <- Protein_Abundance_PCA_DF[19,]
Protein_Abundance_PCA_DF <- Protein_Abundance_PCA_DF[20:nrow(Protein_Abundance_PCA_DF),]
Protein_Abundance_PCA_DF <- Protein_Abundance_PCA_DF[,2:ncol(Protein_Abundance_PCA_DF)]

## PPI intensity data from Tarassov et al. 2008 (including contaminants)
Protein_Interaction_PCA <- read_xls("../data/PCA_data/Tarassov_2008_PPI.xls")

# List of PPIs from Tarassov et al. 2008 (excluding contaminants)
Tarassov_only_Int_DF <- read.csv("../data/PCA_data/network_1.tsv", 
                                 sep = "\t",
                                 header = FALSE)

#### Format Data ####

### Format Compendium of Protein abundances
# Rename columns
colnames(Protein_Abundances_comp_DF) <- Protein_Abundances_comp_colnames

# Select appropriate rows
Protein_Abundances_comp_DF <- Protein_Abundances_comp_DF[3:nrow(Protein_Abundances_comp_DF),]

## Only keep useful columns
Protein_Abundances_comp_DF <-
  Protein_Abundances_comp_DF %>%
  select(c("Systematic Name", 
           #"Mean molecules per cell",
           "Median molecules per cell"))

## Remove NAs
Protein_Abundances_comp_DF <-
  Protein_Abundances_comp_DF[!is.na(Protein_Abundances_comp_DF$`Systematic Name`),]

Protein_Abundances_comp_DF <-
  Protein_Abundances_comp_DF[!is.na(Protein_Abundances_comp_DF$`Median molecules per cell`),]

## Set values to numeric (they are character type when imported)
Protein_Abundances_comp_DF$`Median molecules per cell` <- 
  as.numeric(Protein_Abundances_comp_DF$`Median molecules per cell`)


### Format Data from abundance PCA

# Unlist column names
Protein_Abundance_PCA_colnames <- unname(unlist(Protein_Abundance_PCA_colnames))

# Remove NA value
Protein_Abundance_PCA_colnames <- Protein_Abundance_PCA_colnames[!is.na(Protein_Abundance_PCA_colnames)]

# Add colnames to dataframe
colnames(Protein_Abundance_PCA_DF) <- Protein_Abundance_PCA_colnames

# Keep relevant columns
Protein_Abundance_PCA_DF <-
  Protein_Abundance_PCA_DF %>%
  dplyr::select(ORF,
                GENE.NAME,
                GROWTH.CYTO,
                AB.CYTO.AGENT)

# Convert to numeric
Protein_Abundance_PCA_DF$GROWTH.CYTO <- as.numeric(Protein_Abundance_PCA_DF$GROWTH.CYTO)

Protein_Abundance_PCA_DF$AB.CYTO.AGENT <- as.numeric(Protein_Abundance_PCA_DF$AB.CYTO.AGENT)

## Remove empty and NA ORFs

# NA containing ORFs
NA_ORF_bool <- grepl(pattern = "NA", Protein_Abundance_PCA_DF$ORF)
Protein_Abundance_PCA_DF <- Protein_Abundance_PCA_DF[!(NA_ORF_bool),]

# Empty ORFs
Empty_ORF_bool <- grepl(pattern = "Empty", Protein_Abundance_PCA_DF$ORF)
Protein_Abundance_PCA_DF <- Protein_Abundance_PCA_DF[!(Empty_ORF_bool),]

Protein_Abundance_PCA_DF <- Protein_Abundance_PCA_DF %>% drop_na()

#### Plot Levy and Baryshnikova ####

## Merge abundances measured by PCA and 
Protein_Abundances_merged_DF <- 
  merge(x = Protein_Abundances_comp_DF,
        y = Protein_Abundance_PCA_DF,
        by.x = "Systematic Name",
        by.y = "ORF",
        all.x = TRUE,
        all.y = FALSE)


# Add a log10 copy number column to the dataframe for the model
Protein_Abundances_merged_DF$log10_Median_molecules_per_cell <- log10(Protein_Abundances_merged_DF$`Median molecules per cell`)

# add log10 growth (intensity) column
Protein_Abundances_merged_DF$log10_GROWTH.CYTO <- log10(Protein_Abundances_merged_DF$GROWTH.CYTO)

# Create regression linear model (about 0.61 r-squared)
lin_model <- lm(log10_Median_molecules_per_cell~GROWTH.CYTO, data=Protein_Abundances_merged_DF)

#summary(lin_model)

#### Format Tarassov data ####

## Combine Intensities with Interaction data 
Protein_Interaction_PCA_1 <-
  merge(x = Protein_Interaction_PCA,
        y = Tarassov_only_Int_DF, 
        by.x = c("MATa gene_name",
                 "MATalpha_gene_name"),
        by.y = c("V1", 
                 "V2"),
        all.x = FALSE,
        all.y = FALSE)

Protein_Interaction_PCA_2 <-
  merge(x = Protein_Interaction_PCA,
        y = Tarassov_only_Int_DF, 
        by.x = c("MATa gene_name",
                 "MATalpha_gene_name"),
        by.y = c("V2", 
                 "V1"),
        all.x = FALSE,
        all.y = FALSE)

# Combine merged dataframes
Protein_Interaction_PCA_filtered <- rbind(Protein_Interaction_PCA_1, Protein_Interaction_PCA_2)

# Run intensity and zscore filter
Protein_Interaction_PCA_filtered <-
  Protein_Interaction_PCA_filtered %>%
  dplyr::filter(intensity > 23000) %>%
  dplyr::filter(`z-score` > 2.4)

# Remove identical entries (i.e. homodimers which are included in both merges)
Protein_Interaction_PCA_filtered <-
  Protein_Interaction_PCA_filtered %>%
  distinct()

## Select useful columns
Protein_Interaction_PCA_filtered <-
  Protein_Interaction_PCA_filtered %>%
  dplyr::select(`MATa orf_name`,
                `MATalpha orf_name`,
                intensity)

### Find instances where interaction A:B and B:A are present, replace with mean of both

## Iteration through dataframe, keep track of indices of duplicates and their values
row.names(Protein_Interaction_PCA_filtered) <- NULL
duplicates_DF <- Protein_Interaction_PCA_filtered[-c(1:nrow(Protein_Interaction_PCA_filtered)),]
index_list <- c()
for (i in 1:(nrow(Protein_Interaction_PCA_filtered)-1)){
  interactor1 <- Protein_Interaction_PCA_filtered$`MATa orf_name`[i]
  interactor2 <- Protein_Interaction_PCA_filtered$`MATalpha orf_name`[i]
  
  #print(paste(interactor1, interactor2))
  #print(i)
  
  if (interactor1 == interactor2){
    next
  }
  
  j = i + 1
  #print(j)
  for (k in 1:j){
    interactor3 <- Protein_Interaction_PCA_filtered$`MATalpha orf_name`[k]
    interactor4 <- Protein_Interaction_PCA_filtered$`MATa orf_name`[k]
    intensity <- Protein_Interaction_PCA_filtered$intensity[k]
    
    isSame_interactor_1_3 <- interactor1 == interactor3
    isSame_interactor_2_4 <- interactor2 == interactor4
    
    if ((isSame_interactor_1_3 & isSame_interactor_2_4)){
      duplicates_DF <- rbind(duplicates_DF, 
                             Protein_Interaction_PCA_filtered[i,],
                             c(interactor3, interactor4, intensity))
      index_list <- c(index_list, i, k)
    }
    }
  }

duplicates_DF$intensity <- as.numeric(duplicates_DF$intensity)

# Take mean intensity of interactions with duplicate measurements
duplicates_mean_DF <-
  duplicates_DF %>%
  dplyr::group_by(`MATa orf_name`, `MATalpha orf_name`) %>%
  dplyr::summarize(intensity = mean(intensity))

# Temporarily remove duplicate measurements from dataset
Protein_Interaction_PCA_filtered <- 
  Protein_Interaction_PCA_filtered[-c(index_list), ]

# Replace duplicates with a single interaction with mean of intensity
Protein_Interaction_PCA_filtered_final <- 
  rbind(Protein_Interaction_PCA_filtered,
        duplicates_mean_DF)

#hist(log10(Protein_Interaction_PCA_filtered_final$intensity), nclass = 30)

#### Apply linear regression model to Tarassov data (ELevy Conc) ####


## Adding abundance value for individual proteins
Protein_Interaction_PCA_filtered_final <-
  merge(x = Protein_Interaction_PCA_filtered_final,
        y = Protein_Abundance_PCA_DF[, c("ORF", "AB.CYTO.AGENT")],
        by.x = "MATa orf_name",
        by.y = "ORF",
        all.x = TRUE,
        all.y = FALSE)

colnames(Protein_Interaction_PCA_filtered_final)[4] <- "MATa_Copy_Number"

Protein_Interaction_PCA_filtered_final <-
  merge(x = Protein_Interaction_PCA_filtered_final,
        y = Protein_Abundance_PCA_DF[, c("ORF", "AB.CYTO.AGENT")],
        by.x = "MATalpha orf_name",
        by.y = "ORF",
        all.x = TRUE,
        all.y = FALSE)

colnames(Protein_Interaction_PCA_filtered_final)[5] <- "MATalpha_Copy_Number"


## Naive rescale using Emmanuel Levy Data
max_range = max(Protein_Abundance_PCA_DF$GROWTH.CYTO, na.rm = TRUE)
min_range = min(Protein_Abundance_PCA_DF$GROWTH.CYTO, na.rm = TRUE)
new_range = max_range-min_range
old_range = (max(Protein_Interaction_PCA_filtered_final$intensity) - min(Protein_Interaction_PCA_filtered_final$intensity))

Protein_Interaction_PCA_filtered_final$normalized_intensity <-
  new_range*((Protein_Interaction_PCA_filtered_final$intensity - min(Protein_Interaction_PCA_filtered_final$intensity))/old_range)# + min_range

Protein_Interaction_PCA_filtered_final_kd <-
  Protein_Interaction_PCA_filtered_final %>%
  drop_na()

## Apply model

## Exponential
#Protein_Interaction_PCA_filtered_final_kd$log10_Copy_number = 7.223e-01*(2)^(3.589e-05*(Protein_Interaction_PCA_filtered_final_kd$normalized_intensity))-(-2.477)

## Linear
Regression_intercept <- lin_model$coefficients[1]
Regression_Slope <- lin_model$coefficients[2]
Protein_Interaction_PCA_filtered_final_kd$log10_Copy_number = (Regression_Slope*Protein_Interaction_PCA_filtered_final_kd$normalized_intensity)+Regression_intercept


## Convert copy numbers to linear space values
Protein_Interaction_PCA_filtered_final_kd$Copy_number = 10^Protein_Interaction_PCA_filtered_final_kd$log10_Copy_number

Protein_Interaction_PCA_filtered_final_kd$MATa_Conc <- (Protein_Interaction_PCA_filtered_final_kd$MATa_Copy_Number/(6.022*10^23))/8.2e-14
Protein_Interaction_PCA_filtered_final_kd$MATalpha_Conc <- (Protein_Interaction_PCA_filtered_final_kd$MATalpha_Copy_Number/(6.022*10^23))/8.2e-14
Protein_Interaction_PCA_filtered_final_kd$Complex_Conc <- (Protein_Interaction_PCA_filtered_final_kd$Copy_number/(6.022*10^23))/8.2e-14


Protein_Interaction_PCA_filtered_final_kd$MATa_FreeConc <- Protein_Interaction_PCA_filtered_final_kd$MATa_Conc - Protein_Interaction_PCA_filtered_final_kd$Complex_Conc
Protein_Interaction_PCA_filtered_final_kd$MATalpha_FreeConc <- Protein_Interaction_PCA_filtered_final_kd$MATalpha_Conc - Protein_Interaction_PCA_filtered_final_kd$Complex_Conc

Protein_Interaction_PCA_filtered_final_kd$Kd <- Protein_Interaction_PCA_filtered_final_kd$MATa_FreeConc*Protein_Interaction_PCA_filtered_final_kd$MATalpha_FreeConc/Protein_Interaction_PCA_filtered_final_kd$Complex_Conc

Protein_Interaction_PCA_filtered_final_kd <- 
  Protein_Interaction_PCA_filtered_final_kd %>%
  dplyr::filter((MATa_FreeConc > 0) & (MATalpha_FreeConc > 0)) %>%
  dplyr::filter(Kd > 0)

row.names(Protein_Interaction_PCA_filtered_final_kd) <- NULL

write.csv(Protein_Interaction_PCA_filtered_final_kd, file = "../results/Kd_and_GI/Estimated_Kd2_PCA_intensitiesNorm_LinFit_MeanCopyNumberUsed.csv")

