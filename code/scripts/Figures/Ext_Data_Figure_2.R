# Xavier Castellanos-Girouard
#
#
# Take Stoichiometry PPI Network and assign proteins to complexes
#
# Date First Created: February 15 2024
# Date Last Modified: May 16 2024

#### Import Libraries ####

library(dplyr)
library(tidyr)
library(org.Sc.sgd.db)
library(ggplot2)

#### Import data ####


ComplexPortal <-
  read.csv("../data/Modules/ComplexPortal_559292.tsv",
           sep = "\t")

Yeast_PPI_Stoich_DF <-
  read.csv("../results/Stoichiometry_and_GI/Yeast_Interactome_Stoich.csv", row.names = 1)

#### Convert Uniprot IDs to systematic ORFs ####


### Convert Uniprot IDs to Systematic ORFs
bimap <- org.Sc.sgdUNIPROT # Get Uniprot to sys ORF bimap

# Get the Systematic ORF IDs that are mapped to a Uniprot ID
mapped_genes <- mappedkeys(bimap)

# Convert map to list. ORFs as entry names, Uniprot IDs as entry values
list_map <- as.list(bimap[mapped_genes])

# Some ORFs have multiple Uniprot IDs, collapse them into a single string
list_map <- sapply(list_map, FUN = paste0, collapse = ",")

# Construct a conversion dataframe
Uniprot_2_ORF_DF <- 
  data.frame(ORF_ID = names(list_map),
             Uniprot_ID = unlist(unname(list_map)))

# Separate ORF that have more than one Uniprot ID into new rows
Uniprot_2_ORF_DF <-
  Uniprot_2_ORF_DF %>%
  tidyr::separate_longer_delim(cols = "Uniprot_ID",
                               delim = ",")

#### ComplexPortal: uniprot to systemic ORF ####

## Extract only uniprot codes, make available as character vector

# Remove stoichiometry data
ComplexPortal$Expanded.participant.list_formatted <- 
  gsub(ComplexPortal$Expanded.participant.list, 
       pattern = "\\(.\\)", 
       replacement = "")


# Split complex elements, i.e. string -> character vector where every element is
# a protein.
ComplexPortal$Expanded.participant.list_formatted <- 
  strsplit(ComplexPortal$Expanded.participant.list_formatted, 
           split="\\|")

# This function takes a vector of Uniprot IDs and converts them into a vector
# of Locus IDs
Uniprot2Locus <-
  function(participants_list_Uniprot){
    
    print(participants_list_Uniprot)
    
    # Output of the match function is a vector of same length as
    # participants_list_Uniprot. Each element is an index of Uniprot_DF
    # where participant matched Entry.
    UniprotID_matchs <- match(participants_list_Uniprot, Uniprot_2_ORF_DF$Uniprot_ID)
    
    LocusIDs <- Uniprot_2_ORF_DF$ORF_ID[UniprotID_matchs]
    
    # Some Uniprot IDs match to multiple loci, currently in a single merged
    # string. The following command separates the multiple loci into 
    # individual strings.
    
    LocusIDs <- unlist(strsplit(LocusIDs, split = ";"))
    
    return(LocusIDs)
  }

Locus2Uniprot <-
  function(participants_list_Locus){
    
    #print(participants_list_Locus)
    
    # Output of the match function is a vector of same length as
    # participants_list_Uniprot. Each element is an index of Uniprot_DF
    # where participant matched Entry.
    UniprotID_matchs <- 
      match(participants_list_Locus, 
            Uniprot_2_ORF_DF$ORF_ID)
    
    Uniprot_IDs <- Uniprot_2_ORF_DF$Uniprot_ID[UniprotID_matchs]
    
    # Some Uniprot IDs match to multiple loci, currently in a single merged
    # string. The following command separates the multiple loci into 
    # individual strings.
    
    #print(Uniprot_IDs)
    
    Uniprot_IDs <- unlist(strsplit(Uniprot_IDs, split = ";"))
    
    return(Uniprot_IDs)
  }

ComplexPortal$Expanded.participant.list_locus <-
  sapply(ComplexPortal$Expanded.participant.list_formatted,
         FUN = Uniprot2Locus)


#### Preclassification of ORFs into modules or connectors ####

# list of loci that are in complexes
Complex_loci <- c(unlist(ComplexPortal$Expanded.participant.list_locus))

# Remove Na values
Complex_loci <- Complex_loci[!is.na(Complex_loci)]

# Split instances of multiple loci in the same string
Complex_loci <- unlist(strsplit(Complex_loci, split = ";"))

# Reduce list to unique loci (remove duplicates)
Complex_loci <- unique(Complex_loci)

# Check which interactorA is in a module (complex)
MS_PPI_target_mod <- Yeast_PPI_Stoich_DF$target %in% Complex_loci

# Check which interactorB is in a module (complex)
MS_PPI_source_mod <- Yeast_PPI_Stoich_DF$source %in% Complex_loci


#### Functions related to Inter-module Network ####

## Map complex to ORF (keeping track of ORF origin)
# Note: To distinguish between inter and intra complex interactions, the complex
# of both ORFs of an interaction must be determined.
map_orf2_complex <-
  function(ORF, complex2orf_table){
    
    # Get boolean vector that checks whether the current ORF is in a complex
    # Vector will indicate what positions have a complex that contains an ORF
    Orf_inComplex <-
      sapply(complex2orf_table$Expanded.participant.list_locus,
             FUN = function(x, ORF_locus){return((ORF_locus %in% unlist(x)))},
             ORF_locus = ORF)
    
    
    Complexes_containing_ORF <- complex2orf_table[Orf_inComplex,]$Complex.ac
    
    return(Complexes_containing_ORF)
  }



#### Inter-module network : Get interactions between Complexes ####

# Store module-module (complex-complex) interactions (Intra and inter)
MS_PPI_module_module <-
  Yeast_PPI_Stoich_DF[c(MS_PPI_source_mod & MS_PPI_target_mod),]

# Find complexes for InteractorA
MS_PPI_module_module$source_complexes <- 
  sapply(MS_PPI_module_module$source,
         FUN = map_orf2_complex, 
         complex2orf_table = ComplexPortal)

# Find complexes for InteractorB
MS_PPI_module_module$target_complexes <- 
  sapply(MS_PPI_module_module$target,
         FUN = map_orf2_complex, 
         complex2orf_table = ComplexPortal)

## Get intercomplex Interactions

# Output is boolean vector
# Note: an assumption made here is that, in a PPI, if InteractorA and 
# InteractorB have a complex in common, the experimental detection of the PPI 
# was most likely made in the context of the complex. Therefore, if InteractorA
# and InteractorB have any complex in common, it is an intra-complex interaction.
# This assumption may be false in some instances, but they they cannot be 
# distinguished.
find_InterComplex_PPI <- 
  function(index, PPI_table, complex_table){
    InteractorA_Complex <- PPI_table$source_complexes[[index]]
    InteractorB_Complex <- PPI_table$target_complexes[[index]]
    
    # Note: If an interactor is in multiple complexes, datatype is vector.
    #   Otherwise (interactor is only in one complex), datatype is string.
    #   To make appropriate comparisons (avoid error), interactor complex
    #   info must be treated differently according to their datatype.
    
    # Find cases where InteractorA is only in a single complex.
    IntA_Single <- (length(InteractorA_Complex) == 1)
    
    # Find cases where InteractorB is only in a single complex.
    IntB_Single <- (length(InteractorB_Complex) == 1)
    
    # Case where InteractorA and InteractorB both are only in one complex
    if(IntA_Single & IntB_Single){
      
      # Check if InteractorA is in the same complex as InteractorB
      compared_complexes <- InteractorA_Complex == InteractorB_Complex
      
      return(!compared_complexes)}
    
    # Case where InteractorA is in one complex, InteractorB is in multiple
    if(IntA_Single & !IntB_Single){
      
      # Check if InteractorA is in any of the same complexes as interactorB
      compared_complexes <- InteractorA_Complex %in% InteractorB_Complex
      
      return(!compared_complexes)}
    
    # Case where InteractorB is in one complex, InteractorA is in multiple
    if(!IntA_Single & IntB_Single){
      
      # Check if InteractorB is in any of the same complexes as interactorA
      compared_complexes <- InteractorB_Complex %in% InteractorA_Complex
      
      return(!compared_complexes)}
    
    
    # Case where both interactors are in multiple complexes
    if(!IntA_Single & !IntB_Single){
      
      # Check if Interactors have any complex in common
      compared_complexes <- (length(base::intersect(InteractorB_Complex, InteractorA_Complex)) >= 1)
      
      return(!compared_complexes)}
    
    print("SLIPPERY CASE, ERROR")
  }



## Get Uniprot IDs for each interactor locus
MS_PPI_module_module$Uniprot_source <- 
  sapply(MS_PPI_module_module$source, 
         FUN = Locus2Uniprot)

MS_PPI_module_module$Uniprot_target <- 
  sapply(MS_PPI_module_module$target, 
         FUN = Locus2Uniprot)


# Apply function to get boolean vector indicating intermodule interactions
intermodule_PPI_bool <- 
  sapply(seq_along(MS_PPI_module_module$Uniprot_source), 
         FUN = find_InterComplex_PPI,
         PPI_table = MS_PPI_module_module,
         complex_table = ComplexPortal)

# Subset with boolean vector to get intermodule interactions
MS_PPI_intermodule <- MS_PPI_module_module[intermodule_PPI_bool,]

row.names(MS_PPI_intermodule) <- NULL

# Only keep useful columns
MS_PPI_intermodule <- 
  MS_PPI_intermodule[,c("source",
                        "target",
                        "source_complexes",
                        "target_complexes")]

colnames(MS_PPI_intermodule) <- 
  c("source_locus", 
    "target_locus",
    "source_Complex", 
    "target_Complex")


#### IntAct (intra-module PPIs) : Get interactions within Complexes ####

## Get intracomplex Interactions

# Note: an assumption made here is that, in a PPI, if InteractorA and 
# InteractorB have a complex in common, the experimental detection of the PPI 
# was most likely made in the context of the complex. Therefore, if InteractorA
# and InteractorB have any complex in common, it is an intra-complex interaction.
# This assumption may be false in some instances, but they they cannot be 
# distinguished.

# Output is boolean vector. TRUE if in the same complex, FALSE otherwise
find_IntraComplex_PPI <- 
  function(index, PPI_table, complex_table){
    InteractorA_Complex <- PPI_table$source_complexes[[index]]
    InteractorB_Complex <- PPI_table$target_complexes[[index]]
    
    # Note: If an interactor is in multiple complexes, datatype is vector.
    #   Otherwise (interactor is only in one complex), datatype is string.
    #   To make appropriate comparisons (avoid error), interactor complex
    #   info must be treated differently according to their datatype.
    
    # Find cases where InteractorA is only in a single complex.
    IntA_Single <- (length(InteractorA_Complex) == 1)
    
    # Find cases where InteractorB is only in a single complex.
    IntB_Single <- (length(InteractorB_Complex) == 1)
    
    # Case where InteractorA and InteractorB both are only in one complex
    if(IntA_Single & IntB_Single){
      
      # Check if InteractorA is in the same complex as InteractorB
      compared_complexes <- InteractorA_Complex == InteractorB_Complex
      
      return(compared_complexes)}
    
    # Case where InteractorA is in one complex, InteractorB is in multiple
    if(IntA_Single & !IntB_Single){
      
      # Check if InteractorA is in any of the same complexes as interactorB
      compared_complexes <- InteractorA_Complex %in% InteractorB_Complex
      
      return(compared_complexes)}
    
    # Case where InteractorB is in one complex, InteractorA is in multiple
    if(!IntA_Single & IntB_Single){
      
      # Check if InteractorB is in any of the same complexes as interactorA
      compared_complexes <- InteractorB_Complex %in% InteractorA_Complex
      
      return(compared_complexes)}
    
    
    # Case where both interactors are in multiple complexes
    if(!IntA_Single & !IntB_Single){
      
      # Check if Interactors have any complex in common
      compared_complexes <- (length(base::intersect(InteractorB_Complex, InteractorA_Complex)) >= 1)
      
      return(compared_complexes)}
    
    print("SLIPPERY CASE, ERROR")
  }

# Apply function to get boolean vector indicating intramodule interactions
intramodule_PPI_bool <- 
  sapply(seq_along(MS_PPI_module_module$Uniprot_source), 
         FUN = find_IntraComplex_PPI,
         PPI_table = MS_PPI_module_module)

# Subset with boolean vector to get intramodule interactions
MS_PPI_intramodule <- MS_PPI_module_module[intramodule_PPI_bool,]

row.names(MS_PPI_intramodule) <- NULL

# Only keep useful columns
MS_PPI_intramodule <- 
  MS_PPI_intramodule[,c("source",
                        "target",
                        "source_complexes",
                        "target_complexes")]

colnames(MS_PPI_intramodule) <- 
  c("source_locus", 
    "target_locus",
    "source_Complex", 
    "target_Complex")

intramodule_PPI <- MS_PPI_intramodule


# find module-connector interactions
MS_PPI_module_connector <- 
  Yeast_PPI_Stoich_DF[xor((!MS_PPI_source_mod), (!MS_PPI_target_mod)),]

# find connector-connector interactions
MS_PPI_connector_connector <-
  Yeast_PPI_Stoich_DF[(!MS_PPI_source_mod & !MS_PPI_target_mod),]

# Remove self-interactions
MS_PPI_connector_connector <- 
  MS_PPI_connector_connector[!(MS_PPI_connector_connector$source == MS_PPI_connector_connector$target),]

## Make list of module genes, and list of connector genes

# list of connector proteins
#note: !IntAct_InteractorX_mod means not in module, so by default is a connector
source_connector <- Yeast_PPI_Stoich_DF$source[!MS_PPI_source_mod]
target_connector <- Yeast_PPI_Stoich_DF$target[!MS_PPI_target_mod]
connector_proteins <- unique(source_connector, target_connector)

## Make full complex-connector network
# Make network where all proteins of a complex are merged into one node.
# Proteins not in a complex remain as an individual node.
is_inComplex <-
  function(index, module_connector_DF){
    interaction <- module_connector_DF[index,]
    
    A_inComplex <- interaction$source %in% Complex_loci # is A in a complex? (ie a module)
    B_inComplex <- interaction$target %in% Complex_loci # is B in a complex? (ie a module)
    
    if (A_inComplex){ # If A in a complex:
      module <- interaction$source # then a is module protein
      connector <- interaction$target} # and b is connector protein
    else { # Otherwise:
      module <- interaction$target # B is module protein
      connector <- interaction$source # A is connector protein
    }
    
    print(paste(module, connector))
    
    in_complex <- # Check what complexes the module protein is in (indexes)
      sapply(ComplexPortal$Expanded.participant.list_locus,
             function(complex_loci, module_prot){
               return(module %in% complex_loci)
             },
             module_prot = module)
    
    # Get identity of the complexes using indexes above.
    Interacting_Complexes <- ComplexPortal$Complex.ac[in_complex]
    
    return(list(module, connector, Interacting_Complexes, "None"))
  }

# get connector-complex interactions; they will be in list of dataframes
connector_complex_DFs <- 
  lapply(seq_along(MS_PPI_module_connector$source),
         FUN = is_inComplex,
         module_connector_DF = MS_PPI_module_connector)

connector_complex_network <- 
  data.frame(source_locus = character(),
             target_locus = character(),
             source_Complex = character(),
             target_Complex = character())

# Concatenate the lists into one dataframe
connector_complex_network <- do.call(rbind, connector_complex_DFs)

# set column names
colnames(connector_complex_network) <- c("source_locus",
                                         "target_locus",
                                         "source_Complex",
                                         "target_Complex")

connector_complex_network <- as.data.frame(connector_complex_network)

## Format connector-connector network
connector_connector_network <- 
  data.frame(source_locus = MS_PPI_connector_connector$source, 
             target_locus = MS_PPI_connector_connector$target, 
             source_Complex = rep("None",nrow(MS_PPI_connector_connector)),
             target_Complex = rep("None",nrow(MS_PPI_connector_connector)))

# For downstream compatibility
connector_complex_network$source_locus <- as.vector(connector_complex_network$source_locus)
connector_complex_network$target_locus <- as.vector(connector_complex_network$target_locus)

#### Make Full Classified Network ####

intramodule_PPI$Interaction_Type <- "IntraModule"
IntraModule_network <- intramodule_PPI

MS_PPI_intermodule$Interaction_Type <- "InterModule"
InterModule_network <- MS_PPI_intermodule

connector_complex_network$Interaction_Type <- "ModuleConnector"
ModuleConnector_network <- connector_complex_network

connector_connector_network$Interaction_Type <- "ConnectorConnector"
ConnectorConnector_network <- connector_connector_network

#View(IntraModule_network)
#View(InterModule_network)
#View(ModuleConnector_network)


Complex_PPI_network <-
  rbind(IntraModule_network,
        InterModule_network,
        ModuleConnector_network,
        ConnectorConnector_network)

row.names(Complex_PPI_network) <- NULL


# Convert any Vector of Complexes into a string with each complex ID separated
# by the "|" character.
Complex_PPI_network$source_Complex <- 
  sapply(Complex_PPI_network$source_Complex,
         FUN = function(x){return(do.call(paste, c(as.list(x), sep = "|")))}
  )

Complex_PPI_network$target_Complex <-
  sapply(Complex_PPI_network$target_Complex,
         FUN = function(x){return(do.call(paste, c(as.list(x), sep = "|")))}
  )

#View(Complex_PPI_network)

Complex_PPI_network$source_locus <- as.vector(unlist(Complex_PPI_network$source_locus))
Complex_PPI_network$source_Complex <- as.vector(unlist(Complex_PPI_network$source_Complex))
Complex_PPI_network$target_locus <- as.vector(unlist(Complex_PPI_network$target_locus))


#### Export ####


## Add Binding values
Complex_PPI_network_bind1 <- 
  merge(x = Yeast_PPI_Stoich_DF,
        y = Complex_PPI_network,
        by.x = c("source", "target"),
        by.y = c("source_locus", "target_locus"),
        all.x = FALSE,
        all.y = FALSE)

## Add binding values
Complex_PPI_network_bind2 <- 
  merge(x = Yeast_PPI_Stoich_DF,
        y = Complex_PPI_network,
        by.x = c("source", "target"),
        by.y = c("target_locus", "source_locus"),
        all.x = FALSE,
        all.y = FALSE)

Complex_PPI_network_bind <-
  rbind(Complex_PPI_network_bind1, Complex_PPI_network_bind2)

#write.csv(Complex_PPI_network_bind, "../data/results/Yeast_Interactome_Stoich_withComplex.csv")

#### Graphs ####


Complex_PPI_network_bind_small_1 <- Complex_PPI_network_bind
Complex_PPI_network_bind_small_1$Label <- "All Interactions"

Complex_PPI_network_bind_small_2 <- Complex_PPI_network_bind[Complex_PPI_network_bind$Interaction_Type == "IntraModule",]
Complex_PPI_network_bind_small_2$Label <- "Interactions Within Modules"

Complex_PPI_network_bind_small <- rbind(Complex_PPI_network_bind_small_1, Complex_PPI_network_bind_small_2)

CPX_Interaction_density_p <-
  ggplot() +
  #geom_point(data = Complex_PPI_network_bind_small,
  #           mapping = aes(x = Interaction_Stoichiometry,
  #                         y = Abundance_Stoichiometry)) +
  geom_density_2d(data = Complex_PPI_network_bind_small_1,
                 mapping = aes(x = Interaction_Stoichiometry,
                               y = Abundance_Stoichiometry,
                               color = 'All Interactions'),
                 bins = 10,
                 linewidth = 1) +
  geom_density_2d(data = Complex_PPI_network_bind_small_2,
                  mapping = aes(x = Interaction_Stoichiometry,
                                y = Abundance_Stoichiometry,
                                color = '     Interactions\nWithin Complexes'),
                  bins = 10,
                  linewidth = 1) +
  geom_hline(yintercept = 1, linetype='dashed', color = "grey60") +
  geom_vline(xintercept = 1, linetype='dashed', color = "grey60") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black"),
        axis.text.x.bottom = element_text(color = "black", size = 12),
        axis.text.y.left = element_text(color = "black", size = 12),
        axis.title.x.bottom = element_text(color = "black", size = 14),
        axis.title.y.left = element_text(color = "black", size = 14),
        plot.title = element_text(hjust = 0.5,
                                  color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 1),
        legend.text=element_text(size=12)) +
  xlab("Interaction Stoichiometry") +
  ylab("Abundance Stoichiometry") +
  #ggtitle("Yeast Stoichiometry Map") +
  scale_x_log10(breaks = c(10^-3, 10^-2, 10^-1, 10^0, 10^1),
                labels = c(expression(10^"-3"), expression(10^"-2"), 
                           expression(10^"-1"), expression(10^0),
                           expression(10^1))) +
  scale_y_log10(breaks = c(10^-1.5, 10^-1, 10^-0.5, 10^0, 10^0.5, 10^1, 10^1.5),
                labels = c(expression(10^"-1.5"), expression(10^"-1"), 
                           expression(10^"-0.5"), expression(10^0), 
                           expression(10^0.5), expression(10^1),
                           expression(10^1.5))) +
  scale_color_manual('', values = c("red", "black"))


png("../results/Figures/Extended_Data_Figures/SupplementaryFigure2.png", units="mm", width=160, height=80, res=300)
CPX_Interaction_density_p
dev.off()
