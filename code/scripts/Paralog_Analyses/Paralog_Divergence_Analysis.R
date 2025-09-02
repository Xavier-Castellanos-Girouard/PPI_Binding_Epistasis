# Xavier Castellanos-Girouard

# Perform analysis related to paralog divergence

# Date First Created: March 14 2024 
# Date Last Modified: May 20 2024


#### Import Libraries ####

library(dplyr)
library(tidyr)
library(seqinr)

#### Import data ####

## Import duplication types of the paralogs
paralog_DF <- read.csv("../results/IntermediateFiles/Paralog_DuplicationType.csv", row.names = 1)

## Import paralogs merged with GI and Kd data
MS_Kd_GI_AllPara_DF <- read.csv("../results/Paralog_Analyses/Kd_GI_ParalogType.csv")


#### Export nucleotide sequences of paralogs for reverse align input ####

### Read individual DNA sequences, and write fasta file containing both paralog sequences
ORF_Nucl_sequences <- 
  read.fasta(file = "../data/Paralogs/orf_coding_all.fasta",
             seqtype = c("DNA"))

ORF_Nucl_sequences <- 
  lapply(ORF_Nucl_sequences,
         FUN = paste0,
         collapse = "")

for (i in 1:nrow(paralog_DF)){
  paralog1 <- paralog_DF$syst1[i]
  paralog2 <- paralog_DF$syst2[i]
  
  paralog1_seq <- ORF_Nucl_sequences[paralog1]
  paralog2_seq <- ORF_Nucl_sequences[paralog2]
  
  #print(paste(paralog1, paralog2))
  
  seqinr::write.fasta(sequences = list(paralog1_seq, paralog2_seq),
                      as.string = TRUE,
                      names = c(paralog1, paralog2),
                      file.out = paste0("../results/IntermediateFiles/reverseAlign_nucl_FastaFiles/", 
                                        paralog1,
                                        "_",
                                        paralog2,
                                        ".fasta"))
}

#### Align nucleic acid sequences using amino acid sequences ####

bug <- c("YDL061C", "YLR388W", "YLR466C-B", "YLL066W-B")


for (i in 1:nrow(paralog_DF)){
  paralog1 <- paralog_DF$syst1[i]
  paralog2 <- paralog_DF$syst2[i]

  #print(paste(paralog1, paralog2))
  
  # Note: Some paralogs raise an error, skip these
  if ((paralog1 %in% bug) & (paralog2 %in% bug)){
    next
  }
  ## Nucleotide-level alignments
  input_dir <- paste0("../results/IntermediateFiles/reverseAlign_nucl_FastaFiles/",
                      paralog1, 
                      "_",
                      paralog2,
                      ".fasta")
  ## Amino-acid-level alignemnts
  prot_aln_dir <- paste0("../results/IntermediateFiles/ClustalW_Output/",
                         paralog1,
                         "_",
                         paralog2,
                         ".aln")
  
  ## Designate Output path and filename
  output_dir <- paste0("../results/IntermediateFiles/reverseAlign_output_Files/",
                       paralog1,
                       "_",
                       paralog2,
                       ".fasta")
  
  ## Execute reverse.align to compute kaks
  seqinr::reverse.align(nucl.file = input_dir,
                        protaln.file = prot_aln_dir,
                        input.format = "clustal",
                        out.file = output_dir)
}



## Import alignments

Nucl_align_list <- list()

for (i in 1:nrow(paralog_DF)){
  paralog1 <- paralog_DF$syst1[i]
  paralog2 <- paralog_DF$syst2[i]
  
  #print(paste(paralog1, paralog2))
  
  # Note: Some paralogs cause a bug, skip these
  if ((paralog1 %in% bug) & (paralog2 %in% bug)){
    next
  }
  
  alignment <- 
    read.alignment(file = paste0("../results/IntermediateFiles/reverseAlign_output_Files/",
                                 paralog1,
                                 "_",
                                 paralog2,
                                 ".fasta"),
                   format = "fasta")
  
  Nucl_align_list[paste0(paralog1, "_", paralog2)] <- list(alignment)
}

### Calculate Ka and Ks (synonymous and non-synonymous substitutions)

kaks_list <- list()

for (paralogs in names(Nucl_align_list)){
  kaks_computation <- kaks(Nucl_align_list[[paralogs]])
  
  kaks_list[paralogs] <- list(kaks_computation)
}


## Make a DF from Ka data
kaks_paralog_names <- names(kaks_list)

paralog_kaks_DF <- do.call(rbind, kaks_list)

paralog_kaks_DF <- as.data.frame(paralog_kaks_DF)

rownames(paralog_kaks_DF) <- NULL

paralog_kaks_DF$paralogs <- kaks_paralog_names
paralog_kaks_DF$ka <- as.numeric(paralog_kaks_DF$ka)
paralog_kaks_DF$ks <- as.numeric(paralog_kaks_DF$ks)
paralog_kaks_DF$vka <- as.numeric(paralog_kaks_DF$vka)
paralog_kaks_DF$vks <- as.numeric(paralog_kaks_DF$vks)

write.csv(paralog_kaks_DF, "../results/Paralog_Analyses/all_paralog_kaks.csv")


#### Assign Sequence Divergences to Kd/GI dataset ####

## Assign sequence divergence to paralogs

MS_Kd_GI_AllPara_DF$source_Ka <- NA
MS_Kd_GI_AllPara_DF$target_Ka <- NA

paralog_list <- unlist(strsplit(paste0(names(kaks_list), collapse = "_"), split = "_"))

for (i in 1:nrow(MS_Kd_GI_AllPara_DF)){
  
  if (MS_Kd_GI_AllPara_DF$Type[i] == "NoParalog"){
    next
  }
  
  gene1 <- MS_Kd_GI_AllPara_DF$source[i]
  gene2 <- MS_Kd_GI_AllPara_DF$target[i]
  
  #print(paste(gene1, gene2))
  
  if ((gene1 %in% bug) | (gene2 %in% bug)){
    next
  }
  
  gene1_paralogs <- names(kaks_list)[which(grepl(gene1, names(kaks_list)) == TRUE)]
  gene2_paralogs <- names(kaks_list)[which(grepl(gene2, names(kaks_list)) == TRUE)]
  
  
  if((length(gene1_paralogs)>0) & (length(gene2_paralogs)>0)){ ## Paralogous heteromers; IntParalog
    if ((gene1_paralogs == gene2_paralogs)){
      MS_Kd_GI_AllPara_DF$source_Ka[i] <- kaks_list[[gene1_paralogs]]$ka
      MS_Kd_GI_AllPara_DF$target_Ka[i] <- kaks_list[[gene2_paralogs]]$ka
    }
  }
  
  if (MS_Kd_GI_AllPara_DF$Type[i] == "OneParalog"){
    if (gene1 %in% paralog_list){
      MS_Kd_GI_AllPara_DF$source_Ka[i] <- kaks_list[[gene1_paralogs]]$ka
      }
    if (gene2 %in% paralog_list){
      MS_Kd_GI_AllPara_DF$target_Ka[i] <- kaks_list[[gene2_paralogs]]$ka
    }
  }
  
  if (MS_Kd_GI_AllPara_DF$Type[i] == "TwoParalog"){
    MS_Kd_GI_AllPara_DF$source_Ka[i] <- kaks_list[[gene1_paralogs]]$ka
    MS_Kd_GI_AllPara_DF$target_Ka[i] <- kaks_list[[gene2_paralogs]]$ka
  }
}


### Export
write.csv(MS_Kd_GI_AllPara_DF, "../results/Paralog_Analyses/Kd_GI_Divergence.csv")
