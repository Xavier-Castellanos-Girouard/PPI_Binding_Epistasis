# Xavier Castellanos-Girouard

# Date First Created: March 18 2024
# Date Last Modified: May 21 2024

# Note: independant script to convert CDT files to dataframes

#### Functions ####

## Read CDT function borrowed from Kevin Coombes
readCDT <- function(filename) {
  fname <- sub('.cdt$', '', filename) # get rid of the extension
  atr <- read.table(paste(fname, 'atr', sep='.'), sep='\t',
                    header=FALSE, as.is=TRUE)
  gtr <- read.table(paste(fname, 'gtr', sep='.'), sep='\t',
                    header=FALSE, as.is=TRUE)
  cdt <- read.table(paste(fname, 'cdt', sep='.'), sep='\t',
                    header=TRUE, row.names=NULL)
  # we only need the first column of the CDT file, which contains
  # the order information, and the third column, which contains the
  # labels. The rest of the file contains the data matrix.
  rown <- as.character(cdt[,"GID"])
  coln <- colnames(cdt)
  firstRow <- 1 + which(rown=="EWEIGHT")
  firstCol <- 1 + which(coln=="GWEIGHT")
  
  gid <- as.character(cdt[,"GID"])[firstRow:nrow(cdt)]
  aid <- cdt[rown=="AID",][firstCol:ncol(cdt)]
  aid <- as.character(as.matrix(aid))
  # names all start with 'GENE' or 'NODE' (or 'ARRY') and end with 'X'
  gene.order <- 1 + as.numeric(substring(gid, 5, nchar(gid)-1))
  arry.order <- 1 + as.numeric(substring(aid, 5, nchar(aid)-1))
  # Because Cluster reorders things and because hclust and plclust wants
  # to do the same, we have to reinvert the ordering during passage from
  # one to the other
  gene.labels <-
    as.character(cdt$NAME)[firstRow:nrow(cdt)][order(gene.order)]
  arry.labels <- coln[firstCol:ncol(cdt)][order(arry.order)]
  
  temp <- as.matrix(cdt[firstRow:nrow(cdt), firstCol:ncol(cdt)])
  temp <- temp[order(gene.order), order(arry.order)]
  data <- matrix(as.numeric(temp), ncol=ncol(temp))
  dimnames(data) <- list(gene.labels, arry.labels)
  
  # The gtr file contains the "distances" in column 4. Actually,
  # Eisen's Cluster program reports similarities instead of
  # distances. This fix assumes that some kind of correlation was
  # the measure of similarity....
  gene.height <- 1 - gtr$V4
  arry.height <- 1 - atr$V4
  # Columns 2 and 3 describe the two branches below each node.
  # Nodes are listed from bottom to top since clustering is
  # agglomerative.
  #
  foo <- function(alt) {
    # Again, we get the numeric part of the label
    base <- as.numeric(substring(alt, 5, nchar(alt)-1))
    # We also need to know whether the label is a "GENE" or a "NODE".
    # The 'hclust' objects use negative integers to indicate nodes.
    type1 <- rep(1, length(base))
    type1[substring(alt, 1, 4) %in% c('GENE', "ARRY")] <- -1
    base <- base*type1  # make nodes negative
    adder <- (type1-1)/2  # offset the negatives to change from starting
    # at 1 to starting at 0.
    base + adder
  }
  gene.merge1 <- foo(gtr$V3)
  gene.merge2 <- foo(gtr$V2)
  arry.merge1 <- foo(atr$V3)
  arry.merge2 <- foo(atr$V2)
  
  # put everything together into a list and make it an hclust object
  gene <- list(merge=as.matrix(cbind(gene.merge1, gene.merge2)),
               height=gene.height,
               order=gene.order,
               labels=gene.labels,
               method='modified centroid',
               call=NULL,
               dist.method='Pearson correlation')
  class(gene) <- 'hclust'
  arry <- list(merge=as.matrix(cbind(arry.merge1, arry.merge2)),
               height=arry.height,
               order=arry.order,
               labels=arry.labels,
               method='modified centroid',
               call=NULL,
               dist.method='Pearson correlation')
  class(arry) <- 'hclust'
  list(gene=gene,arry=arry, data=data)
}


#### Main ####
CDTlist <- readCDT("../results/IntermediateFiles/Kd_inferred_GI.cdt")

## Convert CDT structure into dataframe
CDT_DF <- CDTlist$data
row.names(CDT_DF) <- CDTlist$gene$labels
colnames(CDT_DF) <- CDTlist$arry$labels # Note!: Here "-" in ORF names were replaced with "."

# Reorder to match clustering
CDT_DF <- CDT_DF[CDTlist$gene$order, CDTlist$arry$order]


CDT_DF <- as.data.frame(CDT_DF)
#CDT_DF$ORF1 <- row.names(CDT_DF)


write.csv(CDT_DF, "../results/Kd_and_GI/Clustered_Kd_inferred_GI.csv")
