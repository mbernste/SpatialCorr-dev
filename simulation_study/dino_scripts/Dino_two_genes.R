#!/usr/bin/env Rscript

##################
## Install Dino ##
##################

# devtools::install_github("JBrownBiostat/Dino", build_vignettes = FALSE)


###############
## Libraries ##
###############
library(Dino)


##########
## Data ##
##########
## Script requires 2 shell arguments:
## 1) location of the file to normalize
## 2) directory including file name in which to save the normalized data
args <- commandArgs(trailingOnly=TRUE)

## source .tsv location ##
rawDat <- read.table(args[1], sep = "\t", header = TRUE)

## save directory ##
saveDir <- args[2]


###############
## Normalize ##
###############
sizeVec <- log(rawDat$size_factors)
sizeVec <- sizeVec - median(sizeVec)
normDat <- rbind(t(rawDat[, c("count_gene_0", "count_gene_1")]), 1)
colnames(normDat) <- rawDat$X
rownames(normDat) <- c("Gene_1", "Gene_2", "Null_Gene")
normDat <- Dino(normDat, depth = sizeVec, slope = 1, nCores = 1)

rawDat$dino_gene_0 <- normDat[1, ]
rawDat$dino_gene_1 <- normDat[2, ]

##########
## Save ##
##########
write.table(rawDat, file = saveDir, row.names=FALSE, sep="\t", quote=FALSE)

