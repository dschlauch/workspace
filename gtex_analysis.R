## 10/20/15
## GTEX tissue vs cell line transition matrix
library(devtools)
install_github("QuackenbushLab/pandaR", username = getOption("github.user"), ref = "master")
library(pandaR)
library(bereR)
library(bptools)
library(reshape2)
library(penalized)

gtext_filename <- "./gd/Harvard/Research/data/GTEx/gtex_cell.rdata"
motif_filename <- "./gd/Harvard/Research/data/GTEx/KG_cisbp_652.txt"
# Load the data
load(gtext_filename)

# Load motif data...
dataset <- list()
dataset$motif    <- cbind(read.table(motif_filename,header=F),1)
dataset$exp      <- exprs(both)
dataset$clinical <- pData(both)

# Data validation and preprocessing
dataset$exp      <- dataset$exp[,order(colnames(dataset$exp))]  # Make sure expression and clinical is in same order
rownames(dataset$exp) <- substring(rownames(dataset$exp),1,15)
phenotypeName <- "our_subtypes"
casesString <- "cells_ebv-transformed_lymphocytes"
controlsString <- "skin"
filterType <- NA

#Get top 20,000 variable genes
rowsds <- sort(apply(dataset$exp, 1, sd), decreasing=T)
genesIncluded <- names(rowsds[1:19000])
dataset$exp <- dataset$exp[genesIncluded,]
dataset$motif <- dataset$motif[dataset$motif[,2]%in%genesIncluded,]

# Run net inf algorithm
dataset$casesNetwork    <- bere(dataset$motif, dataset$exp[,casesFilter], verbose=T)
dataset$controlsNetwork <- bere(dataset$motif,dataset$exp[,controlsFilter], verbose=T)
# Run TM


# Plot results