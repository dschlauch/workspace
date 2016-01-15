library(bereR)
library(pandaR)
library(bptools)
library(reshape2)
library(penalized)
library(Biobase)
library(org.Hs.eg.db)
library(foreach)
library(doParallel)
library(limma)
library(igraph)

setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2016/ECLIPSE_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2016/COPDGene_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2016/LGRC_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2016/LTCOPD_combined_runs/')

setwd('~/gd/Harvard/Research/TM_outputs/CISPB/BERE/ECLIPSE_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/CISPB/BERE/COPDGene_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/CISPB/BERE/LGRC_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/CISPB/BERE/LTCOPD_combined_runs/')

fileList <- sort(list.files('./tms/', full.names=T))
#remove observed except 1
observedIndex <- grep("_1.rds",fileList)[1]
nullIndices <- which(!grepl("_1.rds",fileList))
transMatrices <- lapply(fileList[c(observedIndex,nullIndices)], readRDS)

load("./activeImage.RData")
analysisCode <-""

outputDir <- '.'
source('~/gd/Harvard/Research/R_workspace/process_TM.R')

source('~/gd/Harvard/Research/R_workspace/consolidateResultTables.R')
