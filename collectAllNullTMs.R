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

setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2014/BERE/ECLIPSE_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2014/BERE/COPDGene_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2014/BERE/LGRC_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2014/BERE/LTCOPD_combined_runs/')

setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2014/PANDA/ECLIPSE_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2014/PANDA/COPDGene_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2014/PANDA/LGRC_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/JASPAR2014/PANDA/LTCOPD_combined_runs/')

setwd('~/gd/Harvard/Research/TM_outputs/CISPB/BERE/ECLIPSE_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/CISPB/BERE/COPDGene_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/CISPB/BERE/LGRC_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/CISPB/BERE/LTCOPD_combined_runs/')

setwd('~/gd/Harvard/Research/TM_outputs/CISPB/PANDA/ECLIPSE_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/CISPB/PANDA/COPDGene_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/CISPB/PANDA/LGRC_combined_runs/')
setwd('~/gd/Harvard/Research/TM_outputs/CISPB/PANDA/LTCOPD_combined_runs/')


#####################################################
load("./activeImage.RData")
rm(all_tms)
save.image(file=file.path("activeImage.RData"))
analysisCode <-""
outputDir <- '.'

fileList <- sort(list.files('./tms/', full.names=T))
#remove observed except 1
observedIndex <- grep("_1.rds",fileList)[1]
nullIndices <- which(!grepl("_1.rds",fileList))
all_tms <- lapply(fileList[c(observedIndex,nullIndices)], readRDS)

source('~/gd/Harvard/Research/R_workspace/process_TM.R')

getwd()
rm(list = ls())
#####################################################







source('~/gd/Harvard/Research/R_workspace/consolidateResultTables.R')
