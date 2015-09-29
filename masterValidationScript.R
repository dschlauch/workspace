setwd('~/gd/Harvard/Research/R_workspace')
source("./validation.R")
datasets <- c("Yeast", "DREAM5a", "DREAM5c", "DREAM5d")
sapply(datasets, validateMethodsOnDataset)
