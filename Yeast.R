yeastDir <- "~/gd/Harvard/Research/data/YeastData/"
data(yeast)
exprData <- read.table(file.path(yeastDir,"YeastData_KOExpression.txt"),stringsAsFactors=F,row.names=1)
chipFeatures <- NA
transFactors <- as.character(unique(yeast$motif[,1]))
goldStandard <- read.table(file.path(yeastDir,"YeastData_ChIP.txt"),stringsAsFactors=F)
motifs <- read.table(file.path(yeastDir,"YeastData_Motif.txt"),stringsAsFactors=F)
motifs <- merge(goldStandard,yeast$motif,by=c("V1","V2"), all = TRUE)[,-3]
motifs[is.na(motifs)] <- 0
