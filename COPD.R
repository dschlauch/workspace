COPDDir <- "~/gd/Harvard/Research/data/"
exprData <- t(read.table(file.path(COPDDir,"./DREAM5_NetworkInferenceChallenge_Data/Network1_Data/Network1_expression_data.tsv"), header=T,stringsAsFactors=F))
chipFeatures <- read.table(file.path(COPDDir,"./DREAM5_NetworkInferenceChallenge_Data/Network1_Data/Network1_chip_features.tsv"),stringsAsFactors=F)
transFactors <- read.table(file.path(COPDDir,"./DREAM5_NetworkInferenceChallenge_Data/Network1_Data/Network1_transcription_factors.tsv"),stringsAsFactors=F)
goldStandard <- read.table(file.path(COPDDir,"./DREAM5_NetworkInference_GoldStandard/DREAM5_NetworkInference_GoldStandard_Network1.tsv"),stringsAsFactors=F)

motifs     <- read.table("ECLIPSE_Blood_Motif.txt",header=F)
exprData      <- read.table("ECLIPSE_Blood_Exp.txt",row.names=1,header=T)
chipFeatures <- read.table("ECLIPSE_blood.txt",header=T,fill = TRUE, sep="\t",row.names=1)
goldStandard <- NA