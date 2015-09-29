dreamDir <- "~/gd/Harvard/Research/data/"
exprData <- t(read.table(file.path(dreamDir,"./DREAM5_NetworkInferenceChallenge_Data/Network1_Data/Network1_expression_data.tsv"), header=T,stringsAsFactors=F))
chipFeatures <- read.table(file.path(dreamDir,"./DREAM5_NetworkInferenceChallenge_Data/Network1_Data/Network1_chip_features.tsv"),stringsAsFactors=F)
transFactors <- read.table(file.path(dreamDir,"./DREAM5_NetworkInferenceChallenge_Data/Network1_Data/Network1_transcription_factors.tsv"),stringsAsFactors=F)[,1]
goldStandard <- read.table(file.path(dreamDir,"./DREAM5_NetworkInference_GoldStandard/DREAM5_NetworkInference_GoldStandard_Network1.tsv"),stringsAsFactors=F)

# Introduce error in gold standard to treat as Motifs priors
sensitivity <- .50
specificity <- .90
motifs <- goldStandard
numTrueEdges <- sum(motifs[,3])
numNonEdges <- nrow(motifs) - numTrueEdges
trueEdges <- motifs[,3]==1
nonEdges <- motifs[,3]==0
motifs[trueEdges,3] <- rbinom(numTrueEdges, 1, sensitivity)
motifs[nonEdges,3] <- rbinom(numNonEdges, 1, 1-specificity)
