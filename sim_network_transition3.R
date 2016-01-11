# For creating gene expression data
library(WGCNA)
library(bptools)
library(gplots)
library(ROCR)
library(penalized)
library(tidyr)
library(bereR)
library(pandaR)
library(corpcor)
library(nettools)

##################################################################################################################
#  This script simulates networks from MV gaussian distributions and attempts to reconstruct them using many
#  commonly used NI methods and applies transition matrix on the results
#  Produces a number of comparison plots
#  Estimated run time of script: 30 minutes
##################################################################################################################

set.seed(123)

numTFs <- 100
numGenes <- 5000
numTransitions <- 1000
numSamples <- 500
geneNames <- paste0("Gene",1:numGenes)
TFNames <- paste0("TF",1:numTFs)
true_edges <- runif(numGenes*numTFs)*rbinom(numGenes*numTFs,1,prob=.2)
matrixA_GS <- matrix(true_edges,nrow=numTFs)
matrixB_GS <- matrixA_GS# + matrix(rnorm(10000)/10,nrow=10)
probs <- rexp(numTFs)
TFAs <- sample(1:numTFs,numTransitions, replace=T)
TFBs <- sample(1:numTFs,numTransitions, replace=T, probs)
TFBs[TFAs==TFBs] <- sample(1:numTFs,sum(TFAs==TFBs), replace=T)
sum(sum(TFAs==TFBs))
strengths <- (c(-(numTransitions/2):-1,1:(numTransitions/2)))/(numTransitions/2)
tm_true <- matrix(0,numTFs,numTFs)
  
for(i in 1:numTransitions){
    matrixB_GS[TFBs[i],] <- matrixB_GS[TFBs[i],] + matrixA_GS[TFAs[i],]*strengths[i]
    tm_true[TFAs[i],TFBs[i]] <- strengths[i]
}



varcovA <- t(matrixA_GS)%*%matrixA_GS
varcovB <- t(matrixB_GS)%*%matrixB_GS
# Increase the noise
diag(varcovA) <- 4*diag(varcovA)
diag(varcovB) <- 4*diag(varcovB)
gexpA <- t(mvrnorm(n=numSamples, mu=rep(0,numGenes), Sigma=varcovA))
gexpB <- t(mvrnorm(n=numSamples, mu=rep(0,numGenes), Sigma=varcovB))
rownames(gexpA) <- geneNames
rownames(gexpB) <- geneNames
tfExpA <- matrixA_GS %*% gexpA
tfExpB <- matrixB_GS %*% gexpB
# Add a bunch of noise to the TF expressions
tfExpA <- t(scale(t(tfExpA))) + 2*matrix(rnorm(length(c(tfExpA))),nrow=nrow(tfExpA))
tfExpB <- t(scale(t(tfExpB))) + 2*matrix(rnorm(length(c(tfExpB))),nrow=nrow(tfExpB))
rownames(tfExpA) <- TFNames
rownames(tfExpB) <- TFNames
pearsonAdjMatA <- abs(cor(t(tfExpA), t(gexpA)))
pearsonAdjMatB <- abs(cor(t(tfExpB), t(gexpB)))
wgcnaA6 <- mat2adj(cbind(t(tfExpA), t(gexpA)), method="WGCNA")[1:numTFs,-(1:numTFs)]
wgcnaB6 <- mat2adj(cbind(t(tfExpB), t(gexpB)), method="WGCNA")[1:numTFs,-(1:numTFs)]
wgcnaA12 <- mat2adj(cbind(t(tfExpA), t(gexpA)), method="WGCNA", P=12)[1:numTFs,-(1:numTFs)]
wgcnaB12 <- mat2adj(cbind(t(tfExpB), t(gexpB)), method="WGCNA", P=12)[1:numTFs,-(1:numTFs)]
tomA <- mat2adj(cbind(t(tfExpA), t(gexpA)), method="TOM")[1:numTFs,-(1:numTFs)]
tomB <- mat2adj(cbind(t(tfExpB), t(gexpB)), method="TOM")[1:numTFs,-(1:numTFs)]

aracneA <- mat2adj(cbind(t(tfExpA), t(gexpA)), method="ARACNE")[1:numTFs,-(1:numTFs)]
aracneB <- mat2adj(cbind(t(tfExpB), t(gexpB)), method="ARACNE")[1:numTFs,-(1:numTFs)]
clrA <- mat2adj(cbind(t(tfExpA), t(gexpA)), method="CLR")[1:numTFs,-(1:numTFs)]
clrB <- mat2adj(cbind(t(tfExpB), t(gexpB)), method="CLR")[1:numTFs,-(1:numTFs)]
# mineA <- mat2adj(cbind(t(tfExpA), t(gexpA)), method="MINE", measure="MIC")[1:numTFs,-(1:numTFs)]
# mineB <- mat2adj(cbind(t(tfExpB), t(gexpB)), method="MINE", measure="MIC")[1:numTFs,-(1:numTFs)]


# 9/28/15
# Adding simulated motif data

# Creating motifs with AUC ~ .574
motifs_1 <- matrix(rbinom(length(true_edges), 1, .1+true_edges/5), nrow=numTFs)
rownames(motifs_1) <- TFNames
colnames(motifs_1) <- geneNames
TFsZeros <- diag(numTFs)
rownames(TFsZeros) <- TFNames
colnames(TFsZeros) <- TFNames
motifs_1 <- cbind(motifs_1,TFsZeros)
motifs_1Melt <- melt(motifs_1)
colnames(motifs_1Melt) <- c("V1","V2","value")

# Creating motifs with AUC ~ .54

motifs_2 <- matrix(rbinom(length(true_edges), 1, .1+true_edges/10), nrow=numTFs)
rownames(motifs_2) <- TFNames
colnames(motifs_2) <- geneNames
TFsZeros <- diag(numTFs)
rownames(TFsZeros) <- TFNames
colnames(TFsZeros) <- TFNames
motifs_2 <- cbind(motifs_2,TFsZeros)
motifs_2Melt <- melt(motifs_2)
colnames(motifs_2Melt) <- c("V1","V2","value")

# Run BERE and PANDA on each of the motifs
gexpAWithTFs <- rbind(gexpA,tfExpA)
gexpBWithTFs <- rbind(gexpB,tfExpB)

bere_motif1_AllA <- bere(motifs_1Melt, gexpAWithTFs, score="notincluded")#[TFNames,geneNames]
bere_motif1_AllB <- bere(motifs_1Melt, gexpBWithTFs, score="notincluded")#[TFNames,geneNames]
bere_motif1_AllA <- bere_motif1_AllA[TFNames,c(TFNames,geneNames)]
bere_motif1_AllB <- bere_motif1_AllB[TFNames,c(TFNames,geneNames)]
bereRes_motif1_A <- bere_motif1_AllA[,geneNames]
bereRes_motif1_B <- bere_motif1_AllB[,geneNames]
panda_motif1_AllA <- panda(motifs_1Melt, gexpAWithTFs, hamming = 1e-02)@regNet
panda_motif1_AllB <- panda(motifs_1Melt, gexpBWithTFs, hamming = 1e-02)@regNet
pandaRes_motif1_A <- panda_motif1_AllA[,geneNames]
pandaRes_motif1_B <- panda_motif1_AllB[,geneNames]

bere_motif2_AllA <- bere(motifs_2Melt, gexpAWithTFs, score="notincluded")#[TFNames,geneNames]
bere_motif2_AllB <- bere(motifs_2Melt, gexpBWithTFs, score="notincluded")#[TFNames,geneNames]
bere_motif2_AllA <- bere_motif2_AllA[TFNames,c(TFNames,geneNames)]
bere_motif2_AllB <- bere_motif2_AllB[TFNames,c(TFNames,geneNames)]
bereRes_motif2_A <- bere_motif2_AllA[,geneNames]
bereRes_motif2_B <- bere_motif2_AllB[,geneNames]
panda_motif2_AllA <- panda(motifs_2Melt, gexpAWithTFs, hamming = 1e-02)@regNet
panda_motif2_AllB <- panda(motifs_2Melt, gexpBWithTFs, hamming = 1e-02)@regNet
pandaRes_motif2_A <- panda_motif2_AllA[,geneNames]
pandaRes_motif2_B <- panda_motif2_AllB[,geneNames]

#######################################################
## Calculating the direct results
## Compare to default strategy of pairwise results
#######################################################

corTFA <- mat2adj(t(tfExpA))
corTFB <- mat2adj(t(tfExpB))
directTMcor <- corTFB-corTFA
wgcnaTFA <- mat2adj(t(tfExpA), method="WGCNA")
wgcnaTFB <- mat2adj(t(tfExpB), method="WGCNA")
directTMwgcna <- wgcnaTFB-wgcnaTFA
tomTFA <- mat2adj(t(tfExpA), method="TOM")
tomTFB <- mat2adj(t(tfExpB), method="TOM")
directTMTOM <- tomTFB-tomTFA
wgcna12TFA <- mat2adj(t(tfExpA), method="WGCNA",P=12)
wgcna12TFB <- mat2adj(t(tfExpB), method="WGCNA",P=12)
directTMwgcna12 <- wgcna12TFB-wgcna12TFA
aracneTFA <- mat2adj(t(tfExpA), method="ARACNE")
aracneTFB <- mat2adj(t(tfExpB), method="ARACNE")
directTMaracne <- aracneTFB-aracneTFA
clrTFA <- mat2adj(t(tfExpA), method="CLR")
clrTFB <- mat2adj(t(tfExpB), method="CLR")
directTMclr <- clrTFB-clrTFA
directTMpanda_motif1 <- panda_motif1_AllB[,-1:-numGenes]-panda_motif1_AllA[,-1:-numGenes]
directTMpanda_motif1 <- directTMpanda_motif1[,rownames(directTMpanda_motif1)]
directTMpanda_motif2 <- panda_motif2_AllB[,-1:-numGenes]-panda_motif2_AllA[,-1:-numGenes]
directTMpanda_motif2 <- directTMpanda_motif2[,rownames(directTMpanda_motif2)]
# Issue is here: Bere does not have TFs in results
directbereRes_motif1 <- bere_motif1_AllB[,TFNames]-bere_motif1_AllA[,TFNames]
directbereRes_motif1 <- directbereRes_motif1[TFNames,TFNames]
directbereRes_motif2 <- bere_motif2_AllB[,TFNames]-bere_motif2_AllA[,TFNames]
directbereRes_motif2 <- directbereRes_motif2[TFNames,TFNames]

networkAUCROC <- function(netA, title=""){
    methodPred  <- prediction(c(netA), c(matrixA_GS)>.5)
    roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
    auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
    plot(roc.methodPred, main=title, col = 1, lwd=3)
    
    p.value <- t.test(c(netA), c(matrixB_GS)>.5)$p.value
    if(p.value<2e-16){
        p.value <- "p < 2e-16"
    } else {
        p.value <- paste0("p = ", p.value)
    }
    abline(0,1)
    legend("bottomright", paste0("AUC = ",round(auc.methodPred,4)," \n(",p.value,")"), lty=1,lwd=5,col=1, cex=2)
}
# getAUCROC <- function(netA){
#     methodPred  <- prediction(c(netA), c(matrixA_GS)>.5)
#     roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
#     c(performance(methodPred, "auc")@y.values[[1]])
# }
png("./TM_manuscript/simulated_network_roc.png", width=1800, height=1000)
par(mfrow=c(3,3),oma=c(0,0,4,0),cex.main = 4)
networkAUCROC(pearsonAdjMatA, "Correlation")
networkAUCROC(wgcnaA6, "WGCNA")
networkAUCROC(aracneA, "ARACNE")
networkAUCROC(clrA, "CLR")
networkAUCROC(tomA, "TOM")
# networkAUCROC(motifs_1[,1:numGenes], "Motifs (1)")
# networkAUCROC(motifs_2[,1:numGenes], "Motifs (2)")
networkAUCROC(pandaRes_motif1_A, "PANDA*")
networkAUCROC(pandaRes_motif2_A, "PANDA**")
networkAUCROC(bereRes_motif1_A, "BERE*")
networkAUCROC(bereRes_motif2_A, "BERE**")
title("ROC plots for simulated networks", outer=TRUE)
dev.off()

#### Calculated transition matrices

tm.gs               <- transformation.matrix(matrixA_GS, matrixB_GS,         remove.diagonal=T, standardize=F, method="ols")
tm.pearson          <- transformation.matrix(pearsonAdjMatA, pearsonAdjMatB, remove.diagonal=T, standardize=F, method="ols")
tm.wgcna6           <- transformation.matrix(wgcnaA6, wgcnaB6,               remove.diagonal=T, standardize=F, method="ols")
tm.wgcna12          <- transformation.matrix(wgcnaA12, wgcnaB12,             remove.diagonal=T, standardize=F, method="ols")
tm.aracne           <- transformation.matrix(aracneA, aracneB,               remove.diagonal=T, standardize=F, method="ols")
tm.clr              <- transformation.matrix(clrA, clrB,                     remove.diagonal=T, standardize=F, method="ols")
tm.tom              <- transformation.matrix(tomA, tomB,                     remove.diagonal=T, standardize=F, method="ols")
tm.panda_motif1     <- transformation.matrix(pandaRes_motif1_A, pandaRes_motif1_B,        remove.diagonal=T, standardize=F, method="ols")
tm.panda_motif2     <- transformation.matrix(pandaRes_motif2_A, pandaRes_motif2_B,        remove.diagonal=T, standardize=F, method="ols")
tm.bere_motif1      <- transformation.matrix(bereRes_motif1_A, bereRes_motif1_B,          remove.diagonal=T, standardize=F, method="ols")
tm.bere_motif2      <- transformation.matrix(bereRes_motif2_A, bereRes_motif2_B,          remove.diagonal=T, standardize=F, method="ols")

#### Plot transition matrices heatmaps

heatmap.2(tm.gs, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)
heatmap.2(tm.pearson, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)
heatmap.2(tm.wgcna6, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)
heatmap.2(tm.panda_motif1, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)
heatmap.2(tm.bere_motif1, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)

plotStrengthVsTransitionRank <- function(tm){
    sortedODM <- sort(abs(c(tm)),decreasing=T)
    transition.results <- cbind(TFAs,TFBs, strengths, diag(tm[TFAs,TFBs]))[order(-strengths),]
    cbind(transition.results[,4],sapply(transition.results[,4], function(x){sum(sortedODM>abs(x))})+1)
    #transition.ranks <- do.call(rbind, sapply(transition.results[,4], function(x){which(abs(x) == sort(abs(c(tm.noisy)),decreasing=T))}))[,1]
    transition.ranks <- sapply(transition.results[,4], function(x){sum(sortedODM>abs(x))})+1
    qplot(strengths, transition.ranks, main="Ranks of off-diagonal mass vs transition strength")+ ylim(0, numTFs^2)    
}

plotStrengthVsTransitionRank(tm.pearson)
plotStrengthVsTransitionRank(tm.wgcna6)
plotStrengthVsTransitionRank(tm.wgcna12)
plotStrengthVsTransitionRank(tm.aracne)
plotStrengthVsTransitionRank(tm.clr)
plotStrengthVsTransitionRank(tm.tom)
plotStrengthVsTransitionRank(tm.panda_motif1)
plotStrengthVsTransitionRank(tm.bere_motif1)

dtfi_true <- apply(tm.gs, 1, function(x){sum(abs(x))})
dtfi_obs <- apply(tm.bere_motif1, 1, function(x){sum(abs(x))})
qplot(dtfi_true, dtfi_obs, main="Differential TF Involvement: Observed vs True")+ stat_smooth()
summary(lm(dtfi_true ~ dtfi_obs))


# ROC for transitions
ROCforTransitions <- function(tm, tm2, alpha=.01,method=""){
    
    goldStandard <- c(tm_true)>alpha
    methodPred  <- prediction(c(tm), goldStandard)
    roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
    auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
    methodPred2  <- prediction(c(tm2), goldStandard)
    roc.methodPred2  <- performance(methodPred2, measure = c("tpr","auc"), x.measure = "fpr")
    auc.methodPred2  <- performance(methodPred2, "auc")@y.values[[1]]
    
    
    p.value <- round(t.test(rank(c(tm))[goldStandard], rank(c(tm))[!goldStandard])$p.value,4)
    if(p.value<2e-16){
        p.value <- "p < 2e-16"
    } else {
        p.value <- paste0("p = ", p.value)
    }
    p.value2 <- round(t.test(rank(abs(c(tm2)))[goldStandard], rank(abs(c(tm2)))[!goldStandard])$p.value,4)
    if(p.value2<2e-16){
        p.value2 <- "p < 2e-16"
    } else {
        p.value2 <- paste0("p = ", p.value2)
    }
    plot(roc.methodPred, main=method, cex.main=2, col = "blue", lwd=3)
    lines(roc.methodPred2@x.values[[1]], roc.methodPred2@y.values[[1]], col = "red", lwd=3)
    abline(0,1)
    legend("bottomright", c(paste0("TM = ",round(auc.methodPred,4)," \n(",p.value,")\n"),paste0("EW = ",round(auc.methodPred2,4)," \n(",p.value2,")")), 
           lty=1,lwd=5,col=c("blue","red"), title="Area under the curve",cex=2)
}
# ROCforTransitions(tm.gs,method="Gold Standard")

cutoff <- .4
png('./TM_manuscript/simulation_transition_predictions.png', width=1800, height=1000)
par(mfrow=c(2,5),oma=c(0,0,4,0),cex.main = 4)
ROCforTransitions(tm.pearson, directTMcor, cutoff, method="Correlation Network")
ROCforTransitions(tm.wgcna6, directTMwgcna, cutoff, method="WGCNA")
ROCforTransitions(tm.wgcna12, directTMwgcna12, cutoff, method="WGCNA (12)")
ROCforTransitions(tm.aracne, directTMaracne, cutoff, method="ARACNE")
ROCforTransitions(tm.clr, directTMclr, cutoff, method="CLR")
ROCforTransitions(tm.tom, directTMTOM, cutoff, method="TOM")
ROCforTransitions(tm.panda_motif1, directTMpanda_motif1, cutoff, method="PANDA*")
ROCforTransitions(tm.bere_motif1, directbereRes_motif2, cutoff, method="BERE*")
ROCforTransitions(tm.panda_motif2, t(directTMpanda_motif2), cutoff, method="PANDA**")
ROCforTransitions(tm.bere_motif2, directbereRes_motif2, cutoff, method="BERE**")
title("ROC plots for recovering transitions", outer=TRUE)
dev.off()


# combinedGEXP <- cbind(gexpA, gexpB)
# combinedTFEXP <- cbind(tfExpA,tfExpB) 
# null.networks <- lapply(1:1000, function(x){
#     groupA <- sample(1:(numSamples*2),numSamples)
#     res1 <- cor(t(combinedTFEXP[,groupA]), t(combinedGEXP[,groupA]))
#     res2 <- cor(t(combinedTFEXP[,-groupA]), t(combinedGEXP[,-groupA]))
#     list(res1,res2) 
# })
# # Calculate the transformation matrix for the null data
# tm.null <- lapply(null.networks, function(x){
#     tm <- transformation.matrix(x[[1]],x[[2]],method="ols",remove.diagonal=T, standardize=F)
#     rownames(tm) <- TFNames
#     colnames(tm) <- TFNames
#     tm
# })
# rownames(tm.noisy) <- TFNames
# colnames(tm.noisy) <- TFNames
# 
# ssodm.plot(tm.noisy, tm.null, rescale=T)


# 
# ### penalized TM
# net1 <- t(pearsonAdjMatA)
# net2 <- t(pearsonAdjMatB)
# colFactor <- colSums(net1)/colSums(net2)
# net2.star <- net2-sweep(net1, 2, colFactor, '*')
# tm.pen <- sapply(1:numTFs, function(i){
#     z <- penalized(net1[,i], net2.star, lambda1=1, model="linear", standardize=T)
#     coefficients(z, "penalized")
# })
