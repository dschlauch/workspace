# For creating gene expression data
library(WGCNA)
library(bptools)
library(gplots)
library(ROCR)
library(penalized)
library(tidyr)
library(bereR)

numTFs <- 100
numGenes <- 5000
numTransitions <- 100
numSamples <- 500

true_edges <- runif(numGenes*numTFs)
matrixA_GS <- matrix(true_edges,nrow=numTFs)
matrixB_GS <- matrixA_GS# + matrix(rnorm(10000)/10,nrow=10)
probs <- rexp(numTFs)
TFAs <- sample(1:numTFs,numTransitions, replace=T)
TFBs <- sample(1:numTFs,numTransitions, replace=T, probs)
TFBs[TFAs==TFBs] <- sample(1:numTFs,sum(TFAs==TFBs), replace=T)
sum(sum(TFAs==TFBs))
strengths <- (c(-50:-1,1:50))/50

for(i in 1:numTransitions){
    matrixB_GS[TFBs[i],] <- matrixB_GS[TFBs[i],] + matrixA_GS[TFAs[i],]*strengths[i]
}



varcovA <- t(matrixA_GS)%*%matrixA_GS
varcovB <- t(matrixB_GS)%*%matrixB_GS
gexpA <- t(mvrnorm(n=numSamples, mu=rep(0,numGenes), Sigma=varcovA))
gexpB <- t(mvrnorm(n=numSamples, mu=rep(0,numGenes), Sigma=varcovB))
rownames(gexpA) <- paste0("Gene",1:numGenes)
rownames(gexpB) <- paste0("Gene",1:numGenes)
tfExpA <- matrixA_GS %*% gexpA
tfExpB <- matrixB_GS %*% gexpB
rownames(tfExpA) <- paste0("TF",1:numTFs)
rownames(tfExpB) <- paste0("TF",1:numTFs)
pearsonAdjMatA <- cor(t(tfExpA), t(gexpA))
pearsonAdjMatB <- cor(t(tfExpB), t(gexpB))


# 9/28/15
# Adding simulated motif data
motifs <- matrix(rbinom(length(true_edges), 1, true_edges/4), nrow=numTFs)
rownames(motifs) <- paste0("TF",1:numTFs)
colnames(motifs) <- paste0("Gene",1:numGenes)
TFsZeros <- matrix(0,nrow=numTFs,ncol=numTFs)
rownames(TFsZeros) <- paste0("TF",1:numTFs)
colnames(TFsZeros) <- paste0("TF",1:numTFs)
motifs <- cbind(motifs,TFsZeros)
motifsMelt <- melt(motifs)
colnames(motifsMelt) <- c("V1","V2","value")
gexpAWithTFs <- rbind(gexpA,tfExpA)
gexpBWithTFs <- rbind(gexpB,tfExpB)
bereResA <- bereFull(motifsMelt, gexpAWithTFs, alpha=1.0, score="notincluded")[,paste0("Gene",1:numGenes)]
bereResB <- bereFull(motifsMelt, gexpBWithTFs, alpha=1.0, score="notincluded")[,paste0("Gene",1:numGenes)]
pandaResA <- panda(motifsMelt, gexpAWithTFs)@regNet[,paste0("Gene",1:numGenes)]
pandaResB <- panda(motifsMelt, gexpBWithTFs)@regNet[,paste0("Gene",1:numGenes)]
tm.bere <- transformation.matrix(bereResA, bereResB, remove.diagonal=T, standardize=F, method="ols")
heatmap.2(tm.bere, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)
tm.panda <- transformation.matrix(pandaResA, pandaResB, remove.diagonal=T, standardize=F, method="ols")
heatmap.2(tm.panda, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)

networkAUCROC <- function(netA, netB){
    methodPred  <- prediction(c(netA,netB), c(matrixA_GS,matrixB_GS)>.5)
    roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
    auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
    plot(roc.methodPred, main="ROC for networks", col = 1, lwd=3)
    legend(.5,.6, paste("AUCROC =",round(auc.methodPred,4)), lty=1,lwd=5,col=1)
}
networkAUCROC(pearsonAdjMatA, pearsonAdjMatB)
networkAUCROC(bereResA, bereResB)
networkAUCROC(pandaResA, pandaResA)
networkAUCROC(pearsonAdjMatA+motifs, pearsonAdjMatA+motifs)

tm.gs <- transformation.matrix(matrixA_GS, matrixB_GS, remove.diagonal=T, standardize=F, method="ols")
tm.noisy <- transformation.matrix(pearsonAdjMatA, pearsonAdjMatB, remove.diagonal=T, standardize=F, method="ols")
heatmap.2(tm.gs, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)
heatmap.2(tm.noisy, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)

################## Penalized matrix
penalizedTM <- function(net1, net2){
    net2.star <- sapply(1:numTFs, function(i,x,y){
        lm(y[,i]~x[,i])$resid
    }, net1, net2)
    tm.penL1 <- sapply(1:numTFs, function(i){
        #     z <- optL1(net2.star[,i], net1, fold=5, minlambda1=.1, maxlambda1=2, model="linear", standardize=T)
        #     coefficients(z$fullfit, "penalized")
        z <- penalized(net2.star[,i], net1, lambda1=1, model="linear", standardize=T)
        coefficients(z, "penalized")
    })
    diag(tm.penL1)<-0
    tm.penL1
}
tm.penL1 <- penalizedTM(t(pearsonAdjMatA),t(pearsonAdjMatB))

heatmap.2(tm.penL1, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)

#######################################

plotStrengthVsTransitionRank <- function(tm){
    sortedODM <- sort(abs(c(tm)),decreasing=T)
    transition.results <- cbind(TFAs,TFBs, strengths, diag(tm[TFAs,TFBs]))[order(-strengths),]
    cbind(transition.results[,4],sapply(transition.results[,4], function(x){sum(sortedODM>abs(x))})+1)
    #transition.ranks <- do.call(rbind, sapply(transition.results[,4], function(x){which(abs(x) == sort(abs(c(tm.noisy)),decreasing=T))}))[,1]
    transition.ranks <- sapply(transition.results[,4], function(x){sum(sortedODM>abs(x))})+1
    qplot(strengths, transition.ranks, main="Ranks of off-diagonal mass vs transition strength")+ ylim(0, numTFs^2)    
}

plotStrengthVsTransitionRank(tm.penL1)
plotStrengthVsTransitionRank(tm.noisy)
plotStrengthVsTransitionRank(tm.bere)
plotStrengthVsTransitionRank(tm.panda)

dtfi_true <- apply(tm.gs, 1, function(x){sum(abs(x))})
dtfi_obs <- apply(tm.panda, 1, function(x){sum(abs(x))})
qplot(dtfi_true, dtfi_obs, main="Differential TF Involvement: Observed vs True")+ stat_smooth()
summary(lm(dtfi_true ~ dtfi_obs))


# ROC for transitions
ROCforTransitions <- function(tm){
    TM_GS <- matrix(0,nrow=numTFs,ncol=numTFs)
    for(i in 1:numTransitions){
        TM_GS[TFAs[i],TFBs[i]]<-1
    }
    
    methodPred  <- prediction(abs(c(tm)), c(TM_GS)>=1)
    roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
    auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
    plot(roc.methodPred, main="ROC for transitions", col = 1, lwd=3)
    legend(.5,.2, paste("AUCROC =",round(auc.methodPred,4)), lty=1,lwd=5,col=1)
    abline(0,1)
}
ROCforTransitions(tm.gs)
ROCforTransitions(tm.noisy)
ROCforTransitions(tm.panda)
ROCforTransitions(tm.bere)

combinedGEXP <- cbind(gexpA, gexpB)
combinedTFEXP <- cbind(tfExpA,tfExpB) 
null.networks <- lapply(1:1000, function(x){
    groupA <- sample(1:(numSamples*2),numSamples)
    res1 <- cor(t(combinedTFEXP[,groupA]), t(combinedGEXP[,groupA]))
    res2 <- cor(t(combinedTFEXP[,-groupA]), t(combinedGEXP[,-groupA]))
    list(res1,res2) 
})
# Calculate the transformation matrix for the null data
tm.null <- lapply(null.networks, function(x){
    tm <- transformation.matrix(x[[1]],x[[2]],method="ols",remove.diagonal=T, standardize=F)
    rownames(tm) <- paste0("TF",1:numTFs)
    colnames(tm) <- paste0("TF",1:numTFs)
    tm
})
rownames(tm.noisy) <- paste0("TF",1:numTFs)
colnames(tm.noisy) <- paste0("TF",1:numTFs)

ssodm.plot(tm.noisy, tm.null, rescale=T)



### penalized TM
net1 <- t(pearsonAdjMatA)
net2 <- t(pearsonAdjMatB)
colFactor <- colSums(net1)/colSums(net2)
net2.star <- net2-sweep(net1, 2, colFactor, '*')
tm.pen <- sapply(1:numTFs, function(i){
    z <- penalized(net1[,i], net2.star, lambda1=1, model="linear", standardize=T)
    coefficients(z, "penalized")
})
