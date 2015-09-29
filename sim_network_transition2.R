library(bptools)
library(gplots)
library(ROCR)

numTFs <- 100
numGenes <- 10000
numTransitions <- 100
sds <- 2
runSim <- function(numTFs, numGenes, numTransitions, sds){
    matrixA_GS <- matrix(rnorm(numGenes*numTFs),nrow=numTFs)
    matrixB_GS <- matrixA_GS# + matrix(rnorm(10000)/10,nrow=10)
    TFAs <- sample(1:numTFs,numTransitions, replace=T)
    TFBs <- sample(1:numTFs,numTransitions, replace=T)
    TFBs[TFAs==TFBs] <- sample(1:numTFs,sum(TFAs==TFBs), replace=T)
    sum(sum(TFAs==TFBs))
    strengths <- (c(-50:-1,1:50))/50
    
    for(i in 1:numTransitions){
        matrixB_GS[TFBs[i],] <- matrixB_GS[TFBs[i],] + matrixA_GS[TFAs[i],]*strengths[i]
    }
    
    # Introduce noise
    matrixA <- matrixA_GS + matrix(rnorm(numGenes*numTFs,sd=sds),nrow=numTFs)
    matrixB <- matrixB_GS + matrix(rnorm(numGenes*numTFs,sd=sds),nrow=numTFs)
    
    methodPred  <- prediction(c(matrixA,matrixB), c(matrixA_GS,matrixB_GS)>1)
    roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
    auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
    plot(roc.methodPred, main="ROC for true-edges after added noise", col = 1, lwd=3)
    legend(.5,.6, paste("AUCROC =",round(auc.methodPred,4)), lty=1,lwd=5,col=1)
    
    tm.gs <- transformation.matrix(matrixA_GS, matrixB_GS, remove.diagonal=T, standardize=F, method="ols")
    tm.noisy <- transformation.matrix(matrixA, matrixB, remove.diagonal=T, standardize=F, method="ols")

################## Penalized matrix
    net1 <- t(matrixA)
    net2 <- t(matrixB)
    net2.star <- sapply(1:numTFs, function(i,x,y){
        lm(y[,i]~x[,i])$resid
    }, net1, net2)
    tm.pen <- sapply(1:numTFs, function(i){
        z <- optL1(net1[,i], net2.star, fold=5, minlambda1=1, maxlambda1=1000, model="linear", standardize=T)
        coefficients(z$fullfit, "penalized")
    })
    heatmap.2(t(tm.pen), col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)

#######################################

    transition.results <- cbind(TFAs,TFBs, strengths, diag(tm.noisy[TFAs,TFBs]))[order(-strengths),]
    heatmap.2(tm.gs, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)
    heatmap.2(tm.noisy, col=colorRampPalette(c("blue", "white", "red"))(n = 1000), density.info="none", trace="none", dendrogram="none", Colv=FALSE, Rowv=FALSE)
    sortedODM <- sort(abs(c(tm.noisy)),decreasing=T)
    cbind(transition.results[,4],sapply(transition.results[,4], function(x){sum(sortedODM>abs(x))})+1)
    #transition.ranks <- do.call(rbind, sapply(transition.results[,4], function(x){which(abs(x) == sort(abs(c(tm.noisy)),decreasing=T))}))[,1]
    transition.ranks <- sapply(transition.results[,4], function(x){sum(sortedODM>abs(x))})+1
    qplot(strengths, transition.ranks, main="Ranks of off-diagonal mass vs transition strength")+ ylim(0, numTFs^2)
    #p + geom_smooth(method = "loess", size = 1)
    
    dtfi_true <- rep(0,numTFs)
    for(i in 1:numTransitions){
        dtfi_true[TFAs[i]] <- dtfi_true[TFAs[i]] + (strengths[i])^2
    }
    dtfi_obs <- apply(tm.noisy, 1, function(x){t(x)%*%x})
    qplot(dtfi_true, dtfi_obs, main="Differential TF Involvement: Observed vs True")+ stat_smooth()
    summary(lm(dtfi_true ~ dtfi_obs))
    
    # ROC for transitions
    TM_GS <- matrix(0,nrow=numTFs,ncol=numTFs)
    for(i in 1:numTransitions){
        TM_GS[TFAs[i],TFBs[i]]<-1
    }
    
    methodPred  <- prediction(abs(c(tm.noisy)), c(TM_GS)>=1)
    roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
    auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
    plot(roc.methodPred, main="ROC for transitions", col = 1, lwd=3)
    legend(.5,.2, paste("AUCROC =",round(auc.methodPred,4)), lty=1,lwd=5,col=1)
    list("roc.methodPred"=roc.methodPred,"numGenes"=numGenes, "signal"=1/((sds^2)+1))
}


res <- sapply(c(500,1000,2000,10000), function(x){
    sapply(c(1,2,3,9), function(y,x){
        list(runSim(numTFs, x, numTransitions, y))
    }, x)
})

lineTypes <- 1:4
names(lineTypes) <- c(500,1000,2000,10000)
colorPal <- c("black","red","blue","yellow")
names(colorPal) <- 1/(c(1,2,3,9)^2+1)
length(res)
plot(res[[1]][["roc.methodPred"]], main="ROC for transitions", col = 1, lwd=3)
lapply(res, function(x){
    lines(x[["roc.methodPred"]]@x.values[[1]], x[["roc.methodPred"]]@y.values[[1]], lty=lineTypes[as.character(x[["numGenes"]])], col = colorPal[as.character(x[["signal"]])], lwd=3)
})
