library(MASS)
library(poweRlaw)
library(actuar)
library(reshape2)
library(penalized)

transformation.matrix <- function(network.1, network.2, by.tfs=T, standardize=T, remove.diagonal=T, method="ols"){
    if(is.list(network.1)&&is.list(network.2)){
        if(by.tfs){
            net1 <- t(network.1$reg.net)
            net2 <- t(network.2$reg.net)
        } else {
            net1 <- network.1$reg.net
            net2 <- network.2$reg.net
        }
    } else if(is.matrix(network.1)&&is.matrix(network.2)){
        if(by.tfs){
            net1 <- t(network.1)
            net2 <- t(network.2)
        } else {
            net1 <- network.1
            net2 <- network.2
        }
    } else {
        stop("Networks must be lists or matrices")
    }
    #gene.trans.matrix <- svd(net2)$v %*% diag(1/svd(net2)$d) %*% t(svd(net2)$u) %*% net1
    if (method == "kabsch"){
        tf.trans.matrix <- kabsch(net1,net2)
    }
    if (method == "old"){
        svd.net2 <- svd(net2)
        tf.trans.matrix <- svd.net2$v %*% diag(1/svd.net2$d) %*% t(svd.net2$u) %*% net1
    }
    if (method == "ols"){
        tf.trans.matrix <- ginv(t(net1)%*%net1)%*%t(net1)%*%net2
        print("Using OLS method")
    }
    if (standardize){
        tf.trans.matrix <- apply(tf.trans.matrix, 1, function(x){
            #   x.zero <- (x-mean(x))
            x/sum(abs(x))
        })
    }
    
    if (remove.diagonal){
        diag(tf.trans.matrix) <- 0
    }
    # Add column labels
    colnames(tf.trans.matrix) <- rownames(tf.trans.matrix)
    tf.trans.matrix
}
bereFull <- function(motifs, exprData, alpha=.5, penalized=T, lambda=10, score="motifincluded"){
    require(MASS)
    exprData <- data.frame(exprData)
    tfdcast <- dcast(motifs,V1~V2,fill=0)
    rownames(tfdcast) <- tfdcast[,1]
    tfdcast <- tfdcast[,-1]
    
    exprData <- exprData[sort(rownames(exprData)),]
    tfdcast <- tfdcast[,sort(colnames(tfdcast)),]
    tfNames <- rownames(tfdcast)[rownames(tfdcast) %in% rownames(exprData)]
    
    ## Filtering
    # filter out the TFs that are not in expression set
    tfdcast <- tfdcast[rownames(tfdcast)%in%tfNames,]
    
    # check that IDs match
    if (prod(rownames(exprData)==colnames(tfdcast))!=1){stop("ID mismatch")}
    
    ## Get direct evidence
    directCor <- t(cor(t(exprData),t(exprData[rownames(exprData)%in%tfNames,]))^2)
    
    ## Get the indirect evidence    
    result <- t(apply(tfdcast, 1, function(x){
        cat(".")
        tfTargets <- as.numeric(x)
        
        # Ordinary Logistic Reg
        #         z <- glm(tfTargets ~ ., data=exprData, family="binomial")
        
        # Penalized Logistic Reg
        z <- penalized(tfTargets, exprData, lambda2=lambda, model="logistic", standardize=T)
        #         z <- optL1(tfTargets, exprData, minlambda1=25, fold=5)
        
        
        predict(z, exprData)
    }))
    
    ## Convert values to ranks
    directCor <- matrix(rank(directCor), ncol=ncol(directCor))
    result <- matrix(rank(result), ncol=ncol(result))
    
    consensus <- directCor*(1-alpha) + result*alpha
    rownames(consensus) <- rownames(tfdcast)
    colnames(consensus) <- rownames(exprData)
    #     if(score=="motifincluded"){
    #         consensus <- as.matrix(result + tfdcast)
    #     }
    consensus
}



numTFs <- 50
numGenes <- 1000
numSamples <- 500
numTransitions <- 5

degrees <- floor(rpareto(numTFs,20,10000))
degrees[degrees<10] <- 10
degrees[degrees>(numGenes/4)] <- numGenes/4
adjMat <- sapply(degrees,function(x){
  sample(c(rep(1,x),rep(0,numGenes-x)))
})
selectedTFs <- sample(1:numTFs,numTransitions)
rownames(adjMat) <- paste("GENE",1:numGenes,sep="")
colnames(adjMat) <- paste("TF",1:numTFs,sep="")
rownames(adjMat)[1:numTFs]<- colnames(adjMat)[1:numTFs]

# Included edges
adjMat[selectedTFs,selectedTFs]<- 1
Sigma <- adjMat%*%t(adjMat)
diag(Sigma) <- diag(Sigma)+2
exprData1 <- t(mvrnorm(numSamples, mu=rep(0,numGenes), Sigma=Sigma))
rownames(exprData1)<- rownames(adjMat)

# Excluded edges
adjMat[selectedTFs,selectedTFs] <- 0
Sigma <- adjMat%*%t(adjMat)
diag(Sigma) <- diag(Sigma)+10
exprData2 <- t(mvrnorm(numSamples, mu=rep(0,numGenes), Sigma=Sigma))
rownames(exprData2)<- rownames(adjMat)

# set motif prior adjMat
adjMat[selectedTFs,selectedTFs]<-rbinom(numTransitions^2,1,.5)
adjMat <- data.frame(adjMat)
adjMat$id <- rownames(adjMat)
adjMatMelt <- melt(adjMat, "id")[,c(2,1,3)]
colnames(adjMatMelt)<-c("V1","V2","V3") 

bereScores1 <- bereFull(adjMatMelt, exprData1)
bereScores2 <- bereFull(adjMatMelt, exprData2)

tm <- transformation.matrix(bereScores1, bereScores2, remove.diagonal=T, method="ols")

exprData <- cbind(exprData1,exprData2)

null.tms <- replicate(1000,{
  partition <- sample(1:(2*numSamples))
  nullbereScores1 <- bereFull(adjMatMelt, exprData[,partition[1:numSamples]])
  nullbereScores2 <- bereFull(adjMatMelt, exprData[,partition[(numSamples+1):(numSamples*2)]])
  transformation.matrix(nullbereScores1, nullbereScores2, remove.diagonal=T, method="ols")
})
null.tms.list <- lapply(seq(dim(null.tms)[3]), function(x) null.tms[ , , x])

ssodm.plot(tm, null.tms.list, plot.title="Simulated Data: dTFI")
ssodm.plot(tm, null.tms.list, rescale=T, plot.title="SSODM observed and null, COPD subjects vs Smoker control (R)")
ssodm.plot(tm, null.tms.list, rescale=T, plot.title="SSODM observed and null, COPD subjects vs Smoker control (R)", 
           highlight.tfs = rownames(exprData2)[selectedTFs])

save.image(file="sim_network_TM")