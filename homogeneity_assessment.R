library(penalized)
library(reshape2)
library(tidyr)
# library(pandaR)
# library(bptools)
copd.filename <- "null.networks_all.rds"
eclipse.filename <- "eclipse.networks.rds"
local.wd <- "~/gd/Harvard/Research/"
data.dir  <- "./data/Eclipse/"
setwd(local.wd)
setwd(data.dir)

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



######################################################
##      Data Loading from ECLIPSE dataset          ###
##                                                 ###
######################################################
eclipse <- list()
eclipse$motif    <- read.table("ECLIPSE_Blood_Motif.txt",header=F)
eclipse$exp      <- read.table("ECLIPSE_Blood_Exp.txt",row.names=1,header=T)
eclipse$ppi      <- read.table("OV_PPI.txt",header=F)
eclipse$clinical <- read.table("ECLIPSE_blood.txt",header=T,fill = TRUE, sep="\t",row.names=1)
eclipse$exp      <- eclipse$exp[,order(colnames(eclipse$exp))]  # Make sure expression and clinical is in same order
eclipse$clinical <- eclipse$clinical[colnames(eclipse$exp),]    # Make sure clinical only contains patients with expression data

# Specify the group partition
filter.vec.1 <- eclipse$clinical$Subject.type=="COPD Subjects"
filter.vec.2 <- eclipse$clinical$Subject.type=="Smoker Controls"

permutations <- 10
resMatrix <- matrix(NA, nrow=permutations*4, ncol=164)
nullMatrix <- matrix(NA, nrow=permutations*2, ncol=164)
hetMatrix <- matrix(NA, nrow=permutations*2, ncol=164)
for (i in 0:(permutations-1)){
  print("********************** iteration ***********************")
  print(i)
  sampleSize <- 40
  homosubset1 <- sample(which(filter.vec.1),sampleSize*2)
  net1 <- bereFull(eclipse$motif,eclipse$exp[,homosubset1[1:sampleSize]])
  net2 <- bereFull(eclipse$motif,eclipse$exp[,homosubset1[(sampleSize+1):(sampleSize*2)]])
  homosubset2 <- sample(which(filter.vec.2),sampleSize*2)
  net3 <- bereFull(eclipse$motif,eclipse$exp[,homosubset2[1:sampleSize]])
  net4 <- bereFull(eclipse$motif,eclipse$exp[,homosubset2[(sampleSize+1):(sampleSize*2)]])
  hetsubset <- sample(c(homosubset1,homosubset2))
  net5 <- bereFull(eclipse$motif,eclipse$exp[,hetsubset[1:sampleSize]])
  net6 <- bereFull(eclipse$motif,eclipse$exp[,hetsubset[(sampleSize+1):(sampleSize*2)]])
  net7 <- bereFull(eclipse$motif,eclipse$exp[,hetsubset[(2*sampleSize+1):(sampleSize*3)]])
  net8 <- bereFull(eclipse$motif,eclipse$exp[,hetsubset[(3*sampleSize+1):(sampleSize*4)]])
  
  
  tm <- transformation.matrix(net1, net3, remove.diagonal=T, method="ols")
  resMatrix[i*4+1,] <- apply(tm,2,function(x){t(x)%*%x})
  tm <- transformation.matrix(net1, net4, remove.diagonal=T, method="ols")
  resMatrix[i*4+2,] <- apply(tm,2,function(x){t(x)%*%x})
  tm <- transformation.matrix(net2, net3, remove.diagonal=T, method="ols")
  resMatrix[i*4+3,] <- apply(tm,2,function(x){t(x)%*%x})
  tm <- transformation.matrix(net2, net4, remove.diagonal=T, method="ols")
  resMatrix[i*4+4,] <- apply(tm,2,function(x){t(x)%*%x})
  
  tm <- transformation.matrix(net1, net2, remove.diagonal=T, method="ols")
  nullMatrix[i*2+1,] <- apply(tm,2,function(x){t(x)%*%x})
  tm <- transformation.matrix(net3, net4, remove.diagonal=T, method="ols")
  nullMatrix[i*2+2,] <- apply(tm,2,function(x){t(x)%*%x})
  
  tm <- transformation.matrix(net5, net6, remove.diagonal=T, method="ols")
  hetMatrix[i*2+1,] <- apply(tm,2,function(x){t(x)%*%x})
  tm <- transformation.matrix(net7, net8, remove.diagonal=T, method="ols")
  hetMatrix[i*2+2,] <- apply(tm,2,function(x){t(x)%*%x})
  
}

saveRDS(list(resMatrix=resMatrix,nullMatrix=nullMatrix,hetMatrix=hetMatrix), "./homogeneity_results.rds")

# transcriptionFactor <- "E2F3"
# values <- c(resMatrix[,transcriptionFactor], nullMatrix[,transcriptionFactor])
# xvalues <- c(rep("obs",length(resMatrix[,transcriptionFactor])), rep("null",length(nullMatrix[,transcriptionFactor])))
# qplot(y=values,x=xvalues, geom = "boxplot", main=paste("Observed vs Null for",transcriptionFactor))

resValues <- sort(c(resMatrix[,3])[c(F,T)])
homValues <- sort(c(nullMatrix[,3]))


pvalsObsVsHet <- sapply(1:ncol(resMatrix), function(i){
  tt <- t.test(resMatrix[,i],hetMatrix[,i])
  tt$p.value
})
pvalsHomVsHet <- sapply(1:ncol(nullMatrix), function(i){
  tt <- t.test(nullMatrix[,i],hetMatrix[,i])
  tt$p.value
})

pvalsHomSMC <- sapply(seq(2,ncol(nullMatrix),2), function(i){
  tt <- t.test(nullMatrix[,i],hetMatrix[,i])
  tt$p.value
})
pvalsHomCOPD <- sapply(seq(1,ncol(nullMatrix),2), function(i){
  tt <- t.test(nullMatrix[,i],hetMatrix[,i])
  tt$p.value
})
library(gridExtra)
plot1 <- qplot(pvalsObsVsHet, binwidth=.04, main="p-value plot for t-test of dTFIs of observed vs null")
plot2 <- qplot(pvalsHomVsHet, binwidth=.04, main="p-value plot for t-test of dTFIs of homogenous vs null")
grid.arrange(plot1, plot2, ncol=2)
plot1 <- qplot(pvalsHomCOPD, binwidth=.04, main="COPD null")
plot2 <- qplot(pvalsHomSMC, binwidth=.04, main="Smoker Control null")
grid.arrange(plot1, plot2, ncol=2)
qplot(sort(pvalsHomCOPD), sort(pvalsHomSMC), main="QQ plot: COPD vs SMC (all null)")+ geom_abline(intercept = 0, slope = 1)

print(qplot(y=resValues, x=homValues, main=paste("QQ plot, homogenous networks vs hetergenous networks (all null)"))+ geom_abline(intercept = 0, slope = 1))

homValues <- sort(c(nullMatrix[,3])[c(F,T)])
hetValues <- sort(c(hetMatrix[,3]))

print(qplot(y=hetValues, x=homValues, main=paste("QQ plot, homogenous networks vs hetergenous networks (all null)"))+ geom_abline(intercept = 0, slope = 1))

qplot(colMeans(hetMatrix) - colMeans(nullMatrix))
qplot(colMeans(resMatrix) - colMeans(hetMatrix))


png(filename="./homogeneityQQ.png")
print(qplot(y=hetValues, x=homValues, main=paste("QQ plot, homogenous networks vs hetergenous networks (all null)"))+ geom_abline(intercept = 0, slope = 1))
dev.off()
