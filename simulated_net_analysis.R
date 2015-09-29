library(MASS)
library(pROC)
library(reshape2)
library(PANDA)
library(bereR)


numGenes   <- 200
numTFs     <- 40
numROA     <- 800
numSamples <- 50
propOfMotifs <- .8

acceptable.start <- F
for(i in 1:20){
  true_binary <- sample(c(rep(1,numROA),rep(0,numGenes*numTFs-numROA)),numGenes*numTFs)
  trueROA  <- matrix(true_binary, nrow=numGenes,ncol=numTFs)
  unmatched.genes <- sum(apply(trueROA,1,function(x){
    prod(x==0)
  }))
  unmatched.TFs <- sum(apply(trueROA,2,function(x){
    prod(x==0)
  }))
  if (unmatched.genes==0 && unmatched.TFs==0){
    acceptable.start <- T
    break
  }
}
if(!acceptable.start){
  # This needed because random matrices so not always meet sufficient starting conditions.
  print("Error: Unmatched genes or TFs exist")
}
trueCoop <- t(trueROA)%*%trueROA

trueCoaf <- trueROA%*%t(trueROA)

#  This line for adding additional noise
# trueCoaf <- trueCoaf + 10*diag(numGenes)

gene_exp_sim <- data.frame(t(mvrnorm(n=numSamples,mu=rep(0,numGenes),trueCoaf)))
row.names(gene_exp_sim) <- paste("Gene",1:numGenes)
row.names(trueROA) <- paste("Gene",1:numGenes)
colnames(trueROA)  <- paste("TF",1:numTFs)
names(gene_exp_sim) <- paste("Sample",1:numSamples)

# Generate randomized error
type2Error = .5
type1Error = .1
# Remove a percentage of them
observed_binary <- true_binary
true_ind1 <- observed_binary==1
true_ind0 <- observed_binary==0
observed_binary[true_ind1] <- sample(c(rep(1,sum(true_ind1)*(1-type2Error)),rep(0,sum(true_ind1)*type2Error)),sum(true_ind1))
observed_binary[true_ind0] <- sample(c(rep(0,sum(true_ind0)*(1-type1Error)),rep(1,sum(true_ind0)*type1Error)),sum(true_ind0))
starting.RegNet  <- matrix(observed_binary,nrow=numGenes,ncol=numTFs)
rownames(starting.RegNet) <- paste("Gene",1:numGenes)
colnames(starting.RegNet) <- paste("TF",1:numTFs)
observedCoop <- t(starting.RegNet)%*%starting.RegNet

motif.data <- melt(t(starting.RegNet))


# Generate binary PPI table from cooperative net
ppiTable <- trueCoop
diag(ppiTable) <- 0
#sort(as.vector(ppiTable))
ppicutoff <- as.numeric(quantile(ppiTable, probs=propOfMotifs))
ppiTable <- matrix(as.numeric(ppiTable>ppicutoff), nrow=numTFs,ncol=numTFs)

# Generate PPI subset
propMissing = .75
observedPPIPrior  <- ppiTable*matrix(sample(c(rep(1,length(ppiTable)*(1-propMissing)),rep(0,length(ppiTable)*(propMissing))),length(ppiTable)),nrow=numTFs,ncol=numTFs)
observedPPIPrior  <- observedPPIPrior+t(observedPPIPrior)
observedPPIPrior[observedPPIPrior>1] <- 1
rownames(observedPPIPrior) <- paste("TF",1:numTFs)
colnames(observedPPIPrior) <- paste("TF",1:numTFs)

ppi.data <- melt(observedPPIPrior)

roc(trueROA,starting.RegNet)

# png(file="sim_all_edges.png")
RegNet <- t(panda(motif.data,gene_exp_sim)@reg.net)[rownames(trueROA),colnames(trueROA)]
RegNet2 <- t(regpredict(motif.data,gene_exp_sim))[rownames(trueROA),colnames(trueROA)]

bere.auc <- plot.roc(trueROA,RegNet2, col="blue",main="Bipartite simulation")$auc
# bere.auc <- plot.roc(trueROA,RegNet2+5*starting.RegNet, col="blue",main="Bipartite simulation")$auc
panda.auc <- lines.roc(trueROA,RegNet, col="red")$auc
motif.auc <- lines.roc(trueROA,starting.RegNet, col="grey")$auc
legend(.6,.2,
       c(paste("'BERE' method",round(bere.auc,4)),paste("PANDA",round(panda.auc,4)),paste("motif",round(motif.auc,4))),
       lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("blue","red","grey"))
# dev.off()


