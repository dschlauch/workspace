
setwd("C:/Users/Dan/Dropbox/Dan's/Documents/Harvard/Research/data/Data")

# Set Program Parameters
alpha=0.1;
updateMethod <- 'tanimoto'
data.source <- 'file'

if(data.source=='file'){
  #  Read data from file
  exp_file   = "../Data/YeastData/YeastData_SRExpression.txt";
  motif_file = "../Data/YeastData/YeastData_Motif.txt"
  ppi_file   = "../Data/YeastData/YeastData_PPI.txt"
  yeastChip_file   = "../Data/YeastData/YeastData_ChIP.txt"
  Exp.sr        = read.table(exp_file, sep="\t", row.names=1);
  # Randomize gene labels
  # rownames(Exp) <- sample(rownames(Exp))
  Motif           = read.table(motif_file, sep="\t", row.names=NULL);
  PPI             = read.table(ppi_file, sep="\t", row.names=NULL);
  gold.standard  <- read.table(yeastChip_file, sep="\t", row.names=NULL);
} 

if(data.source=='sim'){
  #  Read data from simulation
  Exp = gene_exp_sim 
  Motif = motif.data
  PPI <- ppi.data
}
NormalizeNetwork<-function(X)
{
  if (updateMethod=='danmethod'){
    X[X<0]<-0
  	return(X);
  }
  if (updateMethod=='tanimoto'){
    
    # overall values
    mu0=mean(X);
    std0=sd(X);
    
    # operations on rows
    mu1=apply(X,1,mean); # operations on rows
    std1=apply(X,1,sd)*sqrt((dim(X)[2]-1)/dim(X)[2]);
    mu1=matrix(rep(mu1, dim(X)[2]), dim(X));
    std1=matrix(rep(std1, dim(X)[2]), dim(X));
    Z1=(X-mu1)/std1;
    
    # operations on columns
    mu2=apply(X,2,mean); # operations on columns
    std2=apply(X,2,sd)*sqrt((dim(X)[1]-1)/dim(X)[1]);
    mu2=matrix(rep(mu2, each=dim(X)[1]), dim(X));
    std2=matrix(rep(std2, each=dim(X)[1]), dim(X));
    Z2=(X-mu2)/std2;
    
    # combine and return
    normMat=Z1/sqrt(2)+Z2/sqrt(2);
    
    # Dan fix to NaN
    normMat[is.na(normMat)]<-0
    return(normMat);
    
  }
}

Tfunction<-function(X,Y)
{
  Amat=(X %*% Y);
	Bmat=apply(Y*Y,2,sum);
	Bmat=matrix(rep(Bmat, each=dim(X)[1]), dim(Amat));
	Cmat=apply(X*X,1,sum);
	Cmat=matrix(rep(Cmat, dim(Y)[2]), dim(Amat));

	Amat=Amat/sqrt(Bmat+Cmat-abs(Amat));

	return(Amat);
}

Dfunction<-function(X,Y)
{
  A <- cor(X,Y)
  A[A<0]<-0
  return(A);
}


UpdateDiagonal<-function(diagMat, num, alpha, step)
{
	diagMat[seq(1, num*num, num+1)]=NaN;
	diagstd=apply(diagMat,2,sd,na.rm=T)*sqrt((num-2)/(num-1));
	diagMat[seq(1, num*num, num+1)]=diagstd*num*exp(2*alpha*step);
	return(diagMat);
}


## Read in Data ##
print('Reading in data!')

# Expression Data


NumConditions=dim(Exp)[2];
GeneNames=row.names(Exp);
NumGenes=length(GeneNames);
if(NumConditions<3) {
	print('Not enough expression conditions detected to calculate correlation. Co-regulation network will be initialized to an identity matrix.');
	GeneCoReg=diag(NumGenes);
} else {
	GeneCoReg=cor(t(Exp), method="pearson", use="pairwise.complete.obs");
}

# Prior Regulatory Network

TFNames=unique(Motif[,1]);
NumTFs=length(TFNames);
Idx1=match(Motif[,1], TFNames);
Idx2=match(Motif[,2], GeneNames);
Idx=(Idx2-1)*NumTFs+Idx1;
RegNet=matrix(data=0, NumTFs, NumGenes);
RegNet[Idx]=as.numeric(Motif[,3])
starting.motifs <- RegNet

if(data.source=='file'){
  
  Idx1=match(gold.standard[,1], TFNames);
  Idx2=match(gold.standard[,2], GeneNames);
  Idx=(Idx2-1)*NumTFs+Idx1;
  gold.net=matrix(data=0, NumTFs, NumGenes);
  gold.net[Idx]=as.numeric(gold.standard[,3])
}

# PPI Data
TFCoop=diag(NumTFs);
Idx1=match(PPI[,1], TFNames);
Idx2=match(PPI[,2], TFNames);
Idx=(Idx2-1)*NumTFs+Idx1;
TFCoop[Idx]=as.numeric(PPI[,3]);
Idx=(Idx1-1)*NumTFs+Idx2;
TFCoop[Idx]=as.numeric(PPI[,3]);



## Run PANDA ##
tic=proc.time()[3];

#RegNet <- RegNet[-51,]
#TFCoop <- TFCoop[-51,-51]
RegNet <- RegNet
TFCoop <- TFCoop
print('Normalizing Networks!');
RegNet=NormalizeNetwork(RegNet);
TFCoop=NormalizeNetwork(TFCoop);
GeneCoReg=NormalizeNetwork(GeneCoReg);

print('Leaning Network!')
step=0;
hamming=1;
if (updateMethod=='danmethod'){
  print("using danmethod similarity")
#   while((hamming>0.00001) && (step<200))
#   {
  	Responsibility=Dfunction(TFCoop, RegNet);
  	Availability=Dfunction(t(RegNet), GeneCoReg);
  	hamming=sum(abs(RegNet-0.5*(Responsibility+Availability)))/(NumTFs*NumGenes);
  	RegNet=(1-alpha)*RegNet+alpha*0.5*(Responsibility+Availability);
  
  	PPI=Dfunction(t(RegNet), t(RegNet));
#   	PPI=UpdateDiagonal(PPI, NumTFs, alpha, step);
  	TFCoop=(1-alpha)*TFCoop+alpha*PPI;
  
  	CoReg2=cor(RegNet, RegNet);
  	CoReg2=UpdateDiagonal(CoReg2, NumGenes, alpha, step);
  	GeneCoReg=(1-alpha)*GeneCoReg+alpha*CoReg2;
  	print(paste("Step #", step, ", hamming =", round(hamming,5)), sep="");
          step=step+1;
#   }
} else if (updateMethod=='tanimoto'){
  print("using tanimoto similarity")
  while(hamming>0.00001)
  {
    Responsibility=Tfunction(TFCoop, RegNet);
    Availability=Tfunction(RegNet, GeneCoReg);
    hamming=sum(abs(RegNet-0.5*(Responsibility+Availability)))/(NumTFs*NumGenes);
    RegNet=(1-alpha)*RegNet+alpha*0.5*(Responsibility+Availability);
    
    PPI=Tfunction(RegNet, t(RegNet));
    PPI=UpdateDiagonal(PPI, NumTFs, alpha, step);
    TFCoop=(1-alpha)*TFCoop+alpha*PPI;
    
    CoReg2=Tfunction(t(RegNet), RegNet);
    CoReg2=UpdateDiagonal(CoReg2, NumGenes, alpha, step);
    GeneCoReg=(1-alpha)*GeneCoReg+alpha*CoReg2;
    print(paste("Step #", step, ", hamming =", round(hamming,5)), sep="");
    step=step+1;
  }
}
toc=proc.time()[3] - tic;
print(paste("Running PANDA on ", NumGenes, " Genes and ", NumTFs, " TFs took ", round(toc,2), " seconds!", sep=""));

library(pROC)
library(gplots)


# # ROC curves
# if (data.source=='file'){
#   par(mfrow=c(2,1)) 
#   print(plot.roc(gold.net,RegNet, col="blue"))
#   print(lines.roc(gold.net,starting.motifs,col="black"))
#   print(plot.roc(gold.net[starting.motifs==1],RegNet[starting.motifs==1], col="blue"))
#   print(lines.roc(gold.net[starting.motifs==0],RegNet[starting.motifs==0], col="blue"))
#   
# }

#RegNet.null1 <- RegNet
#RegNet.null2 <- RegNet
#RegNet.SR <- RegNet
#RegNet.CC <- RegNet
RegNet.KO <- RegNet


RegNet.SR.t <- t(RegNet.SR) 
RegNet.CC.t <- t(RegNet.CC)
RegNet.KO.t <- t(RegNet.KO)

## Find the linear transformation from SR to CC

net1 <- RegNet.SR
net2 <- RegNet.CC
gene.trans.matrix <- svd(net2)$v %*% diag(1/svd(net2)$d) %*% t(svd(net2)$u) %*% net1
tf.trans.matrix   <- svd(t(net2))$v %*% diag(1/svd(t(net2))$d) %*% t(svd(t(net2))$u) %*% t(net1)
tf.trans.matrix.norm <- apply(tf.trans.matrix, 1, function(x){
  #   x.zero <- (x-mean(x))
  x/sum(abs(x))
})
heatmap(tf.trans.matrix.norm)
heatmap.2(tf.trans.matrix.norm, dendrogram = "none", Rowv = FALSE, Colv = FALSE) 
tf.trans.matrix.norm.nondiag <- tf.trans.matrix.norm
diag(tf.trans.matrix.norm.nondiag) <- 0
heatmap.2(tf.trans.matrix.norm.nondiag, dendrogram = "none", Rowv = FALSE, Colv = FALSE) 
heatmap.2(tf.trans.matrix.norm.nondiag, dendrogram = "both",col=bluered) 


gene.trans.matrix.norm <- apply(gene.trans.matrix, 1, function(x){
  x/sum(abs(x))
})
heatmap(gene.trans.matrix.norm)
heatmap.2(gene.trans.matrix.norm , dendrogram = "none", Rowv = FALSE, Colv = FALSE) 
gene.trans.matrix.norm.nondiag <- gene.trans.matrix.norm
diag(gene.trans.matrix.norm.nondiag) <- 0
heatmap.2(gene.trans.matrix.norm.nondiag, dendrogram = "none", Rowv = FALSE, Colv = FALSE) 
heatmap.2(gene.trans.matrix.norm.nondiag, dendrogram = "both") 

hclust(tf.trans.matrix)

# Get off-diagonal mass
off.diags <- apply(tf.trans.matrix.norm.nondiag,1,function(x){sum(abs(x))})
cbind(off.diags[order(off.diags)],as.character(TFNames[order(off.diags)]))

max.off.diagonal <- apply(tf.trans.matrix.norm.nondiag,1,function(x){
  TFNames[which.max(x)]
})

predicted.interactions <- cbind(as.character(TFNames),as.character(max.off.diagonal))
predicted.interactions[20,]


# 12/16/14
# PCA
tf.pca <- prcomp(tf.trans.matrix.norm.nondiag, center = TRUE, scale. = TRUE)

sample.groups <- c(rep(0,105))
scores <- data.frame(sample.groups, tf.pca$x[,1:3])
pc1.2 <- qplot(x=PC1, y=PC2, data=scores, colour=factor(sample.groups),label=rownames(regnet1), geom="text") +
  theme(legend.position="none")
pc1.3 <- qplot(x=PC1, y=PC3, data=scores, colour=factor(sample.groups)) +
  theme(legend.position="none")
pc2.3 <- qplot(x=PC2, y=PC3, data=scores, colour=factor(sample.groups)) +
  theme(legend.position="none")
plot(pc1.2)


library(devtools)
install_github("ggbiplot", "vqv")





library(ggbiplot)
g <- ggbiplot(tf.pca, obs.scale = 1, var.scale = 1, 
              groups = ir.species, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)