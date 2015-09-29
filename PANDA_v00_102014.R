# Set Program Parameters
alpha=0.1;
exp_file="../Data/ToyExpressionData.txt"; #"../YeastData/YeastData_SRExpression.txt";
motif_file="../Data/ToyMotifData.txt"; #"../YeastData/YeastData_Motif.txt"
ppi_file="../Data/ToyPPIData.txt";

NormalizeNetwork<-function(X)
{
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

Exp=read.table(exp_file, sep="\t", row.names=1);
#Exp = gene_exp_sim
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

# Dan edit to use simuated data
Motif=read.table(motif_file, sep="\t", row.names=NULL);
#Motif = motif.data
TFNames=unique(Motif[,1]);
NumTFs=length(TFNames);
Idx1=match(Motif[,1], TFNames);
Idx2=match(Motif[,2], GeneNames);
Idx=(Idx2-1)*NumTFs+Idx1;
RegNet=matrix(data=0, NumTFs, NumGenes);
RegNet[Idx]=as.numeric(Motif[,3])
starting.RegNet <- RegNet
# PPI Data
TFCoop=diag(NumTFs);
if(ppi_file != "")
{
	PPI=read.table(ppi_file, sep="\t", row.names=NULL);
  # Dan using simulated data
  #PPI <- ppi.data
	Idx1=match(PPI[,1], TFNames);
	Idx2=match(PPI[,2], TFNames);
	Idx=(Idx2-1)*NumTFs+Idx1;
	TFCoop[Idx]=PPI[,3];
	Idx=(Idx1-1)*NumTFs+Idx2;
	TFCoop[Idx]=PPI[,3];
}


## Run PANDA ##
tic=proc.time()[3];

print('Normalizing Networks!');
RegNet=NormalizeNetwork(RegNet);
TFCoop=NormalizeNetwork(TFCoop);
GeneCoReg=NormalizeNetwork(GeneCoReg);

print('Leaning Network!')
step=0;
hamming=1;
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
toc=proc.time()[3] - tic;
print(paste("Running PANDA on ", NumGenes, " Genes and ", NumTFs, " TFs took ", round(toc,2), " seconds!", sep=""));
