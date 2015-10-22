# Example transition
P <- matrix(c(1,1,1,0,0,0,0,0,0, 0,1,1,1,1,0,0,0,0, 0,0,0,0,1,1,1,0,0,
                 0,0,0,0,0,0,1,1,1, 0,0,0,0,0,0,0,1,1),nrow=9)
Q <- matrix(c(1,1,1,0,0,0,0,0,0, 0,1,1,1,1,0,0,0,0, 0,0,0,0,1,1,1,0,0,
                 0,0,0,0,1,1,1,0,0, 0,0,0,0,0,0,0,1,1),nrow=9)
P
Q
W.kabsch <- kabsch(P,Q)
W.old <- transformation.matrix(P,Q,by.tfs=F,remove.diagonal=F,method="old")

crossprod(c(P%*%W.kabsch - Q))
crossprod(c(P%*%W.old - Q))

P<- matrix(rnorm(24),nrow=6)
Q<- matrix(rnorm(24),nrow=6)

Q<- t(apply(matrix(rnorm(16),nrow=4),1,function(x){
  x-mean(x)
}))

Q<- scale(matrix(rnorm(8),nrow=2))
W.kabsch <- kabsch(P,Q)
W.old    <- transformation.matrix(P,Q,remove.diagonal=F,by.tfs=F,method='old')
W.solve  <- solve(P,Q)
Q.est.k <- P%*%W.kabsch
Q.est.o <- P%*%W.old
round(P,3)
round(Q,3)

plot(P,Q%*%W.kabsch)
plot(Q,P%*%W.kabsch)
plot(W.kabsch%*%P,Q)
plot(W.kabsch%*%Q,P)
plot(P%*%W.solve,Q)

plot(P,Q%*%W.old)
plot(Q,P%*%W.old)
plot(W.old%*%P,Q)
plot(W.old%*%Q,P)
summary(lm(c(Q)~c(Q.est.k)))

transformation.matrix(diag(5),A,remove.diagonal=F)
kabsch(A,B) %*% B

# http://en.wikipedia.org/wiki/Kabsch_algorithm
# http://math.stackexchange.com/questions/77462/finding-transformation-matrix-between-two-2d-coordinate-frames-pixel-plane-to-w

kabsch <- function(P,Q){
#   P <- scale(P)
#   Q <- scale(Q)
  covmat <- t(Q) %*% P
  num.TFs <- ncol(covmat)
  svd.res <- svd(covmat)
  
  # Note the scalar multiplier in the middle.
  # NOT A MISTAKE!
  c.k <- colSums(P %*% svd.res$v * Q %*% svd.res$u)
  
  E <- diag(num.TFs)
  
  W <- svd.res$v %*% E %*% t(svd.res$u)
  W
}

rot <- function(P,Q){
  print(paste(sep="","nrow=",nrow(P),", ncol=",ncol(P)))
  covmat<-cov(t(P),t(Q))
  P%*%t(Q)
  svd(covmat)
  svd(P%*%t(Q))
  E<- diag(svd(covmat)$d)
  R <- svd(covmat)$v%*%E%*%t(svd(covmat)$u)
  plot(R%*%P,Q)
  plot(t(P))
  points(t(R%*%P),pch="x")
  summary(lm(c(Q)~c(R%*%P)))
}
euclid <- function(A,B){
  sum(apply(A - B,1,function(x){
    x%*%x
  }))
}

rot(P,Q)

a<-matrix(rnorm(16),nrow=4)
a <- t(a)%*%a
a.eig <- eigen(a)
a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)

transMat <- function(P,Q){
  P <- apply(P,2,function(x){
    x - mean(x)
  })
  Q <- apply(Q,2,function(x){
    x - mean(x)
  })
  A <- t(P) %*% Q
  A <- cov(P,Q)
  a.eigen <- eigen(t(A)%*%A)
  a.sqrt <- a.eigen$vectors %*% diag(sqrt(a.eigen$values)) %*% solve(a.eigen$vectors)
  U <- a.sqrt %*% solve(A)
  U
}
P
Q
U <- transMat(P,Q)
t(U%*%t(P))




library(ggplot2)
P<-matrix(c(1,2,3,1,2,3),nrow=3)
Q<-matrix(c(1,2,3,-1,-2,-3),nrow=3)
seed<-matrix(rnorm(100),nrow=50)
P<-seed+matrix(rnorm(100),nrow=50)
Q<-seed+matrix(rnorm(100),nrow=50)
P<-matrix(c(1,2,1,2.1),nrow=2)
Q<-matrix(c(1,2.3,-.5,-1),nrow=2)
P<-matrix(c(rnorm(50,50),rnorm(50,25)),nrow=50)
Q<-matrix(c(rnorm(50,25),rnorm(50,50)),nrow=50)

Q.hat <- P%*%W
Q.hat <- P%*%kabsch(P,Q)
Q.hat <- P%*%transformation.matrix(P,Q,by.tfs=F,method='old',standardize=F,remove.diagonal=F)
beta.matrix <- solve(t(P)%*%P)%*%t(P)%*%Q
Q.hat <- P%*%beta.matrix

plot.transformation(P,Q,Q.hat)

plot.transformation <- function(P,Q,Q.hat){
  df <- data.frame(rbind(P,Q.hat,Q))
  df$net <- c(rep("P",nrow(P)),rep("Q.hat",nrow(Q.hat)),rep("Q",nrow(Q)))
  df$tf  <- rep(1:nrow(P),3) 
  ggplot(df, aes(X1, X2, group = tf)) +   
    geom_point(size=8,aes(color=factor(net),shape = factor(net))) +
    geom_line() +
    geom_text(data = df,aes(x=X1,y=X2, label=tf))
}


beta.matrix

nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
nba$Name <- with(nba, reorder(Name, PTS))
nba.m <- melt(nba)
nba.m <- ddply(nba.m, .(variable), transform,rescale = rescale(value))

melt.tm <- melt(lasso)
ggplot(melt.tm, aes(y=Var1,x=Var2)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="white", high="darkblue") + 
  xlab("") + 
  ylab("") +
  theme(axis.text.x = element_text(angle = 90,hjust=1))


apply(tm.observed,2,)



nn.adj.mat.COPD <- eclipse$COPD.network
nn.adj.mat.COPD[nn.adj.mat.COPD<0] <-0

nn.adj.mat.SmCo <- eclipse$SmCo.network
nn.adj.mat.SmCo[nn.adj.mat.SmCo<0] <-0

nn.tm <- transformation.matrix(nn.adj.mat.COPD, nn.adj.mat.SmCo,remove.diagonal=F,method="ols")
mean(nn.tm)*189

tm.cor.mat <- cor(t(eclipse$COPD.network),t(eclipse$SmCo.network))
tm.cor.mat.melt <- melt(tm.cor.mat)

ggplot(tm.cor.mat.melt, aes(y=Var1,x=Var2)) + 
  geom_tile(aes(fill=value)) + 
  scale_fill_gradient(low="white", high="darkblue") + 
  xlab("") + 
  ylab("") +
  theme(axis.text.x = element_text(angle = 90,hjust=1))

rank.tm.cor.mat <- apply(tm.cor.mat,2,rank)

vals <- c(rnorm(50,1),rnorm(50))
groups <- c(rep("A",50),rep("B",50))
predict(lda(groups~vals),data.frame(vals))

bereLDA <- function(motif.data, 
                            expr.data,
                            verbose=T,
                            randomize="none",
                            cpp=T){
  if(verbose)
    print('Initializing and validating')
  # Create vectors for TF names and Gene names from Motif dataset
  tf.names   <- sort(unique(motif.data[,1]))
  num.TFs    <- length(tf.names)
  if (is.null(expr.data)){
    stop("Error: Expression data null")
  } else {
    # Use the motif data AND the expr data (if provided) for the gene list
    gene.names <- sort(intersect(motif.data[,2],rownames(expr.data)))
    num.genes  <- length(gene.names)
    
    # Filter out the expr genes without motif data
    expr.data <- expr.data[rownames(expr.data) %in% gene.names,]
    
    # Keep everything sorted alphabetically
    expr.data      <- expr.data[order(rownames(expr.data)),]
    num.conditions <- ncol(expr.data);
    if (randomize=='within.gene'){
      expr.data <- t(apply(expr.data, 1, sample))
      if(verbose)
        print("Randomizing by reordering each gene's expression")
    } else if (randomize=='by.genes'){
      rownames(expr.data) <- sample(rownames(expr.data))
      expr.data           <- expr.data[order(rownames(expr.data)),]
      if(verbose)
        print("Randomizing by reordering each gene labels")
    }
  }
  
  # Bad data checking
  if (num.genes==0){
    stop("Error validating data.  No matched genes.\n  Please ensure that gene names in expression file match gene names in motif file.")
  }
  
  strt<-Sys.time()
  if(num.conditions==0) {
    stop("Error: Number of samples = 0")
    gene.coreg <- diag(num.genes)
  } else if(num.conditions<3) {
    stop('Not enough expression conditions detected to calculate correlation.')
  } else {
    if(verbose)
      print('Verified adequate samples, calculating correlation matrix')
    if(cpp){
      # C++ implementation
      gene.coreg <- rcpp_ccorr(t(apply(expr.data, 1, function(x)(x-mean(x))/(sd(x)))))
      rownames(gene.coreg)<- rownames(expr.data)
      colnames(gene.coreg)<- rownames(expr.data)
      
    } else {
      # Standard r correlation calculation
      gene.coreg <- cor(t(expr.data), method="pearson", use="pairwise.complete.obs")
    }
  }
  
  print(Sys.time()-strt)
  
  if(verbose)
    print('More data cleaning')
  # Convert 3 column format to matrix format
  colnames(motif.data) <- c('TF','GENE','value')
  regulatory.network <- tidyr::spread(motif.data, GENE, value, fill=0)
  rownames(regulatory.network) <- regulatory.network[,1]
  # sort the TFs (rows), and remove redundant first column
  regulatory.network <- regulatory.network[order(rownames(regulatory.network)),-1]
  # sort the genes (columns)
  regulatory.network <- as.matrix(regulatory.network[,order(colnames(regulatory.network))])
  
  # Filter out any motifs that are not in expr dataset (if given)
  if (!is.null(expr.data)){
    regulatory.network <- regulatory.network[,colnames(regulatory.network) %in% gene.names]
  }
  
  # store initial motif network (alphabetized for rows and columns)
  #   starting.motifs <- regulatory.network
  
  
  if(verbose)
    print('Main calculation')
  ########################################
  
  strt<-Sys.time()
  result <- apply(regulatory.network,1,function(x){
    sapply(1:nrow(gene.coreg),function(i){
     predict(lda(x~gene.coreg[i,]),data.frame(gene.coreg[i,]))$x
    })
  })
#   correlation.dif <- sweep(regulatory.network,1,rowSums(regulatory.network),`/`)%*%gene.coreg-sweep(1-regulatory.network,1,rowSums(1-regulatory.network),`/`)%*%gene.coreg
#   result <- sweep(correlation.dif, 2, apply(correlation.dif, 2, sd),'/')
  #   regulatory.network <- ifelse(res>quantile(res,1-mean(regulatory.network)),1,0)
  
  print(Sys.time()-strt)
  ########################################
  
  return(result)
}
n<-5;m<-5
sum(sapply(0:min(m,n), function(x){choose(m+n-2*x,n-x)*choose(m+n-x,x)}))


##  Start of lasso for transition matrix
library(genlasso)
lassoTM <- function(net1,net2, unpenalized=F){
  lasso.res <- t(sapply(1:53, function(x){
    D<-diag(53)
    if (unpenalized){
      D[x,x]<-0
    }
    coef(genlasso(net1[x,],t(net2),D=D))$beta[,6]
  }))
  colnames(lasso.res) <- rownames(net1)
  rownames(lasso.res) <- rownames(net2)
  lasso.res
}
OLSTM <- function(net1,net2, unpenalized=F){
  tf.trans.matrix <- ginv(net1%*%t(net1))%*%net1%*%t(net2)
  colnames(tf.trans.matrix) <- rownames(net1)
  rownames(tf.trans.matrix) <- rownames(net2)
  tf.trans.matrix
}

tm.old <- OLSTM(bere.res.cc,bere.res.ko)

library(ggplot2)

melt.tm.ccko.wmotif.old <- melt(transformation.matrix(bere.res.cc.wmotif,bere.res.ko.wmotif,method="old",standardize=T,remove.diagonal=T))
melt.tm.ccko.wmotif <- melt(lassoTM(bere.res.cc.wmotif,bere.res.ko.wmotif))
melt.tm.ccko.wmotif.ols <- melt(OLSTM(bere.res.cc.wmotif,bere.res.ko.wmotif))
melt.tm.ccko <- melt(lassoTM(bere.res.cc,bere.res.ko))
melt.tm.ccko.ols <- melt(OLSTM(bere.res.cc,bere.res.ko))

plottm(melt.tm.ccko.wmotif, title="Lasso TM CC vs KO, motif readded")
plottm(melt.tm.ccko.wmotif.ols, title="OLS TM CC vs KO, motif readded")
plottm(melt.tm.ccko.wmotif.old, title="SVD method TM CC vs KO, motif readded")

plottm(melt.tm.ccko, title="Lasso TM CC vs KO, nomotif readded")
plottm(melt.tm.ccko.ols, title="OLS TM CC vs KO, nomotif readded")

melt.tm.COPD_SMC <- melt(tm.observed)
plottm(melt.tm.COPD_SMC, title="OLS TM COPD vs SMC, nomotif readded")

plottm <- function(melt.tm,title="Transition Matrix"){
  ggplot(melt.tm, aes(y=Var1,x=Var2)) + 
    ggtitle(title) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient(low="white", high="darkblue") + 
    xlab("") + 
    ylab("") +
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))
}

plotmasses <- function(melt.tm, title="Mass distribution",xmin=0,xmax=1,gausscurve=F){
    ggplot(melt.tm, aes(value)) + 
        ggtitle(title) + 
        xlim(xmin,xmax) +
        geom_histogram(colour="black", fill="white",binwidth=(xmax-xmin)/100,aes(y = ..density..)) +
        stat_function(fun=dnorm,size=as.numeric(gausscurve), args=list(mean=mean(melt.tm$value), sd=sd(melt.tm$value)),colour = 'red')    
}

plottm(melt.tm.COPD_COPD)
plottm(melt.tm.COPD_SMC)

melt.tm.null <- melt(tm.null)
plottm(melt.tm.null)

plotmasses(melt.tm.COPD_COPD)
plotmasses(melt.tm.COPD_SMC)

plotmasses(melt.tm.null, gausscurve=F)
melt.tm.null.ODM <- removeDiagonal(melt.tm.null)
plotmasses(melt.tm.null.ODM, gausscurve=T, xmin=-.02,xmax=.02, title="Mass distribution for off-diagonal NULL")

melt.tm.COPD_SMC.ODM <- melt.tm.COPD_SMC[melt.tm.COPD_SMC[,1]!=melt.tm.COPD_SMC[,2],]
plotmasses(melt.tm.COPD_SMC.ODM, gausscurve=T, xmin=-.02,xmax=.02, title="Mass distribution for off-diagonal COPD vs Smoker Control")


#
removeDiagonal <- function(x){
    x[x[,1]!=x[,2],]
}

## creating null networks for yeast
null.exp <- cbind(yeast$exp.ko, yeast$exp.cc)
tm.nulls <- lapply(1:200,function(x){
  group <- c(rep("A",106),rep("B",50))
  rownames(null.exp) <- rownames(null.exp)[sample(1:nrow(null.exp))]
  
  null.A.wmotif <- bere(yeast$motif, null.exp[,group=="A"], cpp=F)
  null.B.wmotif <- bere(yeast$motif, null.exp[,group=="B"], cpp=F)
  res <- OLSTM(null.A.wmotif,null.B.wmotif)
  diag(res)<-0  
  res
})
tm.observed <- OLSTM(bere.res.cc.wmotif,bere.res.ko.wmotif)
diag(tm.observed) <- 0
# Do the sum of sq ODM plot versus null
ssodm.plot(tm.observed, tm.nulls, plot.title="SSODM observed and null, Yeast CC vs KO")
ssodm.plot(tm.observed, tm.nulls, rescale=T, plot.title="SSODM observed and null, Yeast CC vs KO",highlight.tfs = c("ELK1","E2F4"))

dim(yeast$ppi)

ppi.table <- regulatory.network
ppi.table <- cbind(ppi.table,0)
ppi.table <- rbind(ppi.table,0)
rownames(ppi.table)[43] <- "YPR104C"
colnames(ppi.table)[43] <- "YBL005W"
ppi.table <- ppi.table[sort(row.names(ppi.table)),sort(colnames(ppi.table))]
small.tm.observed <- tm.observed[sort(row.names(ppi.table)),sort(colnames(ppi.table))]
diag(ppi.table) <- 0
diag(small.tm.observed) <-0

onlydiag <- 1-diag(43)

summary(lm(c(ppi.table) ~ c(small.tm.observed) + c(onlydiag)))
ppi.auc <- plot.roc(ppi.table, small.tm.observed, col="blue", main="Yeast PPI predictions")$auc
baselineppi.auc <- lines.roc(ppi.table, onlydiag, col="red")$auc
legend(.6,.2,
       c(paste("PPI predictions",round(ppi.auc,4)),paste("Baseline (diag only) ",round(baselineppi.auc,4))),
       lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("blue","red","grey"))



## Create PPI ROC curve for COPD data (4/9/15)
colnames(eclipse$ppi) <- colnames(melt.tm.COPD_SMC.ODM)
ppi.table <- merge(melt.tm.COPD_SMC.ODM, eclipse$ppi, by=c("Var1","Var2"),all.x=T)
ppi.table[is.na(ppi.table)] <- 0
makeROCPlot(ppi.table,"COPD PPI predictions")

ppi.table.null <- merge(melt.tm.null.ODM, eclipse$ppi, by=c("Var1","Var2"),all.x=T)
ppi.table.null[is.na(ppi.table.null)] <- 0
makeROCPlot(ppi.table.null,"NULL PPI predictions")
library(pROC)
makeROCPlot <- function(x, title="PPI predictions"){
    pVal <- t.test(x[x[,4]==1,3], x[x[,4]==0,3])$p.value
    ppi.auc <- plot.roc(x[,4], x[,3], col="blue", main=title)$auc
    legend(.8,.1,
           c(paste("ROC PPI predictions",round(ppi.auc,4),"; p=",round(pVal,4))),
           lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","grey"))
    
}


# lda
ldaExp <- data.frame(eclipse$exp)
tfdcast <- dcast(eclipse$motif,V1~V2,fill=0)
rownames(tfdcast) <- tfdcast[,1]
tfdcast <- tfdcast[,-1]

expData <- eclipse$exp[sort(rownames(eclipse$exp)),]
# check that IDs match
if (prod(rownames(expData)==colnames(tfdcast))!=1){
    stop("ID mismatch")
}

ldaBERE <- function(tfdcast,expData){
    t(apply(tfdcast, 1, function(x){
        cat(".")
        tfTargets <- as.numeric(x)
        z <- lda(tfTargets ~ ., expData)
    #     ldaRes <- cbind(predict(z, expData)$class, predict(z, expData)$posterior,tf1)
        predict(z, expData)$posterior[,2]
    }))
}
ldaPredictions.COPD <- ldaBERE(tfdcast, expData[,filter.vec.1])
ldaPredictions.SMC  <- ldaBERE(tfdcast, expData[,filter.vec.2])

berePredictions <- bere(eclipse$motif, eclipse$exp, cpp=F, score="none")
data <- data.frame(c(c(berePredictions),c(ldaPredictions)))
data$method <- c(rep("bere", nrow(data)/2),rep("lda", nrow(data)/2))
colnames(data) <- c("edge","method")
ggplot(data, aes(x=edge, fill=method)) + ggtitle("Edgeweight distribution") +
    geom_histogram(binwidth=.01, alpha=.5, position="identity") + 
    xlim(0, 3)

#################################
# Validation of TM in DREAM data
# 4/26/15
#################################

# Get BERE for with/out gene upregulated

library(bptools)
dataset <- "DREAM5c"
source(paste("./",dataset,".R",sep=""))

bereMelt_background <- bere(motifs, exprData[, -grep("G313",chipFeatures$V5)], cpp=F, verbose=F)

bereMelt_G313 <- bere(motifs, exprData[, grep("G313",chipFeatures$V5)], cpp=F, verbose=F)

tm.G313.obs <- transformation.matrix(bereMelt_background, bereMelt_G313, remove.diagonal=F,by.tfs=T)


# Copy expression data for null network generation
null.exp <- exprData
tm.G313.null <- lapply(1:200, function(iteration){
    rownames(null.exp) <- rownames(null.exp)[sample(1:nrow(null.exp))]
    res1 <- bere(motifs, null.exp[,-grep("G2388",chipFeatures$V6)], cpp=F)
    res2 <- bere(motifs, null.exp[,grep("G2388",chipFeatures$V6)], cpp=F)
    list(res1,res2) 
})

# Calculate the transformation matrix for the null data
tm.G313.null <- lapply(tm.G313.null, function(x){
    transformation.matrix(x[[1]],x[[2]],method="ols",remove.diagonal = F)
})

# Do the sum of sq ODM plot versus null
ssodm.plot(tm.G313.obs, tm.G313.null, plot.title="Gene 313 over expressed vs Not")
ssodm.plot(tm.G313.obs, tm.G313.null, rescale=T, plot.title="Gene 313 over expressed vs Not", highlight.tfs = c("G313"))

plottm(melt(transformation.matrix(bereMelt_background, bereMelt_G313, remove.diagonal=T,by.tfs=T)))

# Networks for Kimbie's boolean
netForKG1 <- read.table("~/network_edges_1",stringsAsFactors=F)
netForKG1 <- read.table("~/network_edges_2",stringsAsFactors=F)
netForKG1 <- dcast(netForKG1,V1~V2, fill=0)
rownames(netForKG1) <- netForKG1[,1]
netForKG1 <- netForKG1[,-1]
netForKG1 <- netForKG1[,colnames(netForKG1) %in% rownames(netForKG1)]
for (i in 1:114){
    if (sum(netForKG1[,i])==0){
        netForKG1[sample(1:114,1),i] <-1
    }
    if (sum(netForKG1[i,])==0){
        netForKG1[i,sample(1:114,1)] <-1
    }
}
netForKG1 <- removeZeroRowsColumns(netForKG1)
dim(netForKG1)
net.melt <- melt(as.matrix(netForKG1))
write.table(net.melt, file="network2_fixed", row.names = F, col.names = F, quote=F)

removeZeroRowsColumns <- function(x){
    x <- x[,colnames(x) %in% rownames(x)]
    x <- x[rownames(x) %in% colnames(x),]
    x[apply(x,1,sum)>0, apply(x,2,sum)>0]
}

###### 5/9/15
## Generating TM network

require(igraph)

tm.sigmas <- transitionPValues(tm.observed, tm.null)
diag(tm.sigmas) <- 0
tm.sigmas.melt <- melt(tm.sigmas)

adjMat <- tm.observed
diag(adjMat) <- 0
adjMat.melt <- melt(adjMat)

adj.combined <- merge(tm.sigmas.melt,adjMat.melt, by=c("Var1","Var2"))
# adjMat.melt <- cbind(adjMat.melt,abs(adjMat.melt[,3])>25)
adj.combined <- adj.combined[abs(adj.combined[,4])>.012,]
tfNet <- graph.data.frame(adj.combined, directed=T)
vSize <- t.values
vSize[vSize<0] <- 0
vSize <- sapply(vSize*3, min, 50)

V(tfNet)$size <- vSize[V(tfNet)$name]
E(tfNet)$width <- (abs(E(tfNet)$value.x)-2)*2
E(tfNet)$color<-ifelse(E(tfNet)$value.x>0, "blue", "red")
plot.igraph(tfNet, edge.arrow.size=1, vertex.label.cex= 1.5, vertex.label.color= "black",main="Transition: SMC -> COPD")
legend(-1.7,1.3, c("Gained features","Lost features"), lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"))

## Calculate p-values for off-diagonals
transitionPValues <- function(tm.observed, tm.null){
    tm.null.mean <- apply(simplify2array(tm.null), 1:2, mean)
    tm.null.sd <- apply(simplify2array(tm.null), 1:2, sd)
    sigmas <- (tm.observed - tm.null.mean)/tm.null.sd
}


plottm(melt(tm.observed), title="Transition matrix: COPD vs Smoker-control (SMC) observed")
plottm(melt(tm.null[[1]]), title="Transition matrix: COPD vs Smoker-control (SMC) null1")
plottm(melt(tm.null[[2]]), title="Transition matrix: COPD vs Smoker-control (SMC) null2")
plottm(melt(tm.null[[3]]), title="Transition matrix: COPD vs Smoker-control (SMC) null3")
plottm(melt(tm.null[[3]]), title="Transition matrix: COPD vs Smoker-control (SMC) null4")


tm.melt <- melt(tm.null[[4]])
tm.melt[,3][abs(tm.melt[,3])<.008]<-0
plottm(tm.melt, title="Transition matrix: COPD vs Smoker-control (SMC) null4")


### Top transition interactions

topInteractions <- tm.sigmas.melt[order(-abs(tm.sigmas.melt[,3])),]
topInteractions <- cbind(topInteractions,"GL"=ifelse(topInteractions[,3]>0,"Gain","Loss"))
topInteractions <- cbind(topInteractions,"pVal"= 1-pnorm(abs(topInteractions[,3])))
topInteractions <- cbind(topInteractions, "FDR"=sapply(1:35721, function(x){
    topInteractions[x,5]*35721/x
}))
colnames(topInteractions) <- c("Var1","Var2","t-stat","Gain/Loss","p-value","FDR")
topInteractions <- merge(topInteractions,adjMat.melt, by=c("Var1","Var2"))
topInteractions <- topInteractions[order(topInteractions[,5]),]
topInteractions[1:20,]


adj.combined <- merge(topInteractions,adjMat.melt, by=c("Var1","Var2"))


## node degree only approach
degreeApproach <- function(motifs){
    tfDegree <- table(motifs[,c(1,3)])
    geneDegree <- table(motifs[,c(2,3)])
    
    tfMatrix <- matrix(rep(tfDegree[,2],nrow(geneDegree)), ncol=nrow(geneDegree))
    geneMatrix <- t(matrix(rep(geneDegree[,2],nrow(tfDegree)), nrow=nrow(geneDegree)))
    
    result <- tfMatrix+geneMatrix
    rownames(result) <- rownames(tfDegree)
    colnames(result) <- rownames(geneDegree)
    result
}

## Compare motif degree to gold standard degree
tfMotifDegree <- table(motifs[,c(1,3)])
tfGSDegree <- table(goldStandard[motifs[,3]==1,c(1,3)])


qplot(tfMotifDegree[,2],tfGSDegree[,2], geom=c("point", "smooth"), main="Motif priors degree vs Gold standard degree")


##  Homogeneity assessment
##  Attempt to measure if dTFI values are inflated due to homogeneity of groups compared to heterogeneity of null

permutations <- 3
resMatrix <- matrix(NA, nrow=permutations*4, ncol=164)
nullMatrix <- matrix(NA, nrow=permutations*2, ncol=164)
hetMatrix <- matrix(NA, nrow=permutations, ncol=164)
for (i in 0:(permutations-1)){
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
    hetMatrix[i+1,] <- apply(tm,2,function(x){t(x)%*%x})
    
}
colnames(resMatrix) <- colnames(tm)
colnames(nullMatrix) <- colnames(tm)

transcriptionFactor <- "E2F3"
values <- c(resMatrix[,transcriptionFactor], nullMatrix[,transcriptionFactor])
xvalues <- c(rep("obs",length(resMatrix[,transcriptionFactor])), rep("null",length(nullMatrix[,transcriptionFactor])))
qplot(y=values,x=xvalues, geom = "boxplot", main=paste("Observed vs Null for",transcriptionFactor))

hetValues <- sort(c(resMatrix)[c(F,T)])
homValues <- sort(c(nullMatrix))
qplot(y=hetValues, x=homValues, main=paste("QQ plot, homogenous networks vs hetergenous networks (all null)"))+ geom_abline(intercept = 0, slope = 1)



# 7/22/15 Plots for demo
samp <- sample(1:100)
geneX<-c(rnorm(50,0),rnorm(50,5))[samp]
geneY<-c(rnorm(50,0),rnorm(50,5))[samp]
col<-c(rep("red",50),rep("blue",50))
plot(geneX,geneY,main="Expr for two genes",col=col)
abline(lm(geneY[1:50]~geneX[1:50]), col="red" )
abline(lm(geneY[51:100]~geneX[51:100]), col="blue" )
legend(-2,7,       c("Null COPD","Null SMC"),       lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("blue","red"))


#8/10/15 Plot for comparison of two TM TF-TF transitions
qplot(c(tm.observed.COPGene[highlight.tfs,highlight.tfs]),c(tm.observed.ECLIPSE[highlight.tfs,highlight.tfs]), size=20, color=c(sapply(highlight.tfs, rep, 9)))


# 10/9/15 exploring properties of extreme p-values
library(MASS)
library(pcalg)

numPerms <- 10000
numSamples <- 10000
probGenos <- c(rep(.02,numSamples/2),rep(.01,numSamples/2))
probPhenos <- c(rep(.1,numSamples/2),rep(.02,numSamples/2))

suspectLoci <- rbinom(numSamples, size=1, prob=probGenos)
correctionVec <- c(rep(0,numSamples/2),rep(1,numSamples/2))
x <- cbind(1, correctionVec)
correctionHat  <- x %*% ginv(t(x)%*%x) %*% t(x)
suspectFitted <- correctionHat %*% suspectLoci
correctedSuspect <- suspectLoci - suspectFitted
suspectPheno <- replicate(numPerms, rbinom(numSamples, size=1, prob=probPhenos))
suspectPhenoFitted <- correctionHat %*% suspectPheno
correctedsuspectPheno <- suspectPheno - suspectPhenoFitted
correlations <- cor(correctedSuspect, correctedsuspectPheno)
uncorrectedCorrelations <- cor(suspectLoci, suspectPheno)

xProb <- suspectFitted/2
xVariances <- 2*xProb*(1-xProb) #variance of Binomial(2,xProb)
 yVariances <- suspectPhenoFitted*(1-suspectPhenoFitted)
# yVariances <- rbind(matrix(.02*.98,nrow=numSamples/2,ncol=numPerms),matrix(.09,nrow=numSamples/2,ncol=numPerms))
#yVariances <- matrix(probPhenos*(1-probPhenos),nrow=numSamples,ncol=numPerms)

zvar <- t(yVariances)%*%xVariances
numerator <- zvar
denominator <- colSums(yVariances)%*%t(colSums(xVariances))
varRsq <- numerator/denominator
varianceFactor <- 1/varRsq
which(is.nan(c(varianceFactor)))

nullnegLogPValue <- -log((1:numPerms)/numPerms)
negLogPValueUncorrected <- -log(1-pt(uncorrectedCorrelations*sqrt(c(numSamples-2)/c(1-uncorrectedCorrelations^2)), numSamples))
negLogPValue <- -log(1-pt(correlations*sqrt(c(varianceFactor-2)/c(1-correlations^2)), varianceFactor))
negLogPValueN <- -log(1-pt(correlations*sqrt(c(numSamples-2)/c(1-correlations^2)), numSamples))
plot(sort(c(nullnegLogPValue)), sort(c(negLogPValueN),na.last=F), col="red")
points(sort(c(nullnegLogPValue)), sort(c(negLogPValue),na.last=F), col="blue")
points(sort(c(nullnegLogPValue)), sort(c(negLogPValueUncorrected),na.last=F), col="black")
legend(4,2, c("Uncorrected","Corrected by pop","Corrected by pop and variance"), lty=1,lwd=5,col=c("black","red","blue"))
abline(0,1)


combinedMat<- cbind(suspectLoci, correctionVec, suspectPheno)
g.pvalue <- sapply(which(negLogPValueN>5 & negLogPValueN<5.3), function(i){
    -log(gSquareBin(1,i+2,2,combinedMat))
    })
g.pvalue
