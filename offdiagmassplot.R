
#####################################################################################
##  This script randomizes each gene within each condition, effectively removing the
##  correlated structure in the expression data.  It then builds a specified number of
##  "null" transition matrices and compares the distribution of off-diagonal masses to
##  that of the observed transition matrix.
######################################################################################

library(ggplot2)
library(reshape2)

num.iterations <- 30

# Generate a number of "null" networks from the first dataset
null.trans.matrices.A <- lapply(1:num.iterations,FUN=function(x){
  panda(yeast$motif, yeast$exp.cc, yeast$ppi, randomize="by.genes", progress=T)$reg.net[-8,]
})

# Generate the same number of "null" networks from the second dataset
null.trans.matrices.B <- lapply(1:num.iterations,FUN=function(x){
  panda(yeast$motif, yeast$exp.sr, yeast$ppi, randomize="by.genes", progress=T)$reg.net[-8,]
})

#store results from null permutations
yeast.null.networks <- list(null.trans.matrices.A,null.trans.matrices.B)
setwd("C:/Users/Dan/Dropbox/Dan's/Documents/Harvard/Research")
save(yeast.null.networks,file="yeast.null.rda")
load(file = "yeast.null.rda")

# Calculate the transition matrices between those "null" networks
tm.null <- lapply(1:num.iterations, function(x){
  transformation.matrix(null.trans.matrices.A[[x]], null.trans.matrices.B[[x]])
})

# Start of PandaAnalysis function
# Calculate the off-diagonal squared mass for each transition matrix
null.SSODM <- lapply(tm.null,function(x){
  apply(x,1,function(y){t(y)%*%y})
})
null.ssodm.matrix <- matrix(unlist(null.SSODM),ncol=num.iterations)
null.ssodm.matrix <- t(apply(null.ssodm.matrix,1,sort))

# Process the data for ggplot2
null.SSODM.melt <- melt(null.SSODM)
null.SSODM.melt$TF<-rep(rownames(tm.null[[1]]),num.iterations)
null.SSODM.melt$data.source<-"null"

# Combine it with the observed values
# Must run panda to get panda.res...  for this part if not done already
tm <- transformation.matrix(panda.res.cc$reg.net[-8,], panda.res.sr$reg.net[-8,],remove.diagonal=T)
ssodm <- apply(tm,1,function(x){t(x)%*%x})
observed.SSODM <- cbind(value=ssodm,L1=0,TF=rownames(tm),data.source="observed")
null.SSODM.melt <- rbind(null.SSODM.melt,observed.SSODM)
null.SSODM.melt[,1] <- as.numeric(null.SSODM.melt[,1])

# Get p-value (rank of observed within null ssodm)
p.values <- sapply(1:length(ssodm),function(i){
  1-findInterval(ssodm[i], null.ssodm.matrix[i,])/num.iterations
})

## Plot the data
ggplot(null.SSODM.melt, aes(x=TF, y=value))+ 
  geom_point(aes(size=1,color=factor(data.source),alpha = .5 + .5*as.numeric(factor(data.source)))) + 
  scale_color_manual(values = c("blue", "red")) +
  scale_x_discrete(limits = rownames(tm.null[[1]]) ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=pmin(15,7-log(p.values)),face="bold")) + 
  ylab("Sum of Squared Off-Diagonal Mass") +
  ggtitle("SSODM observed and null")
