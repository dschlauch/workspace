# Do the sum of sq ODM plot versus null
pdf(file.path(outputDir,paste0('SSODMplot_unscaled',analysisCode,'.pdf')), width=8)
sPlot <- ssodm.plot(transMatrices[[1]], transMatrices[-1],plot.title=paste("SSODM observed and null, ",casesString," vs ",controlsString,' : ', networkInferenceName, ' : ', analysisName, sep=""))
print(sPlot)
dev.off()
pdf(file.path(outputDir,paste0('SSODMplot_scaled',analysisCode,'.pdf')), width=24)
sPlot <- ssodm.plot(transMatrices[[1]], transMatrices[-1], rescale=T,plot.title="")#, plot.title=paste("SSODM observed and null, ",casesString," vs ",controlsString,' : ', networkInferenceName, ' : ', analysisName, sep=""))
print(sPlot)
dev.off()

# This is for converting TF IDs to their gene names

# transMatrices <- lapply(transMatrices, function(x){
#     rownames(x) <- mappings[match(rownames(x), mappings[,1]),2]
#     colnames(x) <- mappings[match(colnames(x), mappings[,1]),2]
#     x
# })

# Top TFs
#highlight.tfs <- c("E2F4","NRF1","GABPA","ELK1","ELK4","E2F1","ZBTB33","ELF1","ZFX")

# #####################################################
# ### Gene Analysis
# #####################################################
# # Calculate the transformation matrix for the observed data
# tm.observed.genes <- transformation.matrix(dataset$casesNetwork, dataset$controlsNetwork,remove.diagonal=T,method="old",by.tfs=T,standardize=F)
# 
# 
# # Calculate the transformation matrix for the null data
# tm.null.genes <- lapply(null.networks, function(x){
#   transformation.matrix(x[[1]],x[[2]],method="old",standardize=F)
# })
# 
# # Do the sum of sq ODM plot versus null
# ssodm.plot(tm.observed, tm.null,plot.title="SSODM observed and null, COPD subjects vs Smoker control (R)")
# ssodm.plot(tm.observed, tm.null, rescale=T, plot.title="SSODM observed and null, COPD subjects vs Smoker control (R)",highlight.tfs = c("ELK1","E2F4"))



######################################################
##  Running null networks with improved algorithm  ###
##                 2/25/15      END                ###
######################################################
# 
# 
# cppFunction('NumericMatrix correl(NumericMatrix x) {
#   int nrow = x.nrow(), ncol = x.ncol();
#   NumericMatrix resultMatrix(nrow,nrow);
#   for (int i = 0; i < nrow; i++) {
#     for (int j = i; j < nrow; j++) {
#       double sumproduct = 0;
#       for (int k = 0; k < ncol; k++){
#         sumproduct+=x(i,k)*x(j,k);
#       }
#       resultMatrix(i,j) = sumproduct/(nrow-1);
#       resultMatrix(j,i) = sumproduct/(nrow-1);
#     }
#   }
#   return resultMatrix;
# }')
# 
# cppFunction('double squarert(double x){
#   return sqrt(x);
#             }')
# set.seed(1014)
# randmat <- t(matrix(rnorm(30000*50), 50))
# t(apply(randmat, 1, function(x)(x-mean(x))/(sd(x))))
# 
# system.time(correl(t(apply(randmat, 1, function(x)(x-mean(x))/(sd(x))))))
# system.time(cor(t(randmat)))
# #>  [1] 458 558 488 458 536 537 488 491 508 528
# rowSumsC(x)
# #>  [1] 458 558 488 458 536 537 488 491 508 528
# 

## Gene expression analysis
design <- model.matrix(~factor(casesFilter))
diff.exp.res <- lmFit(dataset$exp, design)
diff.exp.res <- ebayes(diff.exp.res)
logfoldchange <- log(rowMeans(dataset$exp[,casesFilter])/rowMeans(dataset$exp[,controlsFilter]))


# 7/28/15 
# create results table
# 10/30/15 updates for GTEx, which does not (or I'm not using) expression values for TFs
obsSsodm <- apply(transMatrices[[1]],1,function(x){t(x)%*%x})
includedTFs <- intersect(names(obsSsodm),rownames(diff.exp.res$p.value))
obsSsodm <- obsSsodm[includedTFs]
dTFI_pVals_All <- 1-2*abs(.5-calculate.tm.p.values(transMatrices[[1]], transMatrices[-1],method="non-parametric"))
dTFI_normalized_scores_All <- 1-2*abs(.5-calculate.tm.p.values(transMatrices[[1]], transMatrices[-1],method="z-score"))
names(dTFI_pVals_All) <- colnames(transMatrices[[1]])
names(dTFI_normalized_scores_All) <- colnames(transMatrices[[1]])
dTFI_pVals <- dTFI_pVals_All[includedTFs]
dTFI_normalized_scores <- dTFI_normalized_scores_All[includedTFs]
negLogZPValues <- -log(dTFI_normalized_scores)
negLogPValues <- -log(dTFI_pVals)
# replace Inf values with max values
negLogPValues[negLogPValues==Inf] <- 35
negLogZPValues[negLogZPValues==Inf] <- 35
labels <- names(obsSsodm)
labels[rank(-negLogPValues)>20 & rank(-obsSsodm)>20]<-""

dTFI_fdr   <- p.adjust(dTFI_pVals, method = 'fdr')
logfoldchangeTF <- logfoldchange[names(dTFI_pVals)]
limma_pVals <- diff.exp.res$p.value[names(dTFI_pVals),2]
limma_fdr <- p.adjust(limma_pVals, method = 'fdr')
limmanegLogPValues <- -log(limma_pVals)
limmanegLogPValues[limmanegLogPValues>10]<-10 # for visual purposes
resultTable <- cbind(obsSsodm,dTFI_pVals, dTFI_normalized_scores, dTFI_fdr, logfoldchangeTF, limmanegLogPValues)
resultTable <- resultTable[order(dTFI_pVals),]

plotDF <- data.frame(obsSsodm, negLogPValues, limmanegLogPValues, logfoldchangeTF, negLogZPValues, "labels"=labels)

colnames(resultTable) <- c("Magnitude","dTFI uncorrected p-value","dTFI normalized scores","dTFI FDR", "log FC", "limma -logp")
write.csv(resultTable,file=file.path(outputDir,paste("resultTable",analysisCode,".csv", sep="")))

tiff(file.path(outputDir,paste('Volcano plot',analysisCode,'.tiff', sep="")), width=1800)
ggplot(data=plotDF,aes(x=obsSsodm, y=negLogZPValues, label=labels, size=100)) + geom_point(aes(col=limmanegLogPValues), alpha=.5) + geom_text(vjust=0) + 
  ylab("-log(dTFI p-value)") + xlab("SSODM") + ggtitle("Signal vs significance") + scale_size_continuous(range = c(0, 12),guide=FALSE) +
  scale_colour_gradientn("LIMMA sig",colours=c("blue","white","red"))
dev.off()

pdf(file.path(outputDir,paste('dTFI vs LIMMA',analysisCode,'.pdf', sep="")), width=9, height=8)
ggplot(data=plotDF, aes(x=limmanegLogPValues, y=negLogZPValues)) + geom_point(aes(col=logfoldchangeTF), size=5, alpha=.6) +
  geom_text_repel(data=plotDF[labels!="",], aes(limmanegLogPValues, negLogZPValues, label=labels)) + 
  ylab("Differential TF Involvement, -log(p-value)") + xlab("Differential Expression,  LIMMA -log(p-value)") + 
  ggtitle("Differential Involvement vs Differential Expression (ECLIPSE)") +
  theme_classic() + scale_colour_continuous(limits=c(-max(abs(logfoldchangeTF)),max(abs(logfoldchangeTF))), name="log(fold-change)", low = "red", high = "blue") +
  theme(plot.title = element_text(size=20,hjust=0))
dev.off()

# periodically save workspace
# save.image(file=file.path(outputDir,paste("activeImage",analysisCode,".RData",sep="")))

## Calculate p-values for off-diagonals
transitionSigmas <- function(tm.observed, tm.null){
  tm.null.mean <- apply(simplify2array(tm.null), 1:2, mean)
  tm.null.sd <- apply(simplify2array(tm.null), 1:2, sd)
  sigmas <- (tm.observed - tm.null.mean)/tm.null.sd
}

tm.sigmas <- transitionSigmas(transMatrices[[1]], transMatrices[-1])
diag(tm.sigmas) <- 0
tm.sigmas.melt <- melt(tm.sigmas)

adjMat <- transMatrices[[1]]
diag(adjMat) <- 0
adjMat.melt <- melt(adjMat)

adj.combined <- merge(tm.sigmas.melt,adjMat.melt, by=c("Var1","Var2"))

# adj.combined[,1] <- mappings[match(adj.combined[,1], mappings[,1]),2]
# adj.combined[,2] <- mappings[match(adj.combined[,2], mappings[,1]),2]

numEdges  <- 40
numTopTFs <- 10
topTFsIncluded <- names(sort(dTFI_pVals_All)[1:numTopTFs])
topTFIndices <- 2>(is.na(match(adj.combined[,1],topTFsIncluded))+is.na(match(adj.combined[,2],topTFsIncluded)))
adj.combined <- adj.combined[topTFIndices,]
adj.combined <- adj.combined[abs(adj.combined[,4])>=sort(abs(adj.combined[,4]),decreasing=T)[numEdges],]
tfNet <- graph.data.frame(adj.combined, directed=T)
vSize <- -log(dTFI_pVals_All)
vSize[vSize<0] <- 0
vSize[vSize>3] <- 3

V(tfNet)$size <- vSize[V(tfNet)$name]*5
E(tfNet)$width <- (abs(E(tfNet)$value.x))*20/max(abs(E(tfNet)$value.x))
E(tfNet)$color<-ifelse(E(tfNet)$value.x>0, "blue", "red")

pdf(file.path(outputDir,paste('Transition plot',analysisCode,'.pdf', sep="")), width=15, height=15)
plot.igraph(tfNet, edge.arrow.size=2, vertex.label.cex= 1.5, vertex.label.color= "black",main="")
legend("bottomleft", c("Gained features","Lost features"), lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"))
dev.off()

saveRDS(list(obsSsodm,dTFI_pVals_All),'dTFI.rdata')
# save.image(file=file.path(outputDir,paste("activeImage",analysisCode,".RData",sep="")))


# ############# Some methylation investigation 11/25/15
# # diffGeneMethPValues <- readRDS("./diffGeneMethPValues.rdata")
# meanTransition <- rowMeans(transMatrices[[1]])
# methDifferences <- diffGeneMethPValues[1,]-diffGeneMethPValues[2,] 
# matched_TF_Meth <- intersect(colnames(diffGeneMethPValues),names(obsSsodm))
# plot(methDifferences[matched_TF_Meth],meanTransition[matched_TF_Meth])
# summary(lm(methDifferences[matched_TF_Meth]~meanTransition[matched_TF_Meth]))

#library(igraph)
# # Visualize networks
# net.triple <- dataset$ppi
# net.triple <- melt(top.binary)
# g <- graph.data.frame(net.triple[net.triple[,3]==1,], directed=F)
# V(g)$color<-ifelse(V(g)$name%in%suspects, 'blue', 'red')
# l <- layout.fruchterman.reingold(g,niter=3,area=vcount(g)^2.3,repulserad=vcount(g)^5.8)
# plot(g,layout=l)
# 
# g <- graph.data.frame(dataset$motif[dataset$motif$V1%in%c(suspects,'ELK1'),], directed=F,types=rep(0:1,590))
# V(g)$color<-ifelse(V(g)$name%in%suspects, 'blue', 'red')
# l <- layout.bipartite (g, types = NULL, hgap = 10, vgap = 1, maxiter = 100) 
# plot(g,layout=l)
# 
# data <- dataset$motif[dataset$motif$V2%in%c(suspects,'ELK1','STAT3'),]
# inc <- spread(data,V2,V3,fill=0)
# rownames(inc) <- inc[,1]
# inc<-inc[,-1]
# g <- graph.incidence(top.binary)
# E(g)$color<-ifelse(E(g)$grade==9, "red", "grey")
# plot(g, layout=layout.bipartite,
#      vertex.color=c("green","cyan")[V(g)$type+1])
# 
# 
# # correlation matrix
# genes <- unique(dataset$motif[,1])
# heatmap.2(cor(t(dataset$exp[genes,])), dendrogram = "both", col=bluered, trace='none',cexCol=1.3,cexRow=1.3,margins=c(12,20))
# 
# 
# 
# # Strangely behaving TFs (because no motif data exists for these TFs)
# # Note 2/26/15  These should be filtered out of the analysis in this version of PandaR
# suspects <- c('HIF1A',
#               'HAND1',
#               'ARID3A',
#               'SOX5',
#               'PRRX2',
#               'NR1H2',
#               'POU5F1',
#               'VDR',
#               'NFE2L1',
#               'NKX3-1',
#               'MAFG',
#               'AHR',
#               'SOX10',
#               'DDIT3',
#               'TAL1',
#               'EWSR1',
#               'NFIC',
#               'TLX1',
#               'ARNT'
# )
# 
# 
# #################################
# ###  Genes approach
# ##################################
# 
# obs.s  <- apply(dataset$controlsNetwork-dataset$casesNetwork,2,function(x){x%*%x})
# null.s <- lapply(null.networks, function(nullnet){
#   apply(nullnet[[1]]-nullnet[[2]],2,function(x){x%*%x})
# })