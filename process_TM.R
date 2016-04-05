# Do the sum of sq ODM plot versus null
pdf(file.path(outputDir,paste0('SSODMplot_unscaled',analysisCode,'.pdf')), width=8)
sPlot <- ssodm.plot(transMatrices[[1]], transMatrices[-1],plot.title=paste("SSODM observed and null, ",casesString," vs ",controlsString,' : ', networkInferenceName, ' : ', analysisName, sep=""))
print(sPlot)
dev.off()
pdf(file.path(outputDir,paste0('SSODMplot_scaled',analysisCode,'.pdf')), width=24)
sPlot <- ssodm.plot(transMatrices[[1]], transMatrices[-1], rescale=T,plot.title="")#, plot.title=paste("SSODM observed and null, ",casesString," vs ",controlsString,' : ', networkInferenceName, ' : ', analysisName, sep=""))
print(sPlot)
dev.off()


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
labels[rank(-negLogZPValues)>20 & rank(-limmanegLogPValues)>20]<-""

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

dTFI_LIMMA_gg <- ggplot(data=plotDF, aes(x=limmanegLogPValues, y=negLogZPValues)) + geom_point(aes(col=logfoldchangeTF), size=7, alpha=.8) +
  geom_text_repel(data=plotDF[labels!="",], aes(limmanegLogPValues, negLogZPValues, label=labels), size = 10,) + 
  ylab("Differential TF Involvement, -log(p-value)") + xlab("Differential Expression,  LIMMA -log(p-value)") + 
  ggtitle(expression(atop("Differential Involvement vs Differential Expression (ECLIPSE)", atop(italic("Smoker Controls to COPD Patients"), ""))))+
  scale_colour_gradient2(limits=c(-max(abs(logfoldchangeTF))/4,max(abs(logfoldchangeTF))/4), name="log(fold-change)", low = "blue", high = "yellow", mid="black") +
  theme_classic() + 
  theme(plot.title = element_text(size=25,hjust=.5), axis.title=element_text(size=22), legend.text=element_text(size=20), legend.title=element_text(size=20))

pdf(file.path(outputDir,paste('dTFI vs LIMMA',analysisCode,'.pdf', sep="")), width=12, height=12)
print(dTFI_LIMMA_gg)
dev.off()

png(file.path(outputDir,paste('dTFI vs LIMMA',analysisCode,'.png', sep="")), width=900, height=800)
print(dTFI_LIMMA_gg)
dev.off()

# Generate heatmap
x <- transMatrices[[1]]
mdf <- melt(x)
mdf[,2] <- factor(mdf[,2],levels=sort(levels(mdf[,2]), decreasing=T))
# Heatmap
p1 <- ggplot(mdf, aes(x=Var1, y=Var2)) +
  geom_tile(aes(fill=value)) + scale_fill_gradient2(name = "dTFI") + theme_bw() + 
  theme(axis.text.y = element_text(size=3), axis.text.x = element_text(size=3,angle = 90, hjust = 1)) + 
  xlab("Transcription Factors") + ylab("Transcription Factors") + ggtitle("ECLIPSE Transition Matrix")

pdf(file.path(outputDir,paste('TM_heatmap.pdf', sep="")), width=9, height=8)
print(p1)
dev.off()

# For the paper with no labels
p1 <- ggplot(mdf, aes(x=Var1, y=Var2)) +
  geom_tile(aes(fill=value)) + 
  xlab("Transcription Factors") + ylab("Transcription Factors") + scale_fill_gradient2(name = "dTFI") + 
  theme_bw() + theme(axis.ticks = element_blank(), axis.title=element_text(size=25), legend.title=element_text(size=20), axis.text.y = element_blank(), axis.text.x = element_blank(), plot.title=element_text(family="Times", face="bold", size=40)) + 
  ggtitle(expression(atop("ECLIPSE Transition Matrix", atop(italic("Smoker Controls to COPD Patients"), ""))))

png(file.path(outputDir,paste('TM_heatmap.png', sep="")), width=800, height=800)
print(p1)
dev.off()
pdf(file.path(outputDir,paste('TM_heatmap.pdf', sep="")), width=8, height=8)
print(p1)
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