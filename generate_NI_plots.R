library(ggplot2)
library(data.table)
setwd('~')
outputDir <- "NI_only_0001"

generateNIDifferencePlot <- function(datasetA, datasetB, niMethod, imageType=png){
  datasetACases <- readRDS(file.path(outputDir, paste0(datasetA, '_cases_', niMethod, '_network.rds')))
  datasetAControls <- readRDS(file.path(outputDir, paste0(datasetA, '_controls_', niMethod, '_network.rds')))
  datasetBCases <- readRDS(file.path(outputDir, paste0(datasetB, '_cases_', niMethod, '_network.rds')))
  datasetBControls <- readRDS(file.path(outputDir, paste0(datasetB, '_controls_', niMethod, '_network.rds')))
  
  matchedTFs <- intersect(rownames(datasetACases),rownames(datasetBCases))
  matchedGenes <- intersect(colnames(datasetACases),colnames(datasetBCases))
  
  datasetACases    <- datasetACases[matchedTFs,matchedGenes]
  datasetAControls <- datasetAControls[matchedTFs,matchedGenes]
  datasetBCases    <- datasetBCases[matchedTFs,matchedGenes]
  datasetBControls <- datasetBControls[matchedTFs,matchedGenes]
  
  
  df <- data.frame(dataA=c(datasetACases)-c(datasetAControls), dataB=c(datasetBCases)-c(datasetBControls), studyA=datasetA, studyB=datasetB, method=niMethod)
  df
}

studies <- c("ECLIPSE", "COPDGENE", "LGRC", "LTCOPD")
niMethods <- c("WGCNA","CLR","MONSTER","ARACNE")
resultsList <- lapply(niMethods, function(meth){
  data.table(do.call(rbind, apply(combn(studies,2), 2, function(x){
    print(x)  
    generateNIDifferencePlot(x[1],x[2],meth)
  })))
})

studyCorrelations <- data.table(do.call(rbind,lapply(resultsList, function(resDT){
  t(apply(combn(studies,2), 2, function(x){
    includerows <- resDT$studyA==x[1] & resDT$studyB==x[2] 
    dt <- subset(resDT, includerows)
    c(x[1],x[2], paste0("r=",round(cor(dt[,dataA],dt[,dataB]),4)))
  }))
})))


colnames(studyCorrelations) <- c("studyA","studyB","corText")
studyCorrelations[['method']] <- rep(niMethods,each=6)
# remove points for plotting purposes
resultsList <- lapply(resultsList, function(x){x[c(rep(F,9),T)]})
studyCorrelations[[1]] <- factor(studyCorrelations[[1]], levels=studies)
studyCorrelations[[2]] <- factor(studyCorrelations[[2]], levels=studies)
lapply(seq_along(niMethods), function(i){
  resultsList[[i]][['studyA']] <- factor(resultsList[[i]][['studyA']], levels=studies)
  resultsList[[i]][['studyB']] <- factor(resultsList[[i]][['studyB']], levels=studies)
  tiff(file.path(outputDir,"plots",paste("all_studies", niMethods[i], 'edgeweight_difference_comparison.tiff', sep="_")),height=2000,width=2000, units = "px", res = 400)
  plot <- ggplot(resultsList[[i]], aes(x=dataA, y=dataB)) + facet_grid(studyA ~ studyB) +
    geom_point(size=.1, alpha=.6, col="blue") + 
    ggtitle(paste(niMethods[i], "edge weight differences between pairs of studies")) + 
    theme_bw() + theme(plot.title = element_text(size=11, face="bold")) + 
    xlab("") + ylab("") +
    geom_text(data=studyCorrelations[method==niMethods[i]], aes(x=min(resultsList[[i]][,dataA]),y=max(min(resultsList[[i]][,dataB])),label=corText), size=3,hjust=0,vjust=0) 
  print(plot)
  dev.off()
})