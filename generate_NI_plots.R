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
#   spearCor <- cor(df[,1],df[,2],method="spearman")
#   corText <- paste0("r[s]==",round(spearCor,4))
#   print(spearCor)
#   imageType(file.path(outputDir,"plots",paste(datasetA, datasetB, niMethod, 'edgeweight_difference_comparison.png', sep="_")))
#   plot <- ggplot(df, aes(x=dataA, y=dataB)) +
#     geom_point(size=.1, alpha=.1) + xlab(datasetA) + ylab(datasetB) + ggtitle(paste0(niMethod, " edge weight differences (", datasetA," vs ", datasetB, ")")) + theme_classic()+
#     annotate("text", x=Inf, y=-Inf, hjust=1, , vjust=-1, label=corText, parse=TRUE, size=8)
#   print(plot)
#   dev.off()
  df
}

studies <- c("ECLIPSE", "COPDGENE", "LGRC", "LTCOPD")
niMethods <- c("WGCNA","CLR","bere","ARACNE")
resultsList <- apply(combn(studies,2), 2, function(x){
  print(x)  
  do.call(rbind, lapply(niMethods, function(meth) generateNIDifferencePlot(x[1],x[2],meth) ))
})
combinedDT <- data.table(do.call(rbind, resultsList))
saveRDS(combinedDT,file.path(outputDir,"allStudiesNIResults.rds"))
# combinedDT <- readRDS(file.path(outputDir,"allStudiesNIResults.rds"))

studyCorrelations <- data.frame(t(apply(combn(studies,2), 2, function(x){
    sapply(niMethods, function(niMethod){
        includerows <- combinedDT$studyA==x[1] & combinedDT$studyB==x[2] & combinedDT$method==niMethod
        dt <- subset(combinedDT, includerows)
        c(x[1],x[2], paste0("r=",round(cor(dt[,dataA],dt[,dataB]),4)),niMethod)
    })
})))
colnames(studyCorrelations) <- c("studyA","studyB","corText")
# remove points for plotting purposes
combinedDTFiltered <- combinedDT[c(rep(F,9),T)]

lapply(niMethods, function(niMethod){
  png(file.path(outputDir,"plots",paste("all_studies", niMethod, 'edgeweight_difference_comparison.png', sep="_")),height=1200,width=1200)
  plot <- ggplot(subset(combinedDTFiltered,method==niMethod), aes(x=dataA, y=dataB)) + facet_grid(studyA ~ studyB) +
    geom_point(size=.1, alpha=1, col="blue") + 
    ggtitle(paste0(niMethod, " edge weight differences between pairs of studies")) + 
    theme_bw() + theme(plot.title = element_text(size=25, face="bold")) + 
    geom_text(data = studyCorrelations, aes(x=min(subset(combinedDTFiltered,method==niMethod)[,dataA]),y=max(min(subset(combinedDTFiltered,method==niMethod)[,dataB])),label =corText), size=10,hjust=0,vjust=0) 
  print(plot)
  dev.off()
})

# generateNIDifferencePlot("ECLIPSE","COPDGENE","bere",png)
# generateNIDifferencePlot("ECLIPSE","LGRC","bere",png)
# generateNIDifferencePlot("ECLIPSE","LTCOPD","bere",png)
# generateNIDifferencePlot("ECLIPSE","COPDGENE","WGCNA",png)
# generateNIDifferencePlot("ECLIPSE","LGRC","WGCNA",png)
# generateNIDifferencePlot("ECLIPSE","LTCOPD","WGCNA",png)
# generateNIDifferencePlot("ECLIPSE","COPDGENE","CLR",png)
# generateNIDifferencePlot("ECLIPSE","LGRC","CLR",png)
# generateNIDifferencePlot("ECLIPSE","LTCOPD","CLR",png)
# generateNIDifferencePlot("ECLIPSE","COPDGENE","ARACNE",png)
# generateNIDifferencePlot("ECLIPSE","LGRC","ARACNE",png)
# generateNIDifferencePlot("ECLIPSE","LTCOPD","ARACNE",png)


# 
# eclipseBERECases <- readRDS(file.path(outputDir,'ECLIPSE_cases_bere_network.rds'))
# eclipseBEREControls <- readRDS(file.path(outputDir,'ECLIPSE_controls_bere_network.rds'))
# COPDGeneBERECases <- readRDS(file.path(outputDir,'COPDGene_cases_bere_network.rds'))
# COPDGeneBEREControls <- readRDS(file.path(outputDir,'COPDGene_controls_bere_network.rds'))
# 
# 
# matchedTFs <- intersect(rownames(COPDGeneBERECases),rownames(eclipseBERECases))
# matchedGenes <- intersect(colnames(COPDGeneBERECases),colnames(eclipseBERECases))
# 
# COPDGeneBERECases <- COPDGeneBERECases[matchedTFs,matchedGenes]
# COPDGeneBEREControls <- COPDGeneBEREControls[matchedTFs,matchedGenes]
# eclipseBERECases <- eclipseBERECases[matchedTFs,matchedGenes]
# eclipseBEREControls <- eclipseBEREControls[matchedTFs,matchedGenes]
# 
# df <- data.frame(cases=c(COPDGeneBERECases),controls=c(COPDGeneBEREControls))
# cor(df)
# png(file.path(outputDir,'COPDGene_edgeweight_comparison.png'))
# ggplot(df, aes(x=controls, y=cases)) +
#   geom_point(size=.1, alpha=.1) + xlab("Controls") + ylab("Cases") + ggtitle("Edgeweights in COPDGene")+ theme_classic() 
# dev.off()
# 
# df <- data.frame(cases=c(eclipseBERECases),controls=c(eclipseBEREControls))
# cor(df)
# png(file.path(outputDir,'ECLIPSE_edgeweight_comparison.png'))
# ggplot(df, aes(x=controls, y=cases)) +
#   geom_point(size=.1, alpha=.1) + xlab("Controls") + ylab("Cases") + ggtitle("Edgeweights in ECLIPSE")+ theme_classic() 
# dev.off()
# 
# df <- data.frame(eclipse=c(eclipseBERECases),copdgene=c(COPDGeneBERECases))
# cor(df)
# png(file.path(outputDir,'COPDGene_vs_ECLIPSE_edgeweight_cases_comparison.png'))
# ggplot(df, aes(x=eclipse, y=copdgene)) +
#   geom_point(size=.1, alpha=.1) + xlab("ECLIPSE") + ylab("COPDGene") + ggtitle("Edgeweights in Cases (ECLIPSE,COPDGene)") + theme_classic() 
# dev.off()
# 
# df <- data.frame(eclipse=c(eclipseBEREControls),copdgene=c(COPDGeneBEREControls))
# cor(df)
# png(file.path(outputDir,'COPDGene_vs_ECLIPSE_edgeweight_controls_comparison.png'))
# ggplot(df, aes(x=eclipse, y=copdgene)) +
#   geom_point(size=.1, alpha=.1) + xlab("ECLIPSE") + ylab("COPDGene") + ggtitle("Edgeweights in Controls (ECLIPSE,COPDGene)") + theme_classic() 
# dev.off()
# 
# # Edgeweight Difference
# 
# df <- data.frame(eclipse=c(eclipseBERECases)-c(eclipseBEREControls), copdgene=c(COPDGeneBERECases)-c(COPDGeneBEREControls))
# cor(df)
# png(file.path(outputDir,'COPDGene_ECLIPSE_edgeweight_difference_comparison.png'))
# ggplot(df, aes(x=eclipse, y=copdgene)) +
#   geom_point(size=.1, alpha=.1) + xlab("ECLIPSE") + ylab("COPDGene") + ggtitle("Edgeweights in Controls (ECLIPSE,COPDGene)") + theme_classic() 
# dev.off()
# 
# 
# 
# # Only TF-TF
# TFsWithGenes <- intersect(matchedTFs,matchedGenes)
# 
# COPDGeneBERECasesTFs <- COPDGeneBERECases[TFsWithGenes,TFsWithGenes]
# COPDGeneBEREControlsTFs <- COPDGeneBEREControls[TFsWithGenes,TFsWithGenes]
# eclipseBERECasesTFs <- eclipseBERECases[TFsWithGenes,TFsWithGenes]
# eclipseBEREControlsTFs <- eclipseBEREControls[TFsWithGenes,TFsWithGenes]
# 
# 
# df <- data.frame(eclipse=c(COPDGeneBERECasesTFs),copdgene=c(COPDGeneBEREControlsTFs))
# cor(df)
# png(file.path(outputDir,'COPDGene_TFTF_edgeweight_controls_comparison.png'))
# ggplot(df, aes(x=eclipse, y=copdgene)) +
#   geom_point(size=.1, alpha=.1) + xlab("ECLIPSE") + ylab("COPDGene") + ggtitle("Edgeweights in Controls (ECLIPSE,COPDGene)") + theme_classic() 
# dev.off()
# 
# df <- data.frame(eclipse=c(eclipseBEREControlsTFs),copdgene=c(COPDGeneBEREControlsTFs))
# cor(df)
# png(file.path(outputDir,'COPDGene_vs_ECLIPSE_TFTF_edgeweight_controls_comparison.png'))
# ggplot(df, aes(x=eclipse, y=copdgene)) +
#   geom_point(size=.1, alpha=.1) + xlab("ECLIPSE") + ylab("COPDGene") + ggtitle("Edgeweights in Controls (ECLIPSE,COPDGene)") + theme_classic() 
# dev.off()
# 
# df <- data.frame(eclipse=c(eclipseBERECasesTFs),copdgene=c(COPDGeneBERECasesTFs))
# cor(df)
# png(file.path(outputDir,'COPDGene_vs_ECLIPSE_TFTF_edgeweight_controls_comparison.png'))
# ggplot(df, aes(x=eclipse, y=copdgene)) +
#   geom_point(size=.1, alpha=.1) + xlab("ECLIPSE") + ylab("COPDGene") + ggtitle("Edgeweights in Controls (ECLIPSE,COPDGene)") + theme_classic() 
# dev.off()
