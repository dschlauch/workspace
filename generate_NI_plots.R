library(ggplot2)
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
  
  
  df <- data.frame(dataA=c(datasetACases)-c(datasetAControls), dataB=c(datasetBCases)-c(datasetBControls))
  spearCor <- cor(df[,1],df[,2],method="spearman")
  print(spearCor)
  imageType(file.path(outputDir,paste(datasetA, datasetB, niMethod, 'edgeweight_difference_comparison.png', sep="_")))
  plot <- ggplot(df, aes(x=dataA, y=dataB)) +
    geom_point(size=.1, alpha=.1) + xlab(datasetA) + ylab(datasetB) + ggtitle(paste0("Edgeweights in Controls (", datasetA, datasetB, ")")) + theme_classic()+
    annotate("text", x = 0.0, y = -Inf, hjust=0, label = paste0("r=",round(spearCor,4)), parse = TRUE, size = 8)
  print(plot)
  dev.off()
  
}

generateNIDifferencePlot("ECLIPSE","COPDGene","bere",png)
generateNIDifferencePlot("ECLIPSE","LGRC","bere",png)
generateNIDifferencePlot("ECLIPSE","LTCOPD","bere",png)
generateNIDifferencePlot("ECLIPSE","COPDGene","WGCNA",png)
generateNIDifferencePlot("ECLIPSE","LGRC","WGCNA",png)
generateNIDifferencePlot("ECLIPSE","LTCOPD","WGCNA",png)
generateNIDifferencePlot("ECLIPSE","COPDGene","CLR",png)
generateNIDifferencePlot("ECLIPSE","LGRC","CLR",png)
generateNIDifferencePlot("ECLIPSE","LTCOPD","CLR",png)

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
