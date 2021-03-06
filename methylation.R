# Open session on interactive node
# qrsh -l vf=10G -l R=1
library("FDb.InfiniumMethylation.hg19")
load('/proj/regeps/regep00/studies/LTCOPD/data/methylation/ECLIPSE/methyECLIPSE_clean_1430052791_V9.RData')
mapping <- read.table('/proj/regeps/regep00/studies/ECLIPSE/analyses/rekrg/PANDA_GENDER_METHYLATION/CompileData/CELs_with_Meth.txt')

################ Kimbie's data
# oldBetas <- read.table('~/methylation/CombinedMethylationBetaValsAll_011315.txt',header=T,row.names=1)
# kimbiepheno <- readlines('/udd/rekrg/ecldata/genetics_phenotype20090810.txt',header=T,sep="\t",fill=T, quote="")
# 
# sampleIDs <- substring(colnames(oldBetas),2)
# copdSamples <- kimbiepheno[match(sampleIDs,kimbiepheno[,1]),4]=="COPD Subjects"
# smcSamples <- kimbiepheno[match(sampleIDs,kimbiepheno[,1]),4]=="Smoker Controls"
# copdSamples[is.na(copdSamples)] <- F
# smcSamples[is.na(smcSamples)] <- F
# 
# diffMethPValues <- apply(oldBetas,1,function(x){
#   t.test(x[copdSamples],x[smcSamples],alternative="two.sided",var.equal=T)$p.value
# })
# sort(diffMethPValues)[1:10]
# sum(diffMethPValues<.05)
# diffMethFDR <- p.adjust(diffMethPValues, method = 'fdr')

# Mapping CpGs to Genes
hm450 <- get450k()
probenames <- rownames(betas)
probes <- hm450[probenames]

genes <- getNearestGene(probes)$nearestGeneSymbol

genesMatrix <- betas
rownames(genesMatrix) <- genes

#Takes ~10min
geneBySampleMeth <- aggregate(genesMatrix, by=list(genes), mean)
rownames(geneBySampleMeth) <- geneBySampleMeth[,1]
geneBySampleMeth <- geneBySampleMeth[,-1]

# cd ~/methylation
phenoData <- read.table("~/methylation/ECLIPSE_phenotypes.txt",header=T,fill = TRUE, sep="\t",row.names=1)
phenoData$GroupAccession <- rownames(phenoData)
rownames(phenoData) <- substring(sub('.*_WB_','',phenoData$Title),1,6)
methylCelFiles <- as.character(mapping[match(rownames(pDat),mapping[,1]),2])
phenoData <- phenoData[methylCelFiles,]

cases <- !is.na(phenoData$Subject.type)&phenoData$Subject.type=="COPD"
controls <- !is.na(phenoData$Subject.type)&phenoData$Subject.type=="Smoker Control"

# Gene level analysis
diffMethPValues <- apply(geneBySampleMeth,1,function(x){
  tt <- t.test(x[cases],x[controls],alternative="two.sided") 
  c(tt$estimate, tt$p.value)
})
saveRDS(diffMethPValues,"~/methylation/diffGeneMethPValues.rdata")

# CpG level analysis
diffMethPValues <- apply(betas,1,function(x){
  tt <- t.test(x[cases],x[controls],alternative="two.sided") 
  c(tt$estimate, tt$p.value)
})
saveRDS(diffMethPValues,"~/methylation/diffCpGMethPValues.rdata")

# To be run locally for plotting, etc...
# scp redsc@rock.bwh.harvard.edu:~/methylation/diffGeneMethPValues.rdata ~/gd/Harvard/Research/R_workspace/
# scp redsc@rock.bwh.harvard.edu:~/methylation/diffCpGMethPValues.rdata ~/gd/Harvard/Research/R_workspace/

negLogQQ <- function(x){
  title <- deparse(substitute(x))
  theoreticalVals <- -log((1:length(x))/length(x))
  obsVals <- -log(sort(x))
  plot(theoreticalVals, obsVals, main=title)
  abline(0,1)  
}
diffCpGMethPValues <- readRDS("./diffCpGMethPValues.rdata")
negLogQQ(diffCpGMethPValues[3,])
diffCpGMethPValues[,order(-abs(diffCpGMethPValues[1,]-diffCpGMethPValues[2,]))[1:10]]
aggregate(diffCpGMethPValues,)

hm450 <- get450k()
probenames <- colnames(diffCpGMethPValues)
probes <- hm450[probenames]

genes <- getNearestGene(probes)$nearestGeneSymbol

# Get dTFI results
library(ggplot2)
library(gridExtra)
library(ggExtra)

load("~/gd/Harvard/Research/TM_outputs/ECLIPSE_bere_bare_55557/activeImage55557.RData")
TFValues1 <- sort(obsSsodm)
cpgs_for_TFs <- diffCpGMethPValues[,genes %in% names(TFValues1)]
tfs1 <- factor(genes[genes %in% names(TFValues1)], levels=names(TFValues1))
negLogPval <- -log(cpgs_for_TFs[3,])
plot1 <- qplot(tfs1, negLogPval) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(aes(colour = negLogPval)) + ggtitle("ordered by effect size")

TFValues2 <- sort(dTFI_pVals, decreasing=T)
cpgs_for_TFs <- diffCpGMethPValues[,genes %in% names(TFValues2)]
tfs2 <- factor(genes[genes %in% names(TFValues2)], levels=names(TFValues2))
negLogPval <- -log(cpgs_for_TFs[3,])
plot2 <- qplot(tfs2, negLogPval) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_point(aes(colour = negLogPval)) + ggtitle("ordered by significance")

png('diffMethylCpG', width=1800)
grid.arrange(plot1, plot2, nrow=2, top="differentially methylated CpG sites near TFs")
dev.off()

diffGeneMethPValues <- readRDS("./diffGeneMethPValues.rdata")
negLogQQ(diffGeneMethPValues[3,matched_TF_Meth])
diffGeneMethPValues[,order(-abs(diffGeneMethPValues[1,]-diffGeneMethPValues[2,]))[1:10]]
