library(ggplot2)
setwd('~')
outputDir <- "NI_only_0001"
eclipseCases <- readRDS(file.path(outputDir,'ECLIPSE_cases_bere_network.rds'))
eclipseControls <- readRDS(file.path(outputDir,'ECLIPSE_controls_bere_network.rds'))
COPDGeneCases <- readRDS(file.path(outputDir,'COPDGene_cases_bere_network.rds'))
COPDGeneControls <- readRDS(file.path(outputDir,'COPDGene_controls_bere_network.rds'))


matchedTFs <- intersect(rownames(COPDGeneCases),rownames(eclipseCases))
matchedGenes <- intersect(colnames(COPDGeneCases),colnames(eclipseCases))

COPDGeneCases <- COPDGeneCases[matchedTFs,matchedGenes]
COPDGeneControls <- COPDGeneControls[matchedTFs,matchedGenes]
eclipseCases <- eclipseCases[matchedTFs,matchedGenes]
eclipseControls <- eclipseControls[matchedTFs,matchedGenes]

df <- data.frame(cases=c(COPDGeneCases),controls=c(COPDGeneControls))
cor(df)
png(file.path(outputDir,'COPDGene_edgeweight_comparison.png'))
ggplot(df, aes(x=controls, y=cases)) +
  geom_point(size=.1, alpha=.1) + xlab("Controls") + ylab("Cases") + ggtitle("Edgeweights in COPDGene")+ theme_classic() 
dev.off()

df <- data.frame(cases=c(eclipseCases),controls=c(eclipseControls))
cor(df)
png(file.path(outputDir,'ECLIPSE_edgeweight_comparison.png'))
ggplot(df, aes(x=controls, y=cases)) +
  geom_point(size=.1, alpha=.1) + xlab("Controls") + ylab("Cases") + ggtitle("Edgeweights in ECLIPSE")+ theme_classic() 
dev.off()

df <- data.frame(eclipse=c(eclipseCases),copdgene=c(COPDGeneCases))
cor(df)
png(file.path(outputDir,'COPDGene_vs_ECLIPSE_edgeweight_cases_comparison.png'))
ggplot(df, aes(x=eclipse, y=copdgene)) +
  geom_point(size=.1, alpha=.1) + xlab("ECLIPSE") + ylab("COPDGene") + ggtitle("Edgeweights in Cases (ECLIPSE,COPDGene)") + theme_classic() 
dev.off()

df <- data.frame(eclipse=c(eclipseControls),copdgene=c(COPDGeneControls))
cor(df)
png(file.path(outputDir,'COPDGene_vs_ECLIPSE_edgeweight_controls_comparison.png'))
ggplot(df, aes(x=eclipse, y=copdgene)) +
  geom_point(size=.1, alpha=.1) + xlab("ECLIPSE") + ylab("COPDGene") + ggtitle("Edgeweights in Controls (ECLIPSE,COPDGene)") + theme_classic() 
dev.off()
