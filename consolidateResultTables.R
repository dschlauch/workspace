library(ggplot2)
library(gridExtra)
library(ggExtra)
library(gtable)
library(VennDiagram)

##################################################################################################################
#  This script takes transition matrix results from multiple datasets and creates tables and plots comparing them.
#  Consolidate different results into single table
##################################################################################################################

# Specify the analysis folders and display names for the analyses
# analysisNames <- c("ECLIPSE_bere_bare_55557","COPDGene_bere_70856", "LGRC_bere_56432","LTCOPD_bere_bare_92540")
analysisNames <- c("ECLIPSE_combined_runs","COPDGene_combined_runs", "LGRC_combined_runs","LTCOPD_combined_runs")
displayNames <- c("ECLIPSE","COPDGene","LGRC","LTCOPD")
baseDir <- "~/gd/Harvard/Research/TM_outputs/JASPAR2014/BERE/"
outputDir <- './'
setwd(baseDir)
# Find the files for comparison
# Read them in.
# Add within-study rank order for magnitude and significance
filenames <- sapply(analysisNames, function(aname){
  dir(paste0(baseDir,aname), full.names = T, recursive = TRUE, all.files = TRUE, pattern="resultTable*") 
})
resultTables <- lapply(filenames, read.csv)
names(resultTables) <- analysisNames
resultTables <- lapply(resultTables, function(x){
    x$rankSig <- rank(x$dTFI.normalized.scores)
    x <- x[order(-x$Magnitude),]
    x$rankMag <- 1:nrow(x)
    x$limma <- exp(-(x$limma))
    x$negLogPValues <- -log(x$dTFI.normalized.scores)
    colnames(x) <- paste(analysisNames[parent.frame()$i], colnames(x), sep="_")
    x
})

# Merge results into single data.frame
merged.data.frame = Reduce(function(...) merge(..., by=1,all=T), resultTables)
merged.data.frame <- merged.data.frame[order(merged.data.frame[,2]),]

# Function for generation of a plot based on an index pair
makeComparisonPlot <- function(pair, plotTopNTFs=15, filterColIndices = c(8,18,28,38), metric="Magnitude", xlimits=c(0,.0325), ylimits=c(0,.0325)){
  # Include labels for any TFs that are in the top 15 of any list
  includedLabels <- apply(merged.data.frame[,filterColIndices[pair]],1,function(...) suppressWarnings(min(...,na.rm=T))) < plotTopNTFs
  merged.data.frame$labels <- as.character(merged.data.frame[,1])
  merged.data.frame$labels[!includedLabels] <- ""
  corValue <- cor(merged.data.frame[[paste(analysisNames[pair[1]], metric, sep="_")]], merged.data.frame[[paste(analysisNames[pair[2]], metric, sep="_")]], use="complete.obs", method="spearman")
  p.value <- cor.test(merged.data.frame[[paste(analysisNames[pair[1]], metric, sep="_")]], merged.data.frame[[paste(analysisNames[pair[2]], metric, sep="_")]], use="complete.obs", method="spearman")$p.value
  pValueText <- ifelse (p.value<1e-16,"p<1e-16",paste0("p=",as.character(p.value)))
  
  corText <- paste0("r[s]==",round(corValue,3),"~\n(",pValueText,")") 
  plot1 <- ggplot(merged.data.frame, aes_string(paste(analysisNames[pair[1]], metric, sep="_"),paste(analysisNames[pair[2]], metric, sep="_"), label="labels"))
  plot1 <- plot1 + geom_point(colour="blue",alpha=.5, size=4) + xlab(displayNames[pair[1]]) + ylab(displayNames[pair[2]]) + 
    geom_text(vjust=0) + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-1, label = corText, parse = TRUE, size = 8)+    
    scale_x_continuous(limits=xlimits, expand = c(0, 0)) +
    scale_y_continuous(limits=ylimits, expand = c(0, 0)) + 
    theme_classic()
  ggMarginal(plot1)
}

# Create the 4 comparison plots for ECLIPSE, LGRC, COPDGene, LTCOPD and combine them

generatePlots <- function(metric, filterColIndices){
  
  plot1 <- suppressWarnings(makeComparisonPlot(c(1,2),metric=metric, filterColIndices=filterColIndices, xlimits=c(0,.0325), ylimits=c(0,.025)))
  plot2 <- suppressWarnings(makeComparisonPlot(c(3,4),metric=metric, filterColIndices=filterColIndices, xlimits=c(0,.0325), ylimits=c(0,.025)))
  plot3 <- suppressWarnings(makeComparisonPlot(c(1,3),metric=metric, filterColIndices=filterColIndices))
  plot4 <- suppressWarnings(makeComparisonPlot(c(2,3),metric=metric, filterColIndices=filterColIndices))
  plot5 <- suppressWarnings(makeComparisonPlot(c(1,4),metric=metric, filterColIndices=filterColIndices))
  plot6 <- suppressWarnings(makeComparisonPlot(c(2,4),metric=metric, filterColIndices=filterColIndices))
  
  # suppressWarnings(grid.arrange(plot1, plot2, ncol=2, top="Comparison of Differential TF Involvement Across Studies of Same Tissue"))
  # suppressWarnings(grid.arrange(plot3, plot4, plot5, plot6, ncol=2, top="Comparison of Differential TF Involvement Across Studies of Different Tissues",
  #                               left = textGrob("Lung Tissue", rot = 90, vjust = 1), bottom = textGrob("Blood")))
  
  # Generate the pdfs for the above plots
  pdf(paste0(outputDir, metric, '_comparison_same_tissue.pdf'), width=16, height=9)
  suppressWarnings(grid.arrange(plot1, plot2, ncol=2, top=textGrob("Comparison of Differential TF Involvement Across Studies of Same Tissue", gp=gpar(fontsize=30))))
  dev.off()
  png(paste0(outputDir, metric, '_comparison_same_tissue.png'), width=1600, height=900)
  suppressWarnings(grid.arrange(plot1, plot2, ncol=2, top=textGrob("Comparison of Differential TF Involvement Across Studies of Same Tissue", gp=gpar(fontsize=30))))
  dev.off()
  pdf(paste0(outputDir, metric, '_comparison_diff_tissue.pdf'), width=16,height=12)
  suppressWarnings(grid.arrange(plot3, plot4, plot5, plot6, ncol=2, top=textGrob("Comparison of Differential TF Involvement Across Studies of Different Tissues", gp=gpar(fontsize=30)),
                                left = textGrob("Lung Tissue", rot = 90, vjust = 1, gp=gpar(fontsize=20)), bottom = textGrob("Blood", gp=gpar(fontsize=20))))
  dev.off()
  png(paste0(outputDir, metric, '_comparison_diff_tissue.png'), width=1600,height=900)
  suppressWarnings(grid.arrange(plot3, plot4, plot5, plot6, ncol=2, top=textGrob("Comparison of Differential TF Involvement Across Studies of Different Tissues", gp=gpar(fontsize=30)),
                                left = textGrob("Lung Tissue", rot = 90, vjust = 1, gp=gpar(fontsize=20)), bottom = textGrob("Blood", gp=gpar(fontsize=20))))
  dev.off()
  
  topTFlist <- lapply(filterColIndices, function(i){
    merged.data.frame[merged.data.frame[,i] %in% 1:20,1]
  })
  names(topTFlist) <- displayNames
  venn.diagram(x = topTFlist,
               filename = "Venn.tiff",height = 2000, width = 3000,
               col = "transparent", fill = c("green","yellow","darkorchid1","lightblue"),
               alpha = 0.50, label.col = rep("white",15), cex = 1.5, fontfamily = "serif", fontface = "bold",
               cat.col = rep("black",4), cat.cex = 1.5,
               cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", rotation.degree = 0,
               margin = 0.2, main="Top 20 differentially involved TFs in each COPD Study",main.pos= c(0.5, .9),main.cex = 1.4)
}
metric <- "negLogPValues"
# metric <- "Magnitude"
generatePlots("negLogPValues",filterColIndices = c(8,18,28,38))
generatePlots("Magnitude", filterColIndices = c(9,19,29,39))

# Create table
keepColnames <- c(t(outer(analysisNames,c("Magnitude","rankMag","dTFI.FDR","limma"), paste, sep="_")))
displayColnames <- rep(c("dTFI","rank","FDR","LIMMA"), 4)
publicationTable <- cbind(merged.data.frame[,1], round(merged.data.frame[,keepColnames],4))
colnames(publicationTable) <- c('TF',displayColnames)
options(scipen=999)
publicationTable <- publicationTable[order(-publicationTable[,2]),]
publicationTable <- publicationTable[which(apply(publicationTable[,c(3,7,11)],1,min)<11),]
publicationTable[publicationTable==0] <- "<.0001"
# Table is now a character table (not numeric)
publicationTable <- data.frame(gsub("0.", ".", as.matrix(publicationTable),fixed=T),check.names=FALSE)
dim(publicationTable)
write.table(publicationTable, file="study_comparison_table.txt", row.names=F, sep="\t", quote=F)
options(scipen=0)