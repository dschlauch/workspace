library(ggplot2)
library(gridExtra)
library(ggExtra)

##################################################################################################################
#  This script takes transition matrix results from multiple datasets and creates tables and plots comparing them.
#  Consolidate different results into single table
##################################################################################################################

# Specify the analysis folders and display names for the analyses
analysisNames <- c("ECLIPSE_bere_bare_55557","LGRC_bere_38678","COPDGene_bere_70856")
displayNames <- c("ECLIPSE","LGRC","COPDGene")

# Find the files for comparison
# Read them in.
# Add within-study rank order for magnitude and significance
filenames <- dir(paste0("~/gd/Harvard/Research/TM_outputs/",analysisNames), full.names = T, recursive = TRUE, all.files = TRUE, pattern="resultTable*") 
resultTables <- lapply(filenames, read.csv)
names(resultTables) <- analysisNames
resultTables <- lapply(resultTables, function(x){
    x$rankSig <- 1:nrow(x)
    x <- x[order(-x$Magnitude),]
    x$rankMag <- 1:nrow(x)
    x$limma <- exp(-(x$limma))
    colnames(x) <- paste(analysisNames[parent.frame()$i], colnames(x), sep="_")
    x
})

# Merge results into single data.frame
merged.data.frame = Reduce(function(...) merge(..., by=1,all=T), resultTables)
merged.data.frame <- merged.data.frame[order(merged.data.frame[,2]),]

# Function for generation of a plot based on an index pair
makeComparisonPlot <- function(pair, plotTopNTFs=15, filterColIndices = c(8,16,24)){
  # Include labels for any TFs that are in the top 15 of any list
  includedLabels <- apply(merged.data.frame[,filterColIndices[pair]],1,function(...) suppressWarnings(min(...,na.rm=T))) < plotTopNTFs
  merged.data.frame$labels <- as.character(merged.data.frame[,1])
  merged.data.frame$labels[!includedLabels] <- ""
  corText <- paste0("R^{2}==",round(cor(merged.data.frame[[paste(analysisNames[pair[1]], "Magnitude", sep="_")]], merged.data.frame[[paste(analysisNames[pair[2]], "Magnitude", sep="_")]], use="complete.obs"),4)) 
  plot1 <- ggplot(merged.data.frame, aes_string(paste(analysisNames[pair[1]], "Magnitude", sep="_"),paste(analysisNames[pair[2]], "Magnitude", sep="_"), label="labels"))
  plot1 <- plot1 + geom_point(colour="blue",alpha=.5, size=4) + xlab(displayNames[pair[1]]) + ylab(displayNames[pair[2]]) + geom_text(vjust=0)  + expand_limits(x=c(0,.025)) + annotate("text", x = .02, y = 0, label = corText, parse = TRUE)
  ggMarginal(plot1)
}

# Create the 3 comparison plots for ECLIPSE, LGRC, COPDGene and combine them
plot1 <- makeComparisonPlot(c(1,2))
plot2 <- makeComparisonPlot(c(2,3))
plot3 <- makeComparisonPlot(c(3,1))
suppressWarnings(grid.arrange(plot1, plot2, plot3, ncol=3, top="Comparison of Differential TF Involvement Across Studies"))

# Generate the png for the above plots
png('eclipse_copdgene_lgrc_comparison.png', width=1800)
suppressWarnings(grid.arrange(plot1, plot2, plot3, ncol=3, top="Comparison of Differential TF Involvement Across Studies"))
dev.off()

# Create table
keepColnames <- c(t(outer(analysisNames,c("Magnitude","rankMag","dTFI.FDR","limma..logp"), paste, sep="_")))
displayColnames <- rep(c("Magnitude","rank","FDR","LIMMA"), 3)
publicationTable <- cbind(merged.data.frame[,1], round(merged.data.frame[,keepColnames],4))
colnames(publicationTable) <- c('TF',displayColnames)
publicationTable <- publicationTable[order(-publicationTable[,2]),]
publicationTable[publicationTable==0] <- "<.0001"
head(publicationTable)
