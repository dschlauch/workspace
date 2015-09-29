yeastDir <- "~/gd/Harvard/Research/data/YeastData/"
motifs <- read.table(file.path(yeastDir,"YeastData_Motif.txt"), stringsAsFactors=F)
transFactors <- as.character(unique(motifs[,1]))

lapply(transFactors, function(tf){
    write.table(cbind(motifs[motifs[,1]==tf,2:3], 100, 1000, "forward"), file=file.path("~/gd/Harvard/Research/SEREND/YeastCC/MotifScores", paste(tf,".txt",sep="")), quote=F, row.names=F, col.names=F, sep="\t")
})
write.table(transFactors, file=file.path("~/gd/Harvard/Research/SEREND/YeastCC/", "TFlist.txt"), quote=F, row.names=F, col.names=F, sep="\n")

directEvidence <- motifs
directEvidence[,3] <- '+-'
write.table(directEvidence, file=file.path("~/gd/Harvard/Research/SEREND/YeastCC/", "directEvidence.txt"), quote=F, row.names=F, col.names=F, sep="\t")
