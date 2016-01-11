library(data.table)
ECLIPSE_Blood_Motif <- data.table(read.delim("~/gd/Harvard/Research/data/Eclipse/ECLIPSE_Blood_Motif.txt", header=FALSE, stringsAsFactors=FALSE))
ECLIPSE_Blood_Motif
ECLIPSE_Blood_Motif[V1=="GABPA",]
print(ECLIPSE_Blood_Motif[V2=="IREB2",])


motifs_new <- read.delim("~/gd/Harvard/Research/data/GTEx/KG_cisbp_652(with1).txt", header=FALSE, stringsAsFactors=F)
TFmappings <- read.table("~/gd/Harvard/Research/data/GTEx/cisbpall_motinf.txt", header=T)
TFmappings[,1] <- substring(TFmappings[,1],0,5) 

library(org.Hs.eg.db)
symbols <- mapIds(org.Hs.eg.db, keys=motifs_new[,2],column="SYMBOL", keytype="ENSEMBL", multiVals="first")
motifs_new[,2] <- symbols
rownames(dataset$exp) <- symbols[!is.na(symbols) & !duplicated(symbols)]

motifs_new[,1] <- TFmappings[match(motifs_new[,1], TFmappings[,1]),2]

write.table(motifs_new,"~/gd/Harvard/Research/data/motifs695.txt", quote=F, col.names = F, row.names = F)
