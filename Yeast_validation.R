library(PANDA)
library(reshape2)
library(pROC)
library(gplots)
library(bereR)
data(yeast)

# data directory where all Yeast results is found (ChIP, C output, YeastNetworks_AllTFxGene.pairs, etc...)
# Actual yeast expression/motif/ppi data is loaded via PANDA package 
local.wd <- "C:/Users/Dan/Google Drive/Harvard/Research/data"
setwd(local.wd)

yeast$gold    <- read.table(file.path("YeastData_ChIP.txt"),header=F)
yeast$all     <- read.table(file.path("Data","YeastData","YeastNetworks_AllTFxGene.pairs"),header=T)
yeast.cpanda  <- read.table(file.path("yeastko_FinalNetwork.pairs"), header=F)

# Run panda algorithm
panda.res.ko <- panda(yeast$motif,yeast$exp.ko,yeast$ppi,hamming=.001,progress=F)
panda.res.cc <- panda(yeast$motif,yeast$exp.cc,yeast$ppi,hamming=.001,progress=F)
panda.res.sr <- panda(yeast$motif,yeast$exp.sr,yeast$ppi,hamming=.001,progress=F)
panda.ko <- melt(panda.res.ko@reg.net)
panda.cc <- melt(panda.res.cc@reg.net)
panda.sr <- melt(panda.res.sr@reg.net)

bere.res.ko <- bere(yeast$motif,yeast$exp.ko,score="nomotif",cpp=F)
bere.res.cc <- bere(yeast$motif,yeast$exp.cc,score="nomotif",cpp=F)
bere.res.sr <- bere(yeast$motif,yeast$exp.sr,score="nomotif",cpp=F)
bere.ko   <- melt(bere.res.ko)
bere.cc   <- melt(bere.res.cc)
bere.sr   <- melt(bere.res.sr)

bere.res.ko.wmotif <- bere(yeast$motif,yeast$exp.ko,cpp=F)
bere.res.cc.wmotif <- bere(yeast$motif,yeast$exp.cc,cpp=F)
bere.res.sr.wmotif <- bere(yeast$motif,yeast$exp.sr,cpp=F)
bere.ko.wmotif   <- melt(bere.res.ko.wmotif)
bere.cc.wmotif   <- melt(bere.res.cc.wmotif)
bere.sr.wmotif   <- melt(bere.res.sr.wmotif)

bereLDA.res.ko <- bereLDA(yeast$motif,yeast$exp.ko,cpp=F)
bereLDA.res.cc <- bereLDA(yeast$motif,yeast$exp.cc,cpp=F)
bereLDA.res.sr <- bereLDA(yeast$motif,yeast$exp.sr,cpp=F)
bereLDA.ko   <- melt(bereLDA.res.ko)
bereLDA.cc   <- melt(bereLDA.res.cc)
bereLDA.sr   <- melt(bereLDA.res.sr)

# Match predicted and gold based on sorting
# levels(yeast.panda[,1]) <- levels(yeast$gold[,1])
# levels(yeast.panda[,2]) <- levels(yeast$gold[,2])
# levels(yeast.dan[,1])   <- levels(yeast$gold[,1])
# levels(yeast.dan[,2])   <- levels(yeast$gold[,2])
# yeast.panda <- with(yeast.panda, yeast.panda[order(Var1, Var2),])

colnames(yeast$gold)   <- c("TF","GENE","gold")
colnames(yeast$motif)  <- c("TF","GENE","motif")
colnames(yeast.cpanda) <- c("TF","GENE","MOTIF","c_panda")
colnames(panda.ko)     <- c("TF","GENE","panda.ko")
colnames(panda.cc)     <- c("TF","GENE","panda.cc")
colnames(panda.sr)     <- c("TF","GENE","panda.sr")
colnames(bere.ko)      <- c("TF","GENE","bere.ko")
colnames(bere.cc)      <- c("TF","GENE","bere.cc")
colnames(bere.sr)      <- c("TF","GENE","bere.sr")

merged.yeast <- merge(yeast$gold,yeast$motif, by=c("TF","GENE"), all = TRUE)
merged.yeast <- merge(merged.yeast,yeast$all, by=c("TF","GENE"), all = TRUE)
merged.yeast <- merge(merged.yeast,yeast.cpanda, by=c("TF","GENE"), all = TRUE)

merged.yeast <- merge(merged.yeast,panda.ko, by=c("TF","GENE"), all = TRUE)
merged.yeast <- merge(merged.yeast,panda.cc, by=c("TF","GENE"), all = TRUE)
merged.yeast <- merge(merged.yeast,panda.sr, by=c("TF","GENE"), all = TRUE)
merged.yeast <- merge(merged.yeast,bere.ko, by=c("TF","GENE"), all = TRUE)
merged.yeast <- merge(merged.yeast,bere.cc, by=c("TF","GENE"), all = TRUE)
merged.yeast <- merge(merged.yeast,bere.sr, by=c("TF","GENE"), all = TRUE)

merged.yeast[is.na(merged.yeast)] <- 0
head(merged.yeast)

# png(file="all_edges.png")
# ROC curves all
dan.auc <- plot.roc(merged.yeast$gold, merged.yeast$bere.ko, col="blue", main="All edges (KO)")$auc
panda.auc <- lines.roc(merged.yeast$gold, merged.yeast$panda.ko, col="red")$auc
motif.auc <- lines.roc(merged.yeast$gold, merged.yeast$motif, col="grey")$auc
legend(.6,.2,
       c(paste("BERE",round(dan.auc,4)),paste("PANDA ",round(panda.auc,4)),paste("motif",round(motif.auc,4))),
       lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("blue","red","grey"))
# dev.off()

# png(file="motif_edges.png")
# ROC curves motif
dan.auc.motif <- plot.roc(merged.yeast$gold[merged.yeast$motif==1], merged.yeast$dan[merged.yeast$motif==1], col="blue",main="Motif edges")$auc
panda.auc.motif <- lines.roc(merged.yeast$gold[merged.yeast$motif==1], merged.yeast$KNOCK.OUT[merged.yeast$motif==1], col="red")$auc
legend(.6,.2,
       c(paste("'Classification' method",round(dan.auc.motif,4)),paste("PANDA",round(panda.auc.motif,4))),
       lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"))
# dev.off()

# png(file="non_motif_edges.png")
# ROC curves non-motif
dan.auc.nomotif <- plot.roc(merged.yeast$gold[merged.yeast$motif==0], merged.yeast$dan[merged.yeast$motif==0], col="blue",main="Non-motif edges")$auc
panda.auc.nomotif <- lines.roc(merged.yeast$gold[merged.yeast$motif==0], merged.yeast$KNOCK.OUT[merged.yeast$motif==0], col="red")$auc
legend(.6,.2,
       c(paste("'Classification' method",round(dan.auc.nomotif,4)),paste("PANDA",round(panda.auc.nomotif,4))),
       lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"))
# dev.off()

glm.motifonly <- glm(gold ~ motif, family=binomial(logit), data=merged.yeast)
glm.motifC <- glm(gold ~ motif + KNOCK.OUT, family=binomial(logit), data=merged.yeast)
glm.motifCBPP <- glm(gold ~ motif + KNOCK.OUT + dan, family=binomial(logit), data=merged.yeast)
glm.motifCR <- glm(gold ~ motif + KNOCK.OUT + Rpanda, family=binomial(logit), data=merged.yeast)
glm.motifBPP <- glm(gold ~ motif + dan, family=binomial(logit), data=merged.yeast)
summary(glm.motifCR)

anova(glm.motifonly, glm.motifC,test="Chisq")
anova(glm.motifC, glm.motifCBPP,test="Chisq")
anova(glm.motifC, glm.motifCR,test="Chisq")
