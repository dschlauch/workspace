
prev.wd <- getwd()
local.wd <- "C:/Users/Dan/Google Drive/Harvard/Research/"
setwd(local.wd)
data.dir  <- "data/Eclipse/"
setwd(data.dir)
eclipse <- list()
eclipse$motif    <- read.table("ECLIPSE_Blood_Motif.txt",header=F)
eclipse$exp      <- read.table("ECLIPSE_Blood_Exp.txt",row.names=1,header=T)
eclipse$ppi      <- read.table("OV_PPI.txt",header=F)

eclipse$clinical <- read.table("ECLIPSE_blood.txt",header=T,fill = TRUE, sep="\t",row.names=1)
eclipse$exp <- eclipse$exp[,order(colnames(eclipse$exp))]
eclipse$clinical <- eclipse$clinical[colnames(eclipse$exp),]

filter.vec.1 <- eclipse$clinical$Subject.type=="COPD Subjects"
filter.vec.2 <- eclipse$clinical$Subject.type=="Smoker Controls"

write.table(eclipse$exp[,filter.vec.1], "ECLIPSE_Blood_Exp_Subset1.txt", quote = FALSE, sep = "\t")
write.table(eclipse$exp[,filter.vec.2], "ECLIPSE_Blood_Exp_Subset2.txt", quote = FALSE, sep = "\t")

#### Note ##### PANDA expects first line to start with a pound sign (#)
# manually edit for now...

# Set up system commands, only change is the expression dataset used, along with output file
panda.exe <- "panda/Version2/PANDA"
panda.command.1 <- paste(sep="", panda.exe,
                         " -e ", data.dir, "ECLIPSE_Blood_Exp_Subset1.txt",
                         " -m ", data.dir, "ECLIPSE_Blood_Motif.txt",
                         " -p ", data.dir, "OV_PPI.txt",
                         " -o ", "eclipse_sub1_COPD")
panda.command.2 <- paste(sep="", panda.exe,
                         " -e ", data.dir, "ECLIPSE_Blood_Exp_Subset2.txt",
                         " -m ", data.dir, "ECLIPSE_Blood_Motif.txt",
                         " -p ", data.dir, "OV_PPI.txt",
                         " -o ", "eclipse_sub2_Controls")

#WARNING, very long process
system(panda.command.1, wait = FALSE)
system(panda.command.2, wait = FALSE)

### After running...

###  Load the results from cluster

regnet.copd    <- file.to.regnet(paste(data.dir,'eclipse_sub1_COPD_FinalNetwork.pairs',sep=""))
regnet.smokers <- file.to.regnet(paste(data.dir,'eclipse_sub2_Controls_FinalNetwork.pairs',sep=""))
regnet.females <- file.to.regnet(paste(data.dir,'eclipse_sub1_Female_FinalNetwork.pairs',sep=""))
regnet.males   <- file.to.regnet(paste(data.dir,'eclipse_sub2_Male_FinalNetwork.pairs',sep=""))

tm.observed <- transformation.matrix(regnet.copd, regnet.smokers)


setwd(local.wd)

# Load the null networks for Male vs Female as a list
tm.null <- load.null.tms("data/Eclipse/Eclipse_sex_null_perms","sub1")

### Same analysis as above, but for COPD vs controls
tm.null <- load.null.tms("data/Eclipse/Eclipse_COPDvsSmoker_null_perms","sub1")


# Do the sum of sq ODM plot versus null
ssodm.plot(tm.observed, tm.null,plot.title="SSODM observed and null, COPD subjects vs Smoker control")
ssodm.plot(tm.observed, tm.null, rescale=T, highlight.tfs=suspects, plot.title="SSODM observed and null, COPD subjects vs Smoker control")


# Do the sum of sq ODM plot versus null
ssodm.plot(tm.null[[11]], tm.null)
ssodm.plot(tm.null[[5]], tm.null, rescale=T)
calculate.tm.p.values(tm.eclipse.sex, tm.null,method="non-parametric")
sort(calculate.tm.p.values(tm.observed, tm.null,method="z-score"))[1:5]

p.values

######################################################
##  Running null networks with improved algorithm  ###
##                 2/25/15    START                ###
######################################################
eclipse$COPD.network <- regpredict(eclipse$motif,eclipse$exp[,filter.vec.1])
eclipse$SmCo.network <- regpredict(eclipse$motif,eclipse$exp[,filter.vec.2])

null.exp <- eclipse$exp


for(i in 1:6){
  rownames(null.exp) <- rownames(null.exp)[sample(1:nrow(null.exp))]
  eclipse[[paste(sep="","nullCOPD_",i)]] <- regpredict(eclipse$motif,null.exp[,filter.vec.1])
  eclipse[[paste(sep="","nullSmCo_",i)]] <- regpredict(eclipse$motif,null.exp[,filter.vec.2])  
}

#Parallel stuff
library(foreach)
library(doParallel)

#setup parallel backend to use 32 processors
cl<-makeCluster(2)
registerDoParallel(cl)

#start time
strt<-Sys.time()
iters <- 2
#loop
ls<-foreach(icount(iters)) %dopar% {
  
  rownames(null.exp) <- rownames(null.exp)[sample(1:nrow(null.exp))]
  res1 <- regpredict(eclipse$motif,null.exp[,filter.vec.1])
  res2 <- regpredict(eclipse$motif,null.exp[,filter.vec.2])
  list(res1,res2)
  
}

print(Sys.time()-strt)
stopCluster(cl)


tm.observed <- transformation.matrix(eclipse$COPD.network, eclipse$SmCo.network)
tm.null <- lapply(ls, function(x){
  transformation.matrix(x[[1]],x[[2]])
})
# Do the sum of sq ODM plot versus null
ssodm.plot(tm.observed, tm.null,plot.title="SSODM observed and null, COPD subjects vs Smoker control (R)")
ssodm.plot(tm.observed, tm.null, rescale=T, highlight.tfs=suspects, plot.title="SSODM observed and null, COPD subjects vs Smoker control (R)")
######################################################
##  Running null networks with improved algorithm  ###
##                 2/25/15      END                ###
######################################################

#Parallel stuff
library(foreach)
library(doParallel)

#number of iterations
iters<-1e3

#setup parallel backend to use 8 processors
cl<-makeCluster(1)
registerDoParallel(cl)

#start time
strt<-Sys.time()

#loop
ls<-foreach(icount(iters)) %dopar% {
  
  to.ls<-rnorm(1e6)
  to.ls<-summary(to.ls)
  to.ls
  
}

print(Sys.time()-strt)
stopCluster(cl)


## Gene expression analysis
library(limma)
?lmFit
design <- cbind(1,as.numeric(!filter.vec.1))
dimnames(design) <- list(colnames(eclipse$exp), c("COPD","SMOKER"))
design <- model.matrix(~factor(filter.vec.1))
diff.exp.res <- lmFit(eclipse$exp, design)
diff.exp.res <- ebayes(diff.exp.res)
diff.exp.res$p.value[order(diff.exp.res$p.value[,2])['ELK1'],]

setwd(prev.wd)


library(igraph)
# Visualize networks
net.triple <- eclipse$ppi
net.triple <- melt(top.binary)
g <- graph.data.frame(net.triple[net.triple[,3]==1,], directed=F)
V(g)$color<-ifelse(V(g)$name%in%suspects, 'blue', 'red')
l <- layout.fruchterman.reingold(g,niter=3,area=vcount(g)^2.3,repulserad=vcount(g)^5.8)
plot(g,layout=l)

g <- graph.data.frame(eclipse$motif[eclipse$motif$V1%in%c(suspects,'ELK1'),], directed=F,types=rep(0:1,590))
V(g)$color<-ifelse(V(g)$name%in%suspects, 'blue', 'red')
l <- layout.bipartite (g, types = NULL, hgap = 10, vgap = 1, maxiter = 100) 
plot(g,layout=l)

data <- eclipse$motif[eclipse$motif$V2%in%c(suspects,'ELK1','STAT3'),]
inc <- spread(data,V2,V3,fill=0)
rownames(inc) <- inc[,1]
inc<-inc[,-1]
g <- graph.incidence(top.binary)
E(g)$color<-ifelse(E(g)$grade==9, "red", "grey")
plot(g, layout=layout.bipartite,
     vertex.color=c("green","cyan")[V(g)$type+1])


# correlation matrix
genes <- unique(eclipse$motif[,1])
heatmap.2(cor(t(eclipse$exp[genes,])), dendrogram = "both", col=bluered, trace='none',cexCol=1.3,cexRow=1.3,margins=c(12,20))


# Strangely behaving TFs (because no motif data exists for these TFs)
suspects <- c('HIF1A',
              'HAND1',
              'ARID3A',
              'SOX5',
              'PRRX2',
              'NR1H2',
              'POU5F1',
              'VDR',
              'NFE2L1',
              'NKX3-1',
              'MAFG',
              'AHR',
              'SOX10',
              'DDIT3',
              'TAL1',
              'EWSR1',
              'NFIC',
              'TLX1',
              'ARNT'
)
