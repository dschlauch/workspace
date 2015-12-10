library(bereR)
library(pandaR)
library(bptools)
library(reshape2)
library(penalized)
library(Biobase)

analysisCode <- sample(100000,1)

# copd.filename <- "~/gd/Harvard/Research/data/Eclipse/null.networks_all.rds"
# dataset.filename <- "~/gd/Harvard/Research/data/Eclipse/eclipse.networks.rds"

motifFile <- "~/gd/Harvard/Research/data/Eclipse/ECLIPSE_Blood_Motif.txt"
# exprFile <- "~/gd/Harvard/Research/data/Ovarian/CombinedOV.txt"
exprFile <- "~/gd/Harvard/Research/data/GTEx/GTEx_expr.txt"
# exprFile <- "~/gd/Harvard/Research/data/COPDGene/COPDGene_GSExpressionData.txt"
#exprFile <- "~/gd/Harvard/Research/data/Eclipse/ECLIPSE_Blood_Exp.txt"
#exprFile <- "~/gd/Harvard/Research/data/LGRC/LGRC_expression.txt"
ppiFile <- "~/gd/Harvard/Research/data/Eclipse/OV_PPI.txt"
#clinicalFile <- "~/gd/Harvard/Research/data/Ovarian/Clinical.txt"
clinicalFile <- "~/gd/Harvard/Research/data/GTEx/GTEx_clinical.txt"
#clinicalFile <- "~/gd/Harvard/Research/data/COPDGene/COPDGene_clinical.txt"
#clinicalFile <- "~/gd/Harvard/Research/data/Eclipse/ECLIPSE_blood.txt"
#clinicalFile <- "~/gd/Harvard/Research/data/LGRC/lgrc.merged.clinical.data.clean.txt"

# exprFile <- "~/gd/Harvard/Research/data/GTEx/gtex_cell.rdata"
exprFile <- "~/gd/Harvard/Research/data/GTEx/gtex_sub_noxymt_qsmooth_cell.rdata"
motifFile <- "~/gd/Harvard/Research/data/GTEx/KG_cisbp_652.txt"
phenotypeName <- "our_subtypes"
casesString <- "cells_ebv-transformed_lymphocytes"
controlsString <- "skin"
filterType <- NA

# casesString <- "COPD Subjects"
# controlsString <- "Smoker Controls"
# phenotypeName <- "Subject.type"

casesString <- "Lung"
controlsString <- "Colon"
phenotypeName <- "SMTS"

analysisName <- "ECLIPSE_matchMF"
nullPerms <- 500
networkInferenceName <- "bere"
filterType <- "Gender"
filterBy <- "M"
permuteGeneLabels <- F

args<-commandArgs(TRUE)
if(length(args)!=0){
    motifFile <- args[1]
    exprFile <- args[2]
    ppiFile <- args[3]
    clinicalFile <- args[4]
    casesString <- args[5]
    controlsString <- args[6]
    phenotypeName <- args[7]
    analysisName <- args[8]
    nullPerms <- as.numeric(args[9])
    networkInferenceName <- args[10]
    # Optional parameters
    if(length(args)>10){
        filterType <- args[11]
        filterBy <- args[12]
        if(filterType=="NA"){
            filterType <- NA
        }
    } else {
        filterType <- NA
    }
    if(length(args)>12){
        permuteGeneLabels <- (args[13]=="T")
    } else {
        permuteGeneLabels <- F
    }
}
outputDir <- file.path("~",paste(analysisName, analysisCode ,sep="_"))
dir.create(outputDir, showWarnings=FALSE)
write.csv(args,file=file.path(outputDir,paste("arguments_",analysisCode,".txt",sep="")))

# Set the network inference method
if(networkInferenceName=="bere"){
    networkInferenceMethod <- bere  
}
if(networkInferenceName=="bereFull"){
    networkInferenceMethod <- bereFull
}
if(networkInferenceName=="panda"){
    networkInferenceMethod <- function(motifs, exp){
        panda(motifs, exp)@regNet
    }
}
if(networkInferenceName=="pandaC"){
    networkInferenceMethod <- function(motifs, exp){
        randindex <- sample(10000000000,1)
        # write the files to disk
        motifFile <- file.path("~/tmp", paste0("_motifs_tmp",randindex,".txt"))
        expFile   <- file.path("~/tmp", paste0("_exp_tmp",randindex,".txt"))
        outFile   <- paste0("output", randindex)
        write.table(motifs, motifFile, sep="\t", quote=F, col.names = F,  row.names = F)
        write.table(exp, expFile, sep="\t", quote=F, col.names = F)
        
        system(paste0("./gd/Harvard/Research/panda/Version2/PANDA ",
                      "-e ", expFile, " ",
                      "-m ", motifFile, " ",
                      "-o ", outFile, " ",
                      "-k 4"
        ))
        
        # read in results
        regnet <- read.table(paste0(outFile, "_FinalNetwork.pairs"), header=F)
        regnet <- dcast(regnet, V1 ~ V2, value.var='V4')
        rownames(regnet) <- regnet[,1]
        regnet <- regnet[,-1]
        
        # clean up
        system(paste0('rm ',motifFile))
        system(paste0('rm ',expFile))
        system(paste0('rm ',outFile, "_FinalNetwork.pairs"))
        as.matrix(regnet)
    }
}
if(networkInferenceName=="pandaM"){
    networkInferenceMethod <- function(motifs, exp){
        randindex <- sample(10000000000,1)
        
        # write the files to disk
        motifFilename <- paste("_motifs_tmp",randindex,".txt",sep="")
        expFilename <- paste("_exp_tmp",randindex,".txt",sep="")
        outFile   <- paste0("output", randindex)
        write.table(motifs, file.path("/scratch/",motifFilename), sep="\t", quote=F, col.names = F,  row.names = F)
        write.table(exp, file.path("/scratch/",expFilename), sep="\t", quote=F, col.names = F)
        
        # Copy master panda script
        system(paste('cp ','~/panda_matlab/RunPANDA.m ','~/panda_matlab/RunPANDA', randindex,'.m',sep=""))
        
        # Set file pointers
        system(paste0("sed -i 's/motifPlaceholder.txt/", motifFilename, "/g' ~/panda_matlab/RunPANDA", randindex,".m"))
        system(paste0("sed -i 's/expPlaceholder.txt/", expFilename, "/g' ~/panda_matlab/RunPANDA", randindex,".m"))
        system(paste0("sed -i 's/PANDAOutputPlaceholder/", outFile, "/g' ~/panda_matlab/RunPANDA", randindex,".m"))
        
        # Run matlab script
        system(paste0("matlab -nojvm -nodesktop -r 'run  ~/panda_matlab/RunPANDA", randindex,".m;quit'"))
        
        # clean up (removed after switching to scrarch space)
        #         system(paste0('rm ',file.path("/scratch", motifFilename)))
        #         system(paste0('rm ',file.path("/scratch", expFilename)))
        system(paste0('rm ',file.path("~/panda_matlab/", paste0("RunPANDA",randindex,".m"))))
        
        # load results back into R
        # read in results
        regnet <- read.table(paste0("/scratch/", outFile, "_FinalNetwork.pairs"), header=T)
        #         system(paste0('rm ',file.path("~/panda_matlab", paste0(outFile, "_FinalNetwork.pairs"))))
        regnet <- dcast(regnet, TF ~ gene, value.var='PANDA.prediction')
        rownames(regnet) <- regnet[,1]
        regnet <- regnet[,-1]
        res <- as.matrix(regnet)
        if(sum(is.na(res))>0){
            saveRDS(list(exp,res,motifs),paste0("FailedNetwork",randindex,".rda"))
        }
        res
    }
}
######################################################
##      Data Loading from ECLIPSE dataset          ###
##                                                 ###
######################################################
dataset <- list()
if (grepl(".txt", exprFile)){
    dataset$motif    <- read.table(motifFile,header=F)
    dataset$exp      <- read.table(exprFile,row.names=1,header=T)
    dataset$clinical <- read.table(clinicalFile,header=T,fill = TRUE, sep="\t",row.names=1)
    
    # Removed this substring line for GTEx data (may need to reinsert for some other dataset)
    colnames(dataset$exp) <- substr(colnames(dataset$exp), 1, 10)
    rownames(dataset$clinical) <- substr(rownames(dataset$clinical), 1, 10)
} else if (grepl(".rdata", exprFile)){
    #GTEx analysis
    load(exprFile)
    dataset$motif    <- cbind(read.table(motifFile,header=F),1)
    # 11/14/15 changed "both" to "obj" for camilla dataset
    dataset$exp      <- exprs(obj)
    dataset$clinical <- pData(obj)
    
    # Remove ensembl decimal and value after
    rownames(dataset$exp) <- substring(rownames(dataset$exp),1,15)
    
    #Get top 20,000 variable genes
    rowsds <- sort(apply(dataset$exp, 1, sd), decreasing=T)
    genesIncluded <- names(rowsds[1:19000])
    dataset$exp <- dataset$exp[genesIncluded,]
    dataset$motif <- dataset$motif[dataset$motif[,2]%in%genesIncluded,]
    
    mappingFile <- "~/gd/Harvard/Research/data/GTEx/cisbpall_motinf.txt"
    mappings <- read.table(mappingFile, header=T)
    mappings[,1] <- substring(mappings[,1],0,5) 
    dataset$motif[,1] <- mappings[match(dataset$motif[,1], mappings[,1]),2]
    dataset$motif <- dataset$motif[!is.na(dataset$motif[,1]),]
    
    library(org.Hs.eg.db)
    symbols <- mapIds(org.Hs.eg.db, keys=row.names(dataset$exp),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
    dataset$exp <- dataset$exp[!is.na(symbols) & !duplicated(symbols),]    
    rownames(dataset$exp) <- symbols[!is.na(symbols) & !duplicated(symbols)]
    
}
dataset$ppi      <- read.table(ppiFile,header=F)
dataset$exp      <- dataset$exp[,order(colnames(dataset$exp))]  # Make sure expression and clinical is in same order


matches <- sort(unique(intersect(rownames(dataset$clinical),colnames(dataset$exp))))
dataset$clinical <- dataset$clinical[matches,]    # Make sure clinical only contains patients with expression data
dataset$exp <- dataset$exp[,matches]    # Make sure expression only contains patients with clinical data

if(permuteGeneLabels){
    print("Permuting gene labels once")
    rownames(dataset$exp) <- sample(rownames(dataset$exp))d
} else {
    print("No gene label permutation (default)")
}

# Specify the group partition
if(is.na(filterType)){
    subsetFilter <- rep(T,nrow(dataset$clinical))
} else {
    subsetFilter <- dataset$clinical[,filterType]==filterBy
}
phenoFilter <- (dataset$clinical[,phenotypeName]==casesString)|(dataset$clinical[,phenotypeName]==controlsString)
allFilter <- subsetFilter&phenoFilter

dataset$exp <- dataset$exp[,allFilter]
dataset$clinical <- dataset$clinical[allFilter,]
casesFilter <- dataset$clinical[,phenotypeName]==casesString
controlsFilter <- dataset$clinical[,phenotypeName]==controlsString

# covariateValues <- dataset$clinical[,covariate]
# table(dataset$clinical$pkyrs>40, casesFilter)
# mean(dataset$clinical$pkyrs[casesFilter])
# mean(dataset$clinical$pkyrs[!casesFilter])
# male <- dataset$clinical$GENDER=="1-Male"
# mean(dataset$clinical$pkyrs[male])
# mean(dataset$clinical$pkyrs[!male])
# 
# table(controlsFilter,dataset$clinical[,"Gold.stage"])

######################################################
##  Running null networks with improved algorithm  ###
##                 2/25/15    START                ###
######################################################

# dataset$casesNetwork <- networkInferenceMethod(dataset$motif,dataset$exp[,casesFilter])
# dataset$controlsNetwork <- networkInferenceMethod(dataset$motif,dataset$exp[,controlsFilter])
# 
# # periodically save workspace
# save.image(file=file.path(outputDir,paste("activeImage",analysisCode,".RData",sep="")))

# Copy expression data for null network generation
null.exp <- dataset$exp

#Parallel stuff
library(foreach)
library(doParallel)

# Calculate the number of cores
num_cores <- detectCores() - 4
num_cores <- min(num_cores, 40)

# Initiate cluster
if(!is.na(num_cores)){
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
}

#start time
strt  <- Sys.time()
iters <- nullPerms*2+2 # Two networks for each partition, plus observed partition
#loop
print("Running null permutations in parallel")
print(paste0(num_cores," cores used"))
print(paste0(iters," networks to be estimated"))
null.networks<-foreach(i=1:iters,.packages=c("bereR","pandaR","reshape2","penalized")) %dopar% {
    print(paste0("Running iteration ", i))
    if(i%%2==0){
        selectedSamples <- casesFilter
    } else {
        selectedSamples <- controlsFilter
    }
    if(i<=2){
        # Observed partition : Don't reorder anything
        null.exp <- dataset$exp
    } else {
        # Null partition, randomly reorder
        ## resample case-control
        null.exp <- dataset$exp[,sample(1:ncol(dataset$exp))]
        ## This line scrambles the gene names (toggle this) 8/18/15
        rownames(null.exp) <- rownames(null.exp)[sample(1:nrow(null.exp))]
    }
    #     null.exp <- null.exp + matrix(rnorm(length(null.exp))/10,nrow=nrow(null.exp),ncol=ncol(null.exp))
    null.exp <- null.exp[,selectedSamples]
    if (sum(rowSums(null.exp)==0)>0){
        zeroGenes <- which(rowSums(null.exp)==0)
        for(gene in zeroGenes){
            null.exp[gene,] <- rnorm(ncol(null.exp))
        }
    }
    tmpNet <- networkInferenceMethod(dataset$motif, null.exp)
    print(paste0("Finished running iteration", i))
    tmpNet
}

print(Sys.time()-strt)
if(!is.na(num_cores)){
    stopCluster(cl)
}




# Add new null permutations to existing list if any
# This step is to allow for skipping the above step and starting with a stored null set
# if (file.exists(copd.filename)){
# null.networks <- append(null.networks, readRDS("null.networks_all.rds"))
# }

# Save the observed and null networks (as separate files)
#saveRDS(dataset,dataset.filename)
#saveRDS(null.networks,copd.filename)

#####################################################
# START HERE TO SKIP PERMUTATIONS.
#####################################################

#null.networks  <-  readRDS(copd.filename)
#dataset        <-  readRDS(dataset.filename)

#####################################################
###  TF analysis
#####################################################

if(!is.na(num_cores)){
    cl <- makeCluster(4)
    registerDoParallel(cl)
}

strt  <- Sys.time()
#loop
print("Running transition calculations in parallel")
print(paste0(num_cores," cores used"))
print(paste0(length(null.networks)/2," transitions to be estimated"))
transMatrices <- foreach(i=1:(length(null.networks)/2),.packages=c("bptools","reshape2","penalized")) %dopar% {
    transformation.matrix(null.networks[[2*i]], null.networks[[2*i-1]],remove.diagonal=T,method="ols")    
}

print(Sys.time()-strt)
if(!is.na(num_cores)){
    stopCluster(cl)
}

# Parallelized this part on 10/30/15
# # Calculate the transformation matrix for the observed data
# tm.observed <- transformation.matrix(null.networks[[2]], null.networks[[1]],remove.diagonal=T,method="ols")
# 
# 
# # Calculate the transformation matrix for the null data
# tm.null <- lapply(1:nullPerms, function(x){
#     transformation.matrix(null.networks[[2*x+2]],null.networks[[2*x+1]],method="ols",remove.diagonal = T)
# })

# dataset$controlsNetwork <- null.networks[1]
# dataset$casesNetwork <- null.networks[2]

# This object will be in the many GB range
rm(null.networks)
gc()

# periodically save workspace
save.image(file=file.path(outputDir,paste("activeImage",analysisCode,".RData",sep="")))


# Do the sum of sq ODM plot versus null
png(file.path(outputDir,paste('SSODMplot_unscaled',analysisCode,'.png', sep="")), width=4800)
ssodm.plot(transMatrices[[1]], transMatrices[-1],plot.title=paste("SSODM observed and null, ",casesString," vs ",controlsString,' : ', networkInferenceName, ' : ', analysisName, sep=""))
dev.off()
png(file.path(outputDir,paste('SSODMplot_scaled',analysisCode,'.png', sep="")), width=4800)
ssodm.plot(transMatrices[[1]], transMatrices[-1], rescale=T, plot.title=paste("SSODM observed and null, ",casesString," vs ",controlsString,' : ', networkInferenceName, ' : ', analysisName, sep=""))
dev.off()

# This is for converting TF IDs to their gene names

# transMatrices <- lapply(transMatrices, function(x){
#     rownames(x) <- mappings[match(rownames(x), mappings[,1]),2]
#     colnames(x) <- mappings[match(colnames(x), mappings[,1]),2]
#     x
# })

# Top TFs
#highlight.tfs <- c("E2F4","NRF1","GABPA","ELK1","ELK4","E2F1","ZBTB33","ELF1","ZFX")

# #####################################################
# ### Gene Analysis
# #####################################################
# # Calculate the transformation matrix for the observed data
# tm.observed.genes <- transformation.matrix(dataset$casesNetwork, dataset$controlsNetwork,remove.diagonal=T,method="old",by.tfs=T,standardize=F)
# 
# 
# # Calculate the transformation matrix for the null data
# tm.null.genes <- lapply(null.networks, function(x){
#   transformation.matrix(x[[1]],x[[2]],method="old",standardize=F)
# })
# 
# # Do the sum of sq ODM plot versus null
# ssodm.plot(tm.observed, tm.null,plot.title="SSODM observed and null, COPD subjects vs Smoker control (R)")
# ssodm.plot(tm.observed, tm.null, rescale=T, plot.title="SSODM observed and null, COPD subjects vs Smoker control (R)",highlight.tfs = c("ELK1","E2F4"))



######################################################
##  Running null networks with improved algorithm  ###
##                 2/25/15      END                ###
######################################################
# 
# 
# cppFunction('NumericMatrix correl(NumericMatrix x) {
#   int nrow = x.nrow(), ncol = x.ncol();
#   NumericMatrix resultMatrix(nrow,nrow);
#   for (int i = 0; i < nrow; i++) {
#     for (int j = i; j < nrow; j++) {
#       double sumproduct = 0;
#       for (int k = 0; k < ncol; k++){
#         sumproduct+=x(i,k)*x(j,k);
#       }
#       resultMatrix(i,j) = sumproduct/(nrow-1);
#       resultMatrix(j,i) = sumproduct/(nrow-1);
#     }
#   }
#   return resultMatrix;
# }')
# 
# cppFunction('double squarert(double x){
#   return sqrt(x);
#             }')
# set.seed(1014)
# randmat <- t(matrix(rnorm(30000*50), 50))
# t(apply(randmat, 1, function(x)(x-mean(x))/(sd(x))))
# 
# system.time(correl(t(apply(randmat, 1, function(x)(x-mean(x))/(sd(x))))))
# system.time(cor(t(randmat)))
# #>  [1] 458 558 488 458 536 537 488 491 508 528
# rowSumsC(x)
# #>  [1] 458 558 488 458 536 537 488 491 508 528
# 

## Gene expression analysis
library(limma)
design <- model.matrix(~factor(casesFilter))
diff.exp.res <- lmFit(dataset$exp, design)
diff.exp.res <- ebayes(diff.exp.res)
logfoldchange <- log(rowMeans(dataset$exp[,casesFilter])/rowMeans(dataset$exp[,controlsFilter]))


# 7/28/15 
# create results table
# 10/30/15 updates for GTEx, which does not (or I'm not using) expression values for TFs
obsSsodm <- apply(transMatrices[[1]],1,function(x){t(x)%*%x})
includedTFs <- intersect(names(obsSsodm),rownames(diff.exp.res$p.value))
obsSsodm <- obsSsodm[includedTFs]
dTFI_pVals_All <- 1-2*abs(.5-calculate.tm.p.values(transMatrices[[1]], transMatrices[-1]))
dTFI_pVals <- dTFI_pVals_All[includedTFs]
negLogPValues <- -log(dTFI_pVals)
# replace Inf values with max values
negLogPValues[negLogPValues==Inf] <- 35
labels <- names(obsSsodm)
labels[rank(-negLogPValues)>20 & rank(-obsSsodm)>20]<-""

dTFI_fdr   <- p.adjust(dTFI_pVals, method = 'fdr')
logfoldchangeTF <- logfoldchange[names(dTFI_pVals)]
limma_pVals <- diff.exp.res$p.value[names(dTFI_pVals),2]
limma_fdr <- p.adjust(limma_pVals, method = 'fdr')
limmanegLogPValues <- -log(limma_pVals)
limmanegLogPValues[limmanegLogPValues>10]<-10 # for visual purposes
resultTable <- cbind(obsSsodm,dTFI_pVals,dTFI_fdr,logfoldchangeTF,limmanegLogPValues)
resultTable <- resultTable[order(dTFI_pVals),]

plotDF <- data.frame(obsSsodm, negLogPValues, limmanegLogPValues, logfoldchangeTF, "labels"=labels)

colnames(resultTable) <- c("Magnitude","dTFI uncorrected p-value","dTFI FDR", "log FC", "limma -logp")
write.csv(resultTable,file=file.path(outputDir,paste("resultTable",analysisCode,".csv", sep="")))

png(file.path(outputDir,paste('Volcano plot',analysisCode,'.png', sep="")), width=1800)
ggplot(data=plotDF,aes(x=obsSsodm, y=negLogPValues, label=labels, size=100)) + geom_point(aes(col=limmanegLogPValues), alpha=.5) + geom_text(vjust=0) + 
    ylab("-log(dTFI p-value)") + xlab("SSODM") + ggtitle("Signal vs significance") + scale_size_continuous(range = c(0, 12),guide=FALSE) +
    scale_colour_gradientn("LIMMA sig",colours=c("blue","white","red"))
dev.off()

png(file.path(outputDir,paste('dTFI vs LIMMA',analysisCode,'.png', sep="")), width=1800)
ggplot(data=plotDF,aes(x=logfoldchangeTF, y=negLogPValues, label=labels, size=100)) + geom_point(aes(col=limmanegLogPValues), alpha=.75) + geom_text(vjust=0) + 
    ylab("significance of dTFI") + xlab("-log(FC)") + ggtitle("Transition vs Fold Change") +
    scale_size_continuous(range = c(0, 12),guide=FALSE) +
    scale_colour_gradientn("LIMMA sig",colours=c("blue","white","red"))
dev.off()

# periodically save workspace
save.image(file=file.path(outputDir,paste("activeImage",analysisCode,".RData",sep="")))

require(igraph)
## Calculate p-values for off-diagonals
transitionSigmas <- function(tm.observed, tm.null){
    tm.null.mean <- apply(simplify2array(tm.null), 1:2, mean)
    tm.null.sd <- apply(simplify2array(tm.null), 1:2, sd)
    sigmas <- (tm.observed - tm.null.mean)/tm.null.sd
}

tm.sigmas <- transitionSigmas(transMatrices[[1]], transMatrices[-1])
diag(tm.sigmas) <- 0
tm.sigmas.melt <- melt(tm.sigmas)

adjMat <- transMatrices[[1]]
diag(adjMat) <- 0
adjMat.melt <- melt(adjMat)

adj.combined <- merge(tm.sigmas.melt,adjMat.melt, by=c("Var1","Var2"))

# adj.combined[,1] <- mappings[match(adj.combined[,1], mappings[,1]),2]
# adj.combined[,2] <- mappings[match(adj.combined[,2], mappings[,1]),2]

numEdges  <- 40
numTopTFs <- 10
topTFsIncluded <- names(sort(dTFI_pVals_All)[1:numTopTFs])
topTFIndices <- 2>(is.na(match(adj.combined[,1],topTFsIncluded))+is.na(match(adj.combined[,2],topTFsIncluded)))
adj.combined <- adj.combined[topTFIndices,]
adj.combined <- adj.combined[abs(adj.combined[,4])>=sort(abs(adj.combined[,4]),decreasing=T)[numEdges],]
tfNet <- graph.data.frame(adj.combined, directed=T)
vSize <- -log(dTFI_pVals_All)
vSize[vSize<0] <- 0
vSize[vSize>3] <- 3

V(tfNet)$size <- vSize[V(tfNet)$name]*5
E(tfNet)$width <- (abs(E(tfNet)$value.x))*20/max(abs(E(tfNet)$value.x))
E(tfNet)$color<-ifelse(E(tfNet)$value.x>0, "blue", "red")

png(file.path(outputDir,paste('Transition plot',analysisCode,'.png', sep="")), width=1500, height=1500)
plot.igraph(tfNet, edge.arrow.size=2, vertex.label.cex= 1.5, vertex.label.color= "black",main=paste0("Transition: ",controlsString," to ",casesString))
legend("bottomleft", c("Gained features","Lost features"), lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"))
dev.off()

saveRDS(list(obsSsodm,dTFI_pVals_All),'dTFI.rdata')


# ############# Some methylation investigation 11/25/15
# # diffGeneMethPValues <- readRDS("./diffGeneMethPValues.rdata")
# meanTransition <- rowMeans(transMatrices[[1]])
# methDifferences <- diffGeneMethPValues[1,]-diffGeneMethPValues[2,] 
# matched_TF_Meth <- intersect(colnames(diffGeneMethPValues),names(obsSsodm))
# plot(methDifferences[matched_TF_Meth],meanTransition[matched_TF_Meth])
# summary(lm(methDifferences[matched_TF_Meth]~meanTransition[matched_TF_Meth]))

#library(igraph)
# # Visualize networks
# net.triple <- dataset$ppi
# net.triple <- melt(top.binary)
# g <- graph.data.frame(net.triple[net.triple[,3]==1,], directed=F)
# V(g)$color<-ifelse(V(g)$name%in%suspects, 'blue', 'red')
# l <- layout.fruchterman.reingold(g,niter=3,area=vcount(g)^2.3,repulserad=vcount(g)^5.8)
# plot(g,layout=l)
# 
# g <- graph.data.frame(dataset$motif[dataset$motif$V1%in%c(suspects,'ELK1'),], directed=F,types=rep(0:1,590))
# V(g)$color<-ifelse(V(g)$name%in%suspects, 'blue', 'red')
# l <- layout.bipartite (g, types = NULL, hgap = 10, vgap = 1, maxiter = 100) 
# plot(g,layout=l)
# 
# data <- dataset$motif[dataset$motif$V2%in%c(suspects,'ELK1','STAT3'),]
# inc <- spread(data,V2,V3,fill=0)
# rownames(inc) <- inc[,1]
# inc<-inc[,-1]
# g <- graph.incidence(top.binary)
# E(g)$color<-ifelse(E(g)$grade==9, "red", "grey")
# plot(g, layout=layout.bipartite,
#      vertex.color=c("green","cyan")[V(g)$type+1])
# 
# 
# # correlation matrix
# genes <- unique(dataset$motif[,1])
# heatmap.2(cor(t(dataset$exp[genes,])), dendrogram = "both", col=bluered, trace='none',cexCol=1.3,cexRow=1.3,margins=c(12,20))
# 
# 
# 
# # Strangely behaving TFs (because no motif data exists for these TFs)
# # Note 2/26/15  These should be filtered out of the analysis in this version of PandaR
# suspects <- c('HIF1A',
#               'HAND1',
#               'ARID3A',
#               'SOX5',
#               'PRRX2',
#               'NR1H2',
#               'POU5F1',
#               'VDR',
#               'NFE2L1',
#               'NKX3-1',
#               'MAFG',
#               'AHR',
#               'SOX10',
#               'DDIT3',
#               'TAL1',
#               'EWSR1',
#               'NFIC',
#               'TLX1',
#               'ARNT'
# )
# 
# 
# #################################
# ###  Genes approach
# ##################################
# 
# obs.s  <- apply(dataset$controlsNetwork-dataset$casesNetwork,2,function(x){x%*%x})
# null.s <- lapply(null.networks, function(nullnet){
#   apply(nullnet[[1]]-nullnet[[2]],2,function(x){x%*%x})
# })
