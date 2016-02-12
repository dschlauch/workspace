library(reshape2)
library(bereR)
library(pandaR)
library(ROCR)
library(dplyr)
library(penalized)

validateMethodsOnDataset <- function(dataset){
    source(paste("~/gd/Harvard/Research/R_workspace/",dataset,".R",sep=""))
    
    # Data Procession
    # Remove non-target genes
    exprData <- exprData[rownames(exprData)%in%unique(goldStandard[,2]),]
    # order Genes
    exprData <- scale(exprData[order(rownames(exprData)),])
    if (ncol(exprData)>200){
        # Keep number of samples under 200
        exprData <- exprData[,1:100]
    }
    # Get subset of expression data that is relevant (possibly keep all)
#     dataset <- "DREAM5c_G313"
#     exprData <- exprData[,!grep("G313",chipFeatures$V5)]
    
    # Run algorithms -  PANDA, BERE (with Corr method, LDA, and weighted correlation diff), straight up corellation
    
    ########## PANDA
#     pandaMelt <- panda(motifs, exprData, progress=F)@regNet %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt

    ########## Degree only approach
#     degreeMelt <- degreeApproach(motifs) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
    
    ##########  BERE
    bereMelt <- bere(motifs, exprData, cpp=F, verbose=F) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
    
    ##########  LDA BERE
#     ldabereMelt <- ldaBERE(motifs, exprData) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt

    ##########  FULL BERE
#     fullbereMelt1 <- bereFull(motifs, exprData, alpha=1, lambda=1) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
#     fullbereMelt2 <- bereFull(motifs, exprData, alpha=1, lambda=2) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
#     fullbereMelt5 <- bereFull(motifs, exprData, alpha=1, lambda=5) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
#     fullbereMelt10 <- bereFull(motifs, exprData, alpha=1, lambda=10) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
#     fullbereMelt15 <- bereFull(motifs, exprData, alpha=1, lambda=15) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
#     fullbereMelt25 <- bereFull(motifs, exprData, alpha=1, lambda=25) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
#     fullbereMelt40 <- bereFull(motifs, exprData, alpha=1, lambda=40) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
#     fullbereMelt1000 <- bereFull(motifs, exprData, alpha=1, lambda=1000) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt

    
    ##########  Straight TF correlation
#     tfCorMelt <- abs(cor(x=t(exprData[rownames(exprData) %in% transFactors,]), y=t(exprData))) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
    
    ##########  Weighted Cor diff
#    weightedCorDiffMelt <- t(weightedCorDiff) %>% melt %>% meltToCharacter %>% removeDiagonal %>% sortMelt
    
    ###################################################  
    ##########  Validate Results with gold standard
    ###################################################
    goldStandard <- goldStandard %>% meltToCharacter %>% removeDiagonal %>% sortMelt
    motifs <- motifs %>% meltToCharacter %>% removeDiagonal %>% sortMelt

    ###################################################  
    ##########  Plot Results against gold standard
    ###################################################
    datalist <- list("Gold Standard"=goldStandard[,3],
#                      "Degree-only"=degreeMelt[,3],
#                      "PANDA"=pandaMelt[,3],
                    "BERE"=bereMelt[,3],
#                     "LDA BERE"=ldabereMelt[,3],
#                     "Full BERE 1"=fullbereMelt1[,3],
#                     "Full BERE 2"=fullbereMelt2[,3],
#                     "Full BERE 5"=fullbereMelt5[,3],
#                     "Full BERE 10"=fullbereMelt10[,3],
#                     "Full BERE 15"=fullbereMelt15[,3],
#                     "Full BERE 25"=fullbereMelt25[,3],
#                     "Full BERE 40"=fullbereMelt40[,3],
#                     "Full BERE 1000"=fullbereMelt1000[,3],
#                     "TF Corr"=tfCorMelt[,3],#+motifs[,3],
                    "Motifs"=motifs[,3])
    TFsubset <- goldStandard[,1] %in% unique(goldStandard[,1])[6]
    
    png(filename=paste("./TM_manuscript/figures/",dataset,"_all.png",sep=""))
    plotROC(datalist, "all", organism=dataset, goldStandard)
    dev.off()
    png(filename=paste("./TM_manuscript/figures/",dataset,"_motif.png",sep=""))
    plotROC(datalist, "motif", organism=dataset, goldStandard)
    dev.off()
    png(filename=paste("./TM_manuscript/figures/",dataset,"_nonmotif.png",sep=""))
    plotROC(datalist, "nonmotif", organism=dataset, goldStandard)
    dev.off()
}



###################################################  
##########  Functions
###################################################
# New ROC method
plotROC <- function(datalist, includeSubset="all", organism="", goldStandard=NA, TFsubset=NA){
    require(ROCR)
    
    methods <- names(datalist)[-1] # Remove gold standard from methods list
    if(includeSubset=="all"){ subset <- rep(T,length(datalist[["Gold Standard"]]))}
    if(includeSubset=="motif"){ subset <- (datalist[["Motifs"]]==1)}
    if(includeSubset=="nonmotif"){ subset <- (datalist[["Motifs"]]==0)}
    
    if(!all(is.na(TFsubset))){
        subset <- subset*TFsubset==1
    }
    plotList <- lapply(methods, function(x){
        methodPred  <- prediction(datalist[[x]][subset], datalist[["Gold Standard"]][subset])
        roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
        auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
        list("roc.methodPred"=roc.methodPred, "auc.methodPred"=auc.methodPred)
    })
    names(plotList) <- methods
    plot(plotList[["BERE"]][["roc.methodPred"]], main=paste(organism, "-", includeSubset, "motifs","ROC"), col = 1, lwd=3)
    mapply(function(x,index){
        lines(plotList[[x]][["roc.methodPred"]]@x.values[[1]], plotList[[x]][["roc.methodPred"]]@y.values[[1]], col = (index), lwd=3)
    }, methods, 1:length(methods))
    legendLabels <- c(sapply(methods, function(x){
        paste(x, round(plotList[[x]][["auc.methodPred"]],4))
    }))
    

    legend(.5,.6, legendLabels, lty=rep(1,length(methods)),lwd=rep(5,length(methods)),col=1:length(methods),title="Area under ROC curve")
    rocPerTFResults <- c()
    if(!all(is.na(goldStandard))){
        # Calculate average ROC per TF
        tfs <- unique(goldStandard[,1])
        goldcounts <- table(goldStandard[subset,c(1,3)])[,2]
        tfs <- tfs[goldcounts>0]
        goldcounts <- goldcounts[goldcounts>0]
        
        TFaucrocs <- sapply(tfs, function(tf){
            tfsub <- goldStandard[,1] %in% tf
            sapply(methods, function(x){
                methodPred  <- prediction(datalist[[x]][subset & tfsub], datalist[["Gold Standard"]][subset & tfsub])
                performance(methodPred, "auc")@y.values[[1]]
            })
        })
        meanROC <- apply(TFaucrocs,1,mean)
        wmeanROC <- apply(TFaucrocs,1,function(x){
            sum(x*goldcounts)/sum(goldcounts)
        })
        rocPerTFResults <- c(rocPerTFResults, paste("Average within TF ROC, weighted by gold standard hits -", organism, includeSubset))
        rocPerTFResults <- c(rocPerTFResults, wmeanROC)
    }
    write.table(rocPerTFResults, file=paste("./output/",organism, includeSubset,"_aucroc_per_TF.txt",sep=""))
}
sortMelt <- function(df){df[order(df[,1],df[,2]),]}
removeDiagonal <- function(x){x[x[,1]!=x[,2],]}
meltToCharacter <- function(x){
    x[,1]<-as.character(x[,1])
    x[,2]<-as.character(x[,2])
    x
}