# Specify study and motif in command line or before sourcing
# study <- "LGRC"
# study <- "COPDGENE"
# study <- "LTCOPD"
# study <- "ECLIPSE"
# 
# motifVersion <- "JASPAR2016"
# motifVersion <- "CISBP"
# motifVersion <- "JASPAR2014"

if (motifVersion=="JASPAR2014"){
  motifFile <- "~/gd/Harvard/Research/data/Eclipse/ECLIPSE_Blood_Motif.txt"  
}
if (motifVersion=="JASPAR2016"){
  motifFile <- "does not exist yet"  
}
if (motifVersion=="CISPB"){
  motifFile <- "~/gd/Harvard/Research/data/motifs695.txt"
}

if (study=="ECLIPSE"){
  exprFile <- "~/gd/Harvard/Research/data/Eclipse/ECLIPSE_Blood_Exp.txt"
  clinicalFile <- "~/gd/Harvard/Research/data/Eclipse/ECLIPSE_blood.txt"
  #ECLIPSE labels
  casesString <- "COPD Subjects"
  controlsString <- "Smoker Controls"
  phenotypeName <- "Subject.type"
} else if (study=="COPDGENE"){
  exprFile <- "~/gd/Harvard/Research/data/COPDGene/COPDGene_GSExpressionData.txt"
  clinicalFile <- "~/gd/Harvard/Research/data/COPDGene/COPDGene_clinical.txt"
  #COPDGene labels
  casesString <- "COPD Subjects"
  controlsString <- "Smoker Controls"
  phenotypeName <- "Subject.type"
} else if (study=="LGRC"){
  exprFile <- "~/gd/Harvard/Research/data/LGRC/LGRC_expression.txt"
  clinicalFile <- "~/gd/Harvard/Research/data/LGRC/lgrc.merged.clinical.data.clean.txt"
  #LGRC labels
  casesString <- "COPD Subjects"
  controlsString <- "Smoker Controls"
  phenotypeName <- "Subject.type"
} else if (study=="LTCOPD"){
  exprFile <- "~/gd/Harvard/Research/data/LTCOPD/LTCOPD_exp.txt"
  clinicalFile <- "~/gd/Harvard/Research/data/LTCOPD/LTCOPD_clinical.txt"
  #LTCOPD labels
  casesString <- "COPD"
  controlsString <- "Control"
  phenotypeName <- "diagnosis"
} else if (study=="GTEX"){
  exprFile <- "~/gd/Harvard/Research/data/GTEx/GTEx_expr.txt"
  clinicalFile <- "~/gd/Harvard/Research/data/GTEx/GTEx_clinical.txt"
  #GTEX labels
  casesString <- "cells_ebv-transformed_lymphocytes"
  controlsString <- "skin"
  phenotypeName <- "our_subtypes"
}

# motifFile <- "~/gd/Harvard/Research/data/motifs695.txt"
# exprFile <- "~/gd/Harvard/Research/data/Ovarian/CombinedOV.txt"
# exprFile <- "~/gd/Harvard/Research/data/GTEx/GTEx_expr.txt"
# 
#
#
#
ppiFile <- "~/gd/Harvard/Research/data/Eclipse/OV_PPI.txt"
#clinicalFile <- "~/gd/Harvard/Research/data/Ovarian/Clinical.txt"
# clinicalFile <- "~/gd/Harvard/Research/data/GTEx/GTEx_clinical.txt"

#
#
#clinicalFile <- "~/gd/Harvard/Research/data/LGRC/lgrc.merged.clinical.data.clean.txt"

# exprFile <- "~/gd/Harvard/Research/data/GTEx/gtex_cell.rdata"
# exprFile <- "~/gd/Harvard/Research/data/GTEx/gtex_sub_noxymt_qsmooth_cell.rdata"
# motifFile <- "~/gd/Harvard/Research/data/GTEx/KG_cisbp_652.txt"
# phenotypeName <- "our_subtypes"
# casesString <- "cells_ebv-transformed_lymphocytes"
# controlsString <- "skin"
# filterType <- NA




# 
# casesString <- "Lung"
# controlsString <- "Colon"
# phenotypeName <- "SMTS"

nullPerms <- 500
networkInferenceName <- "bere"
# filterType <- "Gender"
# filterBy <- "M"
filterType <- NA
permuteGeneLabels <- F
analysisName <- paste(study,motifVersion,networkInferenceName,sep="_")

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
    numMaxCores <- as.numeric(args[14])
  } else {
    permuteGeneLabels <- F
    numMaxCores <- 40
  }
  # Create new dir if calling from command line
  outputDir <- file.path("~",paste(analysisName, analysisCode ,sep="_"))
  dir.create(outputDir, showWarnings=FALSE)
  write.csv(args,file=file.path(outputDir,paste("arguments_",analysisCode,".txt",sep="")))
}

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
    # remove genes that do not appear in exp dataset, excluding header
    motifs <- motifs[motifs[,2] %in% rownames(exp[-1,]),]
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
#     system(paste0('rm ',file.path("~/panda_matlab/", paste0("RunPANDA",randindex,".m"))))
    
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
  
  symbols <- mapIds(org.Hs.eg.db, keys=row.names(dataset$exp),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  dataset$exp <- dataset$exp[!is.na(symbols) & !duplicated(symbols),]    
  rownames(dataset$exp) <- symbols[!is.na(symbols) & !duplicated(symbols)]
  
}
dataset$ppi      <- read.table(ppiFile,header=F)
dataset$exp      <- dataset$exp[,order(colnames(dataset$exp))]  # Make sure expression and clinical is in same order

dataset$motif <- dataset$motif[!duplicated(dataset$motif),]

matches <- sort(unique(intersect(rownames(dataset$clinical),colnames(dataset$exp))))
dataset$clinical <- dataset$clinical[matches,]    # Make sure clinical only contains patients with expression data
dataset$exp <- dataset$exp[,matches]    # Make sure expression only contains patients with clinical data

if(permuteGeneLabels){
  print("Permuting gene labels once")
  rownames(dataset$exp) <- sample(rownames(dataset$exp))
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
phenoFilter <- phenoFilter & !is.na(phenoFilter)
allFilter <- subsetFilter&phenoFilter

dataset$exp <- dataset$exp[,allFilter]
dataset$clinical <- dataset$clinical[allFilter,]
casesFilter <- dataset$clinical[,phenotypeName]==casesString
controlsFilter <- dataset$clinical[,phenotypeName]==controlsString

save.image(file=file.path(outputDir,paste0(analysisName,".RData")))
