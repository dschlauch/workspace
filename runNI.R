
library(nettools)
library(bereR)
library(foreach)
library(doParallel)
# ECLIPSE


runNImethods <- function(data="~/NI_only_0001/readyToGoCOPDGene0001.RData", dataName="COPDGene", niNames=c("WGCNA","CLR","ARACNE")){
  load(data)
  exp.cases <- dataset$exp[,casesFilter]
  exp.controls <- dataset$exp[,controlsFilter]
  # Some QC for sparse data
  if (sum(rowSums(dataset$exp)==0)>0){
    zeroGenes <- which(rowSums(dataset$exp)==0)
    for(gene in zeroGenes){
      dataset$exp[gene,] <- rnorm(ncol(dataset$exp))
    }
  }
  
  TFs <- as.character(unique(dataset$motif[,1]))
  # BERE
  if ("BERE"%in%niNames){
    netCases <- bere(dataset$motif, exp.cases, score="no-adjust")
    netControls <- bere(dataset$motif, exp.controls, score="no-adjust")
    saveRDS(netCases,file.path(outputDir,paste0(dataName, '_cases_bere_network.rds')))
    saveRDS(netControls,file.path(outputDir,paste0(dataName, '_controls_bere_network.rds')))
  }
  if("PANDA"%in%niNames){
    # Possibly implement for PANDA, maybe use matlab script
  }
  nettoolsNames <- niNames["BERE"!=niNames]
  
  sapply(nettoolsNames, function(niName){
    print(niName)
    netCases <- mat2adj(t(exp.cases), method=niName)[TFs,]
    saveRDS(netCases,file.path(outputDir,paste0(dataName,'_cases_',niName,'_network.rds')))
    netControls <- mat2adj(t(exp.controls), method=niName)[TFs,]
    saveRDS(netControls,file.path(outputDir,paste0(dataName,'_controls_',niName,'_network.rds')))
  })
}

print("Begin loading data...")
# load the data, save each study as RData
studies <- c("ECLIPSE", "COPDGENE", "LGRC","LTCOPD")
motifVersion <- "JASPAR2014"
outputDir <- "~/NI_only_0001"
sapply(studies, function(study){
  print(study)
  study <<- study
  source("~/gd/Harvard/Research/R_workspace/load_TM_data.R")
})
print("Finished loading data.")

cl <- makeCluster(4)
registerDoParallel(cl)

#start time
strt  <- Sys.time()
foreach(i=studies,.packages=c("bereR","nettools")) %dopar% {
  runNImethods(paste0(outputDir,"/",i,"_JASPAR2014_bere.RData"), i, c("WGCNA","CLR","ARACNE"))
}
print(Sys.time()-strt)
stopCluster(cl)
