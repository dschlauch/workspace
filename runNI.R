
library(nettools)
library(bereR)
# ECLIPSE

# load the data, save each study as RData
studies <- c("ECLIPSE", "COPDGENE", "LGRC","LTCOPD")
motifVersion <- "JASPAR2014"
outputDir <- "~/NI_only_0001"
sapply(studies, function(study){
  study <<- study
  source("~/gd/Harvard/Research/R_workspace/load_TM_data.R")
})

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
  
  TFs <- unique(dataset$motif[,1])
  # BERE
  if ("BERE"%in%niNames){
    netCases <- bere(dataset$motif, exp.cases, score="no-adjust")
    netControls <- bere(dataset$motif, exp.controls, score="no-adjust")
    saveRDS(netCases,file.path(outputDir,paste0(dataName, '_cases_bere_network.rds')))
    saveRDS(netControls,file.path(outputDir,paste0(dataName, '_controls_bere_network.rds')))
  }
  
  nettoolsNames <- niNames["BERE"!=niNames]
  
  sapply(nettoolsNames, function(niName){
    print(niName)
    netCases <- mat2adj(t(exp.cases), method=niName)[TFs,]
    netControls <- mat2adj(t(exp.controls), method=niName)[TFs,]
    saveRDS(netCases,file.path(outputDir,paste0(dataName,'_cases_',niName,'_network.rds')))
    saveRDS(netControls,file.path(outputDir,paste0(dataName,'_controls_',niName,'_network.rds')))
  })
}


num_cores <- 4

# Initiate cluster
if(!is.na(num_cores)){
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
}

#start time
strt  <- Sys.time()
# Changed to run two Networks and calculate transition on each iteration 1/13/16
foreach(i=studies,.packages=c("bereR","nettools")) %dopar% {
  runNImethods(paste0(outputDir,"/",i,"_JASPAR2014_bere.RData"), i, c("BERE","WGCNA","CLR"))
}

print(Sys.time()-strt)
if(!is.na(num_cores)){
  stopCluster(cl)
}
