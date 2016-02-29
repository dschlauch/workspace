
library(nettools)
library(bereR)
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
runNImethods("~/NI_only_0001/readyToGoECLIPSE0001.RData", "ECLIPSE", c("CLR","ARACNE"))
runNImethods("~/NI_only_0001/readyToGoCOPDGene0001.RData", "COPDGene", c("CLR","ARACNE"))
