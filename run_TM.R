# Copy expression data for null network generation
null.exp <- dataset$exp

# Parallelize

# Calculate the number of cores
num_cores <- detectCores() - 4
num_cores <- min(num_cores, numMaxCores)

# Initiate cluster
if(!is.na(num_cores)){
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
}

#start time
strt  <- Sys.time()
iters <- nullPerms+1 # Two networks for each partition, plus observed partition
#loop
print("Running null permutations in parallel")
print(paste0(num_cores," cores used"))
print(paste0(iters," network transitions to be estimated"))
dir.create(file.path(outputDir,"tms"))

# Changed to run two Networks and calculate transition on each iteration 1/13/16
transMatrices <- foreach(i=1:iters,.packages=c("bereR","pandaR","reshape2","penalized","bptools")) %dopar% {
  print(paste0("Running iteration ", i))
  if(i==1){
    null.exp <- dataset$exp
  } else {
    null.exp <- dataset$exp[,sample(1:ncol(dataset$exp))]
    rownames(null.exp) <- rownames(null.exp)[sample(1:nrow(null.exp))]
  }
  null.exp.cases <- null.exp[,casesFilter]
  null.exp.controls <- null.exp[,controlsFilter]
  # Some QC for sparse data
  if (sum(rowSums(null.exp)==0)>0){
    zeroGenes <- which(rowSums(null.exp)==0)
    for(gene in zeroGenes){
      null.exp[gene,] <- rnorm(ncol(null.exp))
    }
  }
  tmpNetCases <- networkInferenceMethod(dataset$motif, null.exp.cases)
  tmpNetControls <- networkInferenceMethod(dataset$motif, null.exp.controls)
  transition.matrix <- transformation.matrix(tmpNetControls, tmpNetCases, remove.diagonal=T,method="ols")    
  print(paste0("Finished running iteration", i))
  saveRDS(transition.matrix,file.path(outputDir,'tms',paste0('tm_',i,'.rds')))
  transition.matrix
}

print(Sys.time()-strt)
if(!is.na(num_cores)){
  stopCluster(cl)
}

gc()

