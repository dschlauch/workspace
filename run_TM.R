
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


####### Moved this code into the main parallelization 1/13/16
# if(!is.na(num_cores)){
#     cl <- makeCluster(4)
#     registerDoParallel(cl)
# }
# 
# strt  <- Sys.time()
# #loop
# print("Running transition calculations in parallel")
# print(paste0(num_cores," cores used"))
# print(paste0(length(null.networks)/2," transitions to be estimated"))
# transMatrices <- foreach(i=1:(length(null.networks)/2),.packages=c("bptools","reshape2","penalized")) %dopar% {
#     transformation.matrix(null.networks[[2*i]], null.networks[[2*i-1]],remove.diagonal=T,method="ols")    
# }
# 
# print(Sys.time()-strt)
# if(!is.na(num_cores)){
#     stopCluster(cl)
# }

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
# rm(null.networks)
gc()

# periodically save workspace
# save.image(file=file.path(outputDir,paste("activeImage",analysisCode,".RData",sep="")))
