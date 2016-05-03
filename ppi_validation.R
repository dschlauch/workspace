library(reshape2)
library(ROCR)

load("~/gd/Harvard/Research/TM_outputs/JASPAR2014/BERE/ECLIPSE_combined_runs/activeImage.RData")

plotList <- lapply(1:401, function(x){
  meltTM <- melt(abs(transMatrices[[x]]))
  colnames(dataset$ppi) <- c("Var1","Var2","ppi")
  ppi <- rbind(dataset$ppi, dataset$ppi[,c(2,1,3)])
  mergedTM <- merge(meltTM, ppi, by=c("Var1","Var2"), all.x=T)
  mergedTM[is.na(mergedTM)] <- 0
  mergedTM <- mergedTM[as.character(mergedTM[,1])!=as.character(mergedTM[,2]),]
  
  methodPred  <- prediction(mergedTM[,3], mergedTM[,4])
  roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
  auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
  list("roc.methodPred"=roc.methodPred, "auc.methodPred"=auc.methodPred)
})


plot(plotList[[1]][["roc.methodPred"]], main="PPI Prediction", col = "red", lwd=3)
mapply(function(x,index){
  lines(plotList[[x]][["roc.methodPred"]]@x.values[[1]], plotList[[x]][["roc.methodPred"]]@y.values[[1]], col = "blue", lwd=.5)
}, 2:length(plotList), 2:length(plotList))
lines(plotList[[1]][["roc.methodPred"]]@x.values[[1]], plotList[[1]][["roc.methodPred"]]@y.values[[1]], col = "red", lwd=3)
p.value <- 1-pnorm(scale(unlist(lapply(plotList, function(x){x$auc.methodPred})))[1])
legend(.5,.2, c("Control to COPD", "Randomized Phenotype"), lty=c(1,1),lwd=c(2,2),col=c("red", "blue"),title="Area under ROC curve")
c("AUC"=plotList[[1]]$auc.methodPred,"p-value="=p.value)


