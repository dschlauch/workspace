nVariants <- 10000000
germline <- rbinom(nVariants,1,.01)
errorRate <- .00002
somaticMutRate <- .00008
errorsPerSample <- errorRate*nVariants
somatPerSample <- somaticMutRate*nVariants

clusterALoci <- sample(nVariants, somatPerSample)
clusterBLoci <- sample(nVariants, somatPerSample)
clusterCLoci <- sample(nVariants, somatPerSample)
clusterDLoci <- sample(nVariants, somatPerSample)
clusterELoci <- sample(nVariants, somatPerSample)
clusterFLoci <- sample(nVariants, somatPerSample)
clusterGLoci <- sample(nVariants, somatPerSample)

clusterA <- germline
clusterA[clusterALoci] <- 1-clusterA[clusterALoci]
clusterB <- clusterA
clusterB[clusterBLoci] <- 1-clusterB[clusterBLoci]
clusterC <- clusterA
clusterC[clusterCLoci] <- 1-clusterC[clusterCLoci]
clusterD <- clusterB
clusterD[clusterDLoci] <- 1-clusterD[clusterDLoci]
clusterE <- clusterB
clusterE[clusterELoci] <- 1-clusterE[clusterELoci]
clusterF <- clusterC
clusterF[clusterFLoci] <- 1-clusterF[clusterFLoci]
clusterG <- clusterF
clusterG[clusterGLoci] <- 1-clusterG[clusterGLoci]


data <- cbind(germline,clusterA,clusterB,clusterC,clusterD,clusterE,clusterF,clusterG)
dataObs <- apply(data,2, function(x){
  errors <- sample(length(x),errorsPerSample)
  x[errors] <- 1-x[errors]
  x
})

rownames(dataObs) <- 1:nVariants

write.table(dataObs,file="./simulatedPhylogeneticData.txt", row.names = T, sep="\t",quote = F)
system('tar -czvf ./simulatedPhylogeneticData.txt.gz ./simulatedPhylogeneticData.txt')


# dm <- dist(t(dataObs))
# plot(hclust(dm))

library(igraph)
adj <- cor2pcor(cor(dataObs))
diag(adj) <- 0
round(adj,3)
adj[adj>.2] <- 1
adj[adj<.2] <- 0
colnames(adj) <- c("Germline",LETTERS[1:7])
rownames(adj) <- c("Germline",LETTERS[1:7])

g <- graph.adjacency(adj, mode="undirected", weighted=T)
plot(g, layout=layout_with_kk)
