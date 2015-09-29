require(network)
require(igraph)
require(sna)
require(ggplot2)
require(sna)
require(ergm)
setwd("C:/Users/Dan/Dropbox/Dan's/Documents/Harvard/Research/panda/Version1")
final.net <- read.table("ToyOutput_FinalNetwork.pairs")
head(final.net)
TFs <- levels(factor(final.net[,1]))
genes <- levels(factor(final.net[,2]))
num.TFs <- nlevels(factor(final.net[,1]))
num.Genes <- nlevels(factor(final.net[,2]))
motif.predictions <- final.net[,3]
panda.weights <- final.net[,4]
dist.matrix <- matrix(panda.weights,length(panda.weights)/num.Genes,length(panda.weights)/num.TFs)


#dist.matrix <- pnorm(dist.matrix)
dist.matrix <- data.frame(dist.matrix)
rownames(dist.matrix) <- TFs
colnames(dist.matrix) <- genes
dist.matrix <- removeUnusedVertices(dist.matrix,cutoff=6)
dim(dist.matrix)
setwd("C:/Users/Dan/Dropbox/Dan's/Documents/Harvard/Research/bipartite_plots-master")
#
# The Nava de las Correhuelas dataset.
#
nch <- read.table("./data/NCH_quant_bmatrix.txt", 
                  header = T, sep = "\t", row.names = 1, 
                  dec = ",", na.strings = "NA")
#
# The Hato Raton dataset.
#
hr <- read.table("./data/HR_quant_bmatrix.txt", 
                 header = T, sep = "\t", row.names = 1, 
                 dec = ".", na.strings = "NA")
#
# Sourcing required functions and initializing the net objects.
#
source("./functions/bip_binplot.R")
source("./functions/bip_gplot.R")
source("./functions/bip_qtplot.R")
source("./functions/vectorize.R")
source("./functions/bip_init_network.R")
source("./functions/bip_init_igraph.R")
net <- bip_init_network(nch)
mymat <- nch

net <- bip_init_network(dist.matrix)
mymat <- dist.matrix

net <- bip_init_network(hr)
mymat <- nch

# Now we source the bip_ggplot2.R file
source("./bip_ggplot2.R")





## Combine existing gene expression data into single txt file
setwd("C:/Users/Dan/Dropbox/Dan's/Documents/Harvard/Research/")

f.expr <- read.table("panda/PANDA_COPD_Files/InputData/ECLSputumData_FExpression.txt",row.names=1)
m.expr <- read.table("panda/PANDA_COPD_Files/InputData/ECLSputumData_MExpression.txt",row.names=1)
write.table(cbind(f.expr,m.expr), file="data/ECLSputumData_AllExpression.txt", quote = FALSE, sep="\t", col.names = FALSE)
num.female <- ncol(f.expr)
num.male <- ncol(m.expr)
original.weights <- c(rep(1,num.female),rep(0,num.male))
write.table(original.weights, file="data/weights/weights_F.txt", sep="\n", row.names = FALSE, col.names = FALSE)
write.table(1-original.weights, file="data/weights/weights_M.txt", sep="\n", row.names = FALSE, col.names = FALSE)
system("panda/Version1/PANDA -e data/ECLSputumData_AllExpression.txt -m data/ECLSputumData_PromoterMotif.txt -p data/ECLSputumData_PPI.txt -w data/weights/weights_F.txt -o data/output/copd_outF_fromR", wait = TRUE)
Sys.time()
system("panda/Version1/PANDA -e data/ECLSputumData_AllExpression.txt -m data/ECLSputumData_PromoterMotif.txt -p data/ECLSputumData_PPI.txt -w data/weights/weights_M.txt -o data/output/copd_outM", wait = TRUE)

final.net.1 <- read.table("data/output/copd_outF_FinalNetwork.pairs")
final.net.2 <- read.table("data/output/copd_outM_FinalNetwork.pairs")
final.net.diff <- final.net.1
final.net.diff[,4] <- final.net.1[,4] - final.net.2[,4] 
final.net.diff.observed <- final.net.diff


for (i in 9:11){
  panda.command.1 <- "panda/Version1/PANDA -e data/ECLSputumData_AllExpression.txt -m data/ECLSputumData_PromoterMotif.txt -p data/ECLSputumData_PPI.txt -w"
  panda.weightsfile1 <- paste("data/weights/weights_iteration_",i,"_net1.txt",sep="")
  panda.weightsfile2 <- paste("data/weights/weights_iteration_",i,"_net2.txt",sep="")
  panda.outfile1 <- paste("-o data/output/COPD_out_iteration_",i,"_net1.txt",sep="")
  panda.outfile2 <- paste("-o data/output/COPD_out_iteration_",i,"_net2.txt",sep="")
  perm.weights <- original.weights[sample(length(original.weights), length(original.weights))]
  write.table(perm.weights, file=panda.weightsfile1, sep="\n", row.names = FALSE, col.names = FALSE)
  write.table(1-perm.weights, file=panda.weightsfile2, sep="\n", row.names = FALSE, col.names = FALSE)
  system(paste(panda.command.1,panda.weightsfile1,panda.outfile1, sep=" "), wait = TRUE)
  system(paste(panda.command.1,panda.weightsfile2,panda.outfile2, sep=" "), wait = TRUE)

}

files <- list.files(path="data/output/net1",pattern=".pairs")
perm.edges.1 <- NULL
for (f in files) {
  dat <- read.table(file.path("data/output/net1",f))
  perm.edges.1 <- cbind(perm.edges.1, dat[,4])
}

files <- list.files(path="data/output/net2",pattern=".pairs")
perm.edges.2 <- NULL
for (f in files) {
  dat <- read.table(file.path("data/output/net2",f))
  perm.edges.2 <- cbind(perm.edges.2, dat[,4])
}

perm.edges.combined <- perm.edges.1 - perm.edges.2

# permutation stats
apply(perm.edges.combined,2, function(x){sum(x^2)})

# Observed stats
sum(final.net.diff.observed[,4]^2)

plot(perm.edges.combined[,1], cex=1, pch='.', col="blue", ylim=c(-2.5, 2.5))
plot(perm.edges.combined[,2], cex=1, pch='.', col="blue", ylim=c(-2.5, 2.5))
plot(perm.edges.combined[,3], cex=1, pch='.', col="blue", ylim=c(-2.5, 2.5))
plot(perm.edges.combined[,4], cex=1, pch='.', col="blue", ylim=c(-2.5, 2.5))
plot(final.net.diff.observed[,4], cex=1, pch='.', col="blue", ylim=c(-2.5, 2.5))
points(final.net.diff.observed[,4], cex=1, pch='.', col="red")



### Toy setup
setwd("C:/Users/Dan/Dropbox/Dan's/Documents/Harvard/Research/panda/Version1")
weights <- c(rep(1,25),rep(0,25))
write.table(c(rep(1,25),rep(0,25)), file="weights/weights1.txt", sep="\n", row.names = FALSE, col.names = FALSE)
write.table(c(rep(0,25),rep(1,25)), file="weights/weights2.txt", sep="\n", row.names = FALSE, col.names = FALSE)
system("PANDA -e ToyData/ToyExpressionData.txt -m ToyData/ToyMotifData.txt -p ToyData/ToyPPIData.txt -w weights/weights1.txt -o output/ToyOutput_R1", wait = TRUE)
system("PANDA -e ToyData/ToyExpressionData.txt -m ToyData/ToyMotifData.txt -p ToyData/ToyPPIData.txt -w weights/weights2.txt -o output/ToyOutput_R2", wait = TRUE)

final.net.1 <- read.table("output/ToyOutput_R1_FinalNetwork.pairs")
final.net.2 <- read.table("output/ToyOutput_R2_FinalNetwork.pairs")
final.net.diff <- final.net.1
final.net.diff[,4] <- final.net.1[,4] - final.net.2[,4] 

points(final.net.diff[,4], cex=1, pch='.', col="black")
plot(final.net.1[,4], cex=1, pch='.', col="blue")
points(final.net.2[,4], cex=1, pch='.', col="red")

head(final.net.diff)
head(final.net.1)
head(final.net.2)