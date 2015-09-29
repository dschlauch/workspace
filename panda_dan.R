panda <- function(expr, motif, ppi, alpha=.2, output='output', metric='correlation'){
  # Get all the proteins that are in PPI and Motif data
  all.proteins <- unique(c(as.vector(ppi[,1]),as.vector(ppi[,2])))
  
  #Make an interaction table that is symmetric
  ppi.table <- table(c(as.vector(ppi[,1]),all.proteins), c(as.vector(ppi[,2]),all.proteins))
  ppi.table <- ppi.table - diag(nrow(ppi.table))
  ppi.table <- ppi.table + t(ppi.table)
  ppi.table <- ifelse(ppi.table>1,1,ppi.table)
  
  # This step needed to account for non-binary PPI (don't know why this happens)
  for(row in 1:nrow(ppi)){
    ppi.table[as.character(ppi[row,1]),as.character(ppi[row,2])] <- as.numeric(as.character(ppi[row,3]))
    ppi.table[as.character(ppi[row,2]),as.character(ppi[row,1])] <- as.numeric(as.character(ppi[row,3]))
  }
  
  # remove the proteins that have no motif data or only interact with self
  keep.protein <- all.proteins %in% as.vector(motif[,1]) & (apply(ppi.table,1,sum)-diag(ppi.table)>0)
  ppi.table <- ppi.table[keep.protein,keep.protein]
  
  #remove interactions with self
  ppi.table <- ppi.table-diag(diag(ppi.table))
  all.proteins <- row.names(ppi.table)
  
  #remove motif data for proteins that we don't use
  keep.motif   <- motif[motif[,1] %in% all.proteins,]
  tf.table <- as.matrix(table(c(as.vector(keep.motif[,1])), c(as.vector(keep.motif[,2]))))
  
  #get the genes that are missing from motif data
  expr.subset <- expr[row.names(expr) %in% colnames(tf.table),]
  all.genes <- row.names(expr.subset)
  expr.cor.matrix   <- cor(t(expr.subset), method='pearson')
#   
#   ppi.table         <- initial.normalize(ppi.table)
#   tf.table          <- initial.normalize(tf.table)
#   expr.cor.table    <- initial.normalize(expr.cor.matrix)
  
  # 6/25 finished initialization step.  next implement tanimoto updating rules.
  
  # Begin updating iterations
  
  #6/25 Have not tested code below this point (will take awhile to run)
  
  #Specify distance metric
  if (metric=='correlation'){
    similarity.metric <- cor
  } else if (metric=='tanimoto'){
    similarity.metric <- tanimoto
  }
  #Update loop (replace with Hamming convergence stop)
  maxiterations <- 100
  hamming <- Inf
  proteins.out <- as.vector(sapply(all.proteins, function (x) rep(x,length(all.genes))))
  genes.out    <- rep(all.genes,length(all.proteins))
  weights.iterations <- matrix(NA,length(proteins.out),maxiterations)
  for (i in 1:2){
    if (hamming < .001){
      break
    }
    #Calculate responsibility and availability
    responsibility <- t(apply(tf.table,2,function(tf.column){
      apply(ppi.table, 1, function(ppi.row, tf.column){
        similarity.metric(ppi.row,tf.column,use="pairwise.complete.obs")
      }, tf.column)
    }))
    availability <- apply(tf.table,1,function(tf.row){
      apply(expr.cor.matrix, 2, function(expr.row, tf.row){
        similarity.metric(expr.row,tf.row,use="pairwise.complete.obs")
      }, tf.row)
    })
    
    #Update weights
    new.weights <- (1-alpha)*tf.table + t(alpha*(responsibility + availability)/2)
    hamming = mean(abs(new.weights-tf.table), na.rm=TRUE)
    tf.table <- new.weights
    print(hamming)
    weights.iterations[,i] <- cbind(as.vector(tf.table))
    
  }
  write.table(cbind(proteins.out,genes.out,as.vector(tf.table)), file='output', sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

tanimoto <- function(x, y, use="pairwise.complete.obs"){
  complete.obs <- as.numeric(!is.na(x))*as.numeric(!is.na(y))==1
  x.c <- x[complete.obs]
  y.c <- y[complete.obs]
  (x.c%*%y.c)/(sqrt(x.c%*%x.c + y.c%*%y.c - 2*abs(x.c%*%y.c)))
}

initial.normalize <-function(mat){
  col.means <- apply(mat, 2, mean, na.rm=TRUE)
  col.sds <- apply(mat, 2, sd, na.rm=TRUE)
  for(i in 1:nrow(mat)){
    row.mean <- mean(mat[i,], na.rm=TRUE)
    row.sd <- sd(mat[i,], na.rm=TRUE)
    mat[i,] <- (1/sqrt(2))*(mat[i,]-row.mean)/(row.sd)+(1/sqrt(2))*(mat[i,]-col.means)/(col.sds)
  }
  return(mat)
}

setwd("C:/Users/Dan/Dropbox/Dan's/Documents/Harvard/Research/")

# expressionA.file <- "data/OvarianNetworkFiles/InputFiles/OV_AngiogenicExpression.txt"
# motif.file       <- "data/OvarianNetworkFiles/InputFiles/OV_TFGene.txt"
# ppi.file         <- "data/OvarianNetworkFiles/InputFiles/OV_PPI.txt"

expressionA.file <- "panda/Version1/YeastData/YeastData_KOExpression.txt"
motif.file       <- "panda/Version1/YeastData/YeastData_Motif.txt"
ppi.file         <- "panda/Version1/YeastData/YeastData_PPI.txt"


expr  <- read.table(expressionA.file,row.names=1)
motif <- read.table(motif.file)
ppi   <- read.table(ppi.file)

alpha <- .2
output <- 'output.panda'
system.time(panda(expr, motif, ppi, alpha, output))



pandaC.result <- read.table("panda/Version1/YeastOutput_FinalNetwork.pairs")
merged.results <- merge(data.frame(weights.iterations),pandaC.result,by=c(1,2))

dim(merged.results)
dim(pandaC.result)