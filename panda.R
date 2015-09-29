panda <- function(expr, motif, ppi, alpha=.2, output='output', metric='correlation', maxiterations=100){
  # Get all the proteins that are in PPI and Motif data
  all.proteins <- unique(c(as.vector(ppi[,1]),as.vector(ppi[,2]),as.vector(motif[,1])))
  all.genes    <- rownames(expr)
  if (sum(!(motif[,2] %in% all.genes))>0)
    print('Some genes in motif data do not have expression values')
  
  #Make an interaction table that is symmetric
  ppi.table <- table(c(as.vector(ppi[,1]),all.proteins), c(as.vector(ppi[,2]),all.proteins))
  ppi.table <- ppi.table - diag(nrow(ppi.table))
  ppi.table <- ppi.table + t(ppi.table)
  ppi.table <- ifelse(ppi.table>1,1,ppi.table)
  

  # remove the proteins that have no motif data or only interact with self
  # 7/21/14 Removed this step, we now want all proteins...
#   keep.protein <- all.proteins %in% as.vector(motif[,1]) & (apply(ppi.table,1,sum)-diag(ppi.table)>0)
#   ppi.table <- ppi.table[keep.protein,keep.protein]
  
  #remove interactions with self
  ppi.table <- ppi.table-diag(diag(ppi.table))
  
  #remove motif data for proteins that we don't use
  keep.motif   <- motif[motif[,2] %in% all.genes,]
  
  tf.table <- as.matrix(table(c(as.vector(keep.motif[,1])), c(as.vector(keep.motif[,2]))))
  zero.col.names <- all.genes[!all.genes %in% motif[,2]]
  zero.columns   <- matrix(0,nrow=nrow(tf.table), ncol=length(zero.col.names))
  colnames(zero.columns) <- zero.col.names
  tf.table <- cbind(tf.table,zero.columns)
  
  #get the genes that are missing from motif data
  # Changed 7/21 to stop excluding genes that don't have starting edges.
  expr.subset <- expr#expr[row.names(expr) %in% colnames(tf.table),]
  all.genes <- row.names(expr.subset)
  expr.cor.table   <- cor(t(expr.subset), method='pearson')
  
  
  ppi.table         <- initial.normalize(ppi.table)
  tf.table          <- initial.normalize(tf.table)
  expr.cor.table    <- initial.normalize(expr.cor.table)
  
  tf.table <- tf.table[,order(colnames(tf.table))]
  expr.cor.table <- expr.cor.table[order(rownames(expr.cor.table)),order(colnames(expr.cor.table))]
  
  
  #Print analysis details
  print(paste("Num regulators:", nrow(ppi.table), "Num genes:", nrow(expr.cor.table),"Num conditions:",ncol(expr), "num interactions:",nrow(motif)))


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
  hamming <- Inf
  proteins.out <- as.vector(sapply(all.proteins, function (x) rep(x,length(all.genes))))
  genes.out    <- rep(all.genes,length(all.proteins))
  weights.iterations <- matrix(NA,length(proteins.out),maxiterations)
  i<-1
  for (i in 1:maxiterations){
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
      apply(expr.cor.table, 2, function(expr.row, tf.row){
        similarity.metric(expr.row,tf.row,use="pairwise.complete.obs")
      }, tf.row)
    })
    responsibility<- initial.normalize(responsibility)
    availability<- initial.normalize(availability)
    
    #Update weights
    new.weights <- (1-alpha)*tf.table + t(alpha*(responsibility + availability)/2)
    hamming = mean(abs(new.weights-tf.table), na.rm=TRUE)
    tf.table <- new.weights
    
    tf.table          <- initial.normalize(tf.table)
    # Next up, update the co-regulatory and PPI networks. 7/22/14
    # To clarify: Use updated weights?  Or weights from previous iteration?
    # From diagram, appears that we should use updated weights here....
    
    coregulation <- t(apply(tf.table,2,function(tf.column1){
      apply(tf.table, 2, function(tf.column2, tf.column1){
        similarity.metric(tf.column1,tf.column2,use="pairwise.complete.obs")
      }, tf.column1)
    }))
    cooperativity <- t(apply(tf.table,1,function(tf.row1){
      apply(tf.table, 1, function(tf.row2, tf.row1){
        similarity.metric(tf.row2,tf.row1,use="pairwise.complete.obs")
      }, tf.row1)
    }))
    
    new.coregulation <- (1-alpha)*expr.cor.table + alpha*coregulation
    expr.cor.table   <- new.coregulation
    new.cooperation  <- (1-alpha)*ppi.table + alpha*cooperativity
    ppi.table        <- new.cooperation
    
    #re-standardize (this is not in the PANDA supplement)
    ppi.table         <- initial.normalize(ppi.table)
    expr.cor.table    <- initial.normalize(expr.cor.table)
    
    print(hamming)
    weights.iterations[,i] <- cbind(as.vector(tf.table))
    
  }
  write.table(cbind(proteins.out,genes.out,as.vector(tf.table)), file=paste('output',Sys.Date(),".txt"), sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
}


tanimoto <- function(x, y, use="pairwise.complete.obs"){
  complete.obs <- as.numeric(!is.na(x))*as.numeric(!is.na(y))==1
  x.c <- x[complete.obs]
  y.c <- y[complete.obs]
  (x.c%*%y.c)/(sqrt(x.c%*%x.c + y.c%*%y.c - abs(x.c%*%y.c)))
}

initial.normalize <-function(mat){
  col.means <- apply(mat, 2, mean, na.rm=TRUE)
  col.sds <- apply(mat, 2, sd, na.rm=TRUE)
  for(i in 1:nrow(mat)){
    row.mean <- mean(mat[i,], na.rm=TRUE)
    row.sd <- sd(mat[i,], na.rm=TRUE)
    if(is.na(row.sd)){
      mat[i,] <- 0
    }
    mat[i,] <- (1/sqrt(2))*(mat[i,]-row.mean)/(row.sd)+(1/sqrt(2))*(mat[i,]-col.means)/(col.sds)
  }
  mat[is.na(mat)]<-0
  return(mat)
}

setwd("C:/Users/Dan/Dropbox/Dan's/Documents/Harvard/Research/")

# expressionA.file <- "data/OvarianNetworkFiles/InputFiles/OV_AngiogenicExpression.txt"
# motif.file       <- "data/OvarianNetworkFiles/InputFiles/OV_TFGene.txt"
# ppi.file         <- "data/OvarianNetworkFiles/InputFiles/OV_PPI.txt"
# 
# expressionA.file <- "panda/Version1/YeastData/YeastData_CCExpression.txt"
# motif.file       <- "panda/Version1/YeastData/YeastData_Motif.txt"
# ppi.file         <- "panda/Version1/YeastData/YeastData_PPI.txt"
#
expressionA.file <- "panda/Version1/ToyData/ToyExpressionData.txt"
motif.file       <- "panda/Version1/ToyData/ToyMotifData.txt"
ppi.file         <- "panda/Version1/ToyData/ToyPPIData.txt"



expr  <- read.table(expressionA.file,row.names=1)
motif <- read.table(motif.file)
ppi   <- read.table(ppi.file)

alpha         <- .25
output        <- 'toy.output.panda'
metric        <- 'tanimoto'
maxiterations <- 8
system.time(panda(expr, motif, ppi, alpha, output, metric, maxiterations))



pandaC.result <- read.table("panda/Version1/YeastOutput_FinalNetwork.pairs")
merged.results <- merge(data.frame(weights.iterations),pandaC.result,by=c(1,2))

dim(merged.results)
dim(pandaC.result)