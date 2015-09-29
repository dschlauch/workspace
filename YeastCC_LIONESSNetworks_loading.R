library(tidyr)
library(dplyr)
setwd("C:/Users/Dan/Dropbox/Dan's/Documents/Harvard/Research/YeastCC_LIONESSNetworks")
lioness.nets <- read.table("YeastCCNetworks.txt",header=T)
for(i in 5:52){
  tmp <- lioness.nets[,c(1,2,i)]
  colnames(tmp)[3] <- "edgeweight"
  regnet <- xtabs(edgeweight~TF+Gene, data=tmp)
}

regnet.list <- apply(lioness.nets[,5:8],2,function(x){
  tmp <- cbind(lioness.nets[,c(1,2)],x)
  colnames(tmp)[3] <- "edgeweight"
  regnet <- xtabs(edgeweight~TF+Gene, data=tmp)
  regnet
})

tmp <- lioness.nets[,c(1,2,6)]
colnames(tmp)[3] <- "edgeweight"
net1.spread <- tmp %>% spread(Gene, edgeweight)
tmp <- lioness.nets[,c(1,2,7)]
colnames(tmp)[3] <- "edgeweight"
net2.spread <- tmp %>% spread(Gene, edgeweight)

library(PANDA)
regnet1 <- list(reg.net=net1)
transformation.matrix(net1, net2)
%>% heatmap.2(dendrogram="both",col=bluered)

data(yeast)
stocks <- data.frame(
  time = as.Date('2009-01-01') + 0:9,
  X = rnorm(10, 0, 1),
  Y = rnorm(10, 0, 2),
  Z = rnorm(10, 0, 4)
)
stocksm <- stocks %>% gather(stock, price, -time)
stocksm %>% spread(stock, price)
stocksm %>% spread(time, price)
