library(ggplot2)
library(reshape2)


makeHeatmap <-  function(){
  expr <- matrix(sapply(rnorm(20), function(x){rnorm(100,mean=x)}),nrow=20)
  heatmap1 <- melt(expr)
  ggHeat <- ggplot(heatmap1, aes(Var1,Var2)) + geom_tile(aes(fill=value)) +
    scale_fill_gradientn(colours=c("blue","black","Yellow"))+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  ggHeat
}
makeTMHeatmap <-  function(){
  expr <- diag(20) +rnorm(400)/6
  expr[2,15]<-1
  expr[15,2]<-1
  expr[14,18]<-1
  expr[18,14]<-1
  expr[2,5]<-1
  expr <- expr[20:1,]
  heatmap1 <- melt(expr)
  ggHeat <- ggplot(heatmap1, aes(Var1,Var2)) + geom_tile(aes(fill=value)) +
    scale_fill_gradientn(colours=c("blue","black","Yellow"))+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  ggHeat
}
png("~/gd/Harvard/Research/R_workspace/TM_manuscript/figures/diagramHeatmap1.png", height=800, width=400)
print(makeHeatmap())
dev.off()
png("~/gd/Harvard/Research/R_workspace/TM_manuscript/figures/diagramHeatmap2.png", height=800, width=400)
print(makeHeatmap())
dev.off()
png("~/gd/Harvard/Research/R_workspace/TM_manuscript/figures/TM_Heatmap.png", height=400, width=400)
print(makeTMHeatmap())
dev.off()

