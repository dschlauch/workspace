hcl.heatmap.plot <- function(x, method="pearson"){
  library(ggplot2)
  library(reshape2)
  library(ggdendro)
  library(grid) 
  library(hyperSpec)
  if(method=="pearson"){
    dist.func <- pearson.dist
  } else {
    dist.func <- dist
  }
  # x <- as.matrix(scale(mtcars))
  x <- scale(x)
  dd.col <- as.dendrogram(hclust(dist.func(x)))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist.func(t(x))))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- x[col.ord, row.ord]
  xx_names <- attr(xx, "dimnames")
  df <- as.data.frame(xx)
  colnames(df) <- xx_names[[2]]
  df$Var1 <- xx_names[[1]]
  df$Var1 <- with(df, factor(Var1, levels=Var1, ordered=TRUE))
  mdf <- melt(df)
  
  
  ddata_x <- dendro_data(dd.row)
  ddata_y <- dendro_data(dd.col)
  
  ### Set up a blank theme
  theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
    #axis.ticks.length = element_blank()
  )
  ### Set up a blank theme
  theme_heatmap <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
    #axis.ticks.length = element_blank()
  )
  ### Create plot components ###    
  # Heatmap
  p1 <- ggplot(mdf, aes(x=variable, y=Var1)) + 
    geom_tile(aes(fill=value)) + scale_fill_gradient2() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Dendrogram 1
  p2 <- ggplot(segment(ddata_x)) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    theme_none + theme(axis.title.x=element_blank())
  
  # Dendrogram 2
  p3 <- ggplot(segment(ddata_y)) + 
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
    coord_flip() + theme_none
  
  ### Draw graphic ###
  
  grid.newpage()
  print(p1, vp=viewport(0.80, 0.8, x=0.400, y=0.40))
  print(p2, vp=viewport(0.73, 0.2, x=0.395, y=0.90))
  print(p3, vp=viewport(0.20, 0.8, x=0.910, y=0.43))
}