## FUNCTIONS
############
get.ab.names <- function(data, metadata){
  ab.names <- c()
  for (item in colnames(data[,-c(1)])){
    ab.names <- c(ab.names, getChannelMarker(metadata, item)$desc)
  }
  return(ab.names)
}

get.underlying.data <- function( subsets, gatingSet, metadata){
  total.data <- c()
  for (i in c(1:length(subsets))){
    flow.data <- lapply(gatingSet,getData,y=subsets[i])
    subset.data <- data.frame()
    for (j in c(1:length(flow.data))){
      x <- flow.data[[j]]
      x <- x[,-grep("FSC",colnames(x))]
      x <- x[,-grep("SSC",colnames(x))]
      subset.data <- rbind(subset.data,exprs(x))
    }
    if (dim(subset.data)[1] > 0){
      subset.data <- cbind(subsets[i],subset.data)
      
      marker.names <- get.ab.names(subset.data, metadata)
      marker.names <- sub("-",".",marker.names)
      marker.names <- sub("/",".",marker.names)
      marker.names <- sub(" ",".",marker.names)
      
      colnames(subset.data) <- c("subset",marker.names)
      total.data <- rbind(total.data,subset.data)
    }
  }
  return(total.data)
  
}

getPCAModel <- function(data){
  
  pca <- prcomp(~., data=data, cor = TRUE, scale=T)
  return(pca)
  
}

pcaPredict <- function(data, labels, pca){
  
  prediction <- predict(pca, data[,-c(1)])
  pred <- data.frame(labels, prediction)
  pred <- as.data.frame(pred)
  colnames(pred) <- c("subset", colnames(pred)[-1])
  pred$subset <- as.factor(pred$subset)
  return(pred)
  
}

plotPCA <- function( pred, pathAndFilename){
  
  g <- ggplot(pred, aes(PC1, PC2)) +
    geom_point(aes(color=subset),size=2)+
    theme_bw()+
    scale_color_discrete(name="Cell Population")+#,
    #                     #breaks=c(1:length(subsets)),
    #                     labels=subsets)+
    guides(colour = guide_legend(override.aes = list(size=4)))
  ggsave(file=paste(pathAndFilename,".png",sep=""), g, scale=1, width=10, height=10, dpi=150)
  
}

getKMeans <- function( data, varNames, pathAndFilename){
  
  res <- flowMeans( data, varNames, Mahalanobis = FALSE )
  print(paste("Number of clusters found: ",length(table(res@Label)),sep=""))
  #png( paste(pathAndFilename,".png",sep=""))
  #plot(data, res, varNames, pch='.', cex=3)
  #dev.off()
  return(res)
  
}

compareClusters <- function( df ){
  result <- c() 
  for (i in unique(df$cluster)){
    for (j in unique(df$subset)){
      count = nrow(df[df$cluster == i & df$subset == j,])
      result <- rbind( result, c(i,j,count))
    }
  }
  result <- as.data.frame(result)
  colnames(result) <- c("cluster", "subset", "count")
  result$cluster <- as.numeric(as.character(result$cluster))
  result$count <- as.numeric(as.character(result$count))
  return(result)
}

plotClusterMembership <- function(clusters, myData, pathAndFilename ){
  
  df <- as.data.frame(cbind( clusters@Label, as.character(myData[,c(1)])))
  colnames(df) <- c("cluster", "subset")
  df$cluster <- as.numeric(as.character(df$cluster))
  compare <- compareClusters( df )
  g <- ggplot( compare, aes(cluster,count,fill=subset))+
    geom_bar(position = "fill", stat = "identity")+
    xlab("Cluster")+
    ylab("Proportion")
  ggsave(file=paste(pathAndFilename,".png",sep=""), g, scale=1, width=10, height=10, dpi=72)
  
}

getClusterCentres <- function( trainData, numClusters ){
  centres <- c()
  for (i in c(1:numClusters)){
    centres <- rbind( centres, colMeans(head(trainData[which(fm@Label==i),-c(1)])))
  }
  return(centres)
}

plotPCAScores <- function(){
  
  #scores <- as.data.frame(pca$rotation)
  #scores <- cbind(ab.names,scores)
  
  ## SHOW PC SCORES
  #################
  # p1 <- ggplot(data=scores,aes(x=ab.names,y=PC1))+
  #   geom_bar(stat="identity",colour="black")+
  #   xlab("")+
  #   ylim(c(-1,1))+
  #   theme(axis.text.x=element_text(angle=90, vjust=0.5, size=13,colour="black",hjust=1),
  #         axis.text.y=element_text(size=13,angle=90,colour="black"))
  # p2 <- ggplot(data=scores,aes(x=ab.names,y=PC2))+
  #   geom_bar(stat="identity",colour="black")+
  #   xlab("")+
  #   ylim(c(-1,1))+
  #   theme(axis.text.x=element_text(angle=90, vjust=0.5, size=13,colour="black",hjust=1),
  #         axis.text.y=element_text(size=13,angle=90,colour="black"))
  # pdf(paste("/Users/niclasthomas/Dropbox/pid/figures/",list.files()[k],"-PCscores",".pdf",sep=""))
  # grid.arrange(p1,p2,ncol=2)
  # dev.off()
  return(0)
}
