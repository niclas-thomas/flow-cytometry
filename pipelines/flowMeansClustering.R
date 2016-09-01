rm(list=ls())

library(ggplot2)
library(openCyto)
library(flowMeans)

source("C:\\Users\\Laura\\Desktop\\Alice\\code\\functions.R")

setwd("C:\\Users\\Laura\\Desktop\\Alice\\gatingSets")
folderFiles <- list.files()

myMarkers <- c("CD19","CD20","CD27","IgM","IgD","CD24","CD38")

## LOOP OVER GATING SETS IN FOLDER
##################################

trainData <- c()
for (k in c(1)){
  
  gs <- load_gs(folderFiles[k])
  metadata <- getData(gs)[[1]]
  gsName <- description(metadata)$GUID
  
  ## DEFINE WHAT SUBSETS ARE DESIRED FOR ANALYSIS
  getNodes( gs, order="bfs")
  subsets <- getNodes( gs, order="bfs", path="auto")[c(8:16)]

  oneSampleData <- get.underlying.data( subsets, gs, metadata, myMarkers)
  trainData <- rbind( trainData, oneSampleData )
  
  ## PCA ON SUBSET POPULATION
  ###########################
  pcaModel <- getPCAModel(trainData)
  pred <- pcaPredict(trainData, trainData[,c(1)], pcaModel)
  ## 2D GGPLOT OF PC1 AND PC2
  ###########################
  plotPCA( pred, subsets, size=2, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\PCA_",gsName,sep=""))
  ## PCA WITH AUTOCLUSTER LABELS
  ##############################
  #autoClusterPred <- pcaPredict(trainData, clusterLabels, pcaModel)
  #plotPCA( autoClusterPred, subsets, size=2, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\PCAclusters_",gsName,sep=""))
  
  
}

## FLOW MEANS CLUSTERS
######################
fm <- getKMeans( trainData[,-c(1)], myMarkers, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\kMeans_",gsName,sep=""))
fm

results <- c()
testStart = 1
testEnd = 15
numClusters = length(table(fm@Label))

centres <- getClusterCentres(trainData,numClusters)
matrixCentres <- as.matrix(centres)

for (k in c(testStart:testEnd)){
  
  gs <- load_gs(folderFiles[k])
  metadata <- getData(gs)[[1]]
  gsName <- description(metadata)$GUID
  
  ## DEFINE WHAT SUBSETS ARE DESIRED FOR ANALYSIS
  getNodes( gs, order="bfs")
  subsets <- getNodes( gs, order="bfs", path="auto")[c(8:16)]
  
  testData <- get.underlying.data( subsets, gs, metadata, myMarkers)
  matrixTestData <- as.matrix(testData[,-c(1)])
  
  clusterLabels <- unlist(lapply(seq_len(nrow(matrixTestData)), function(i) which.min(sqrt(colSums((matrixTestData[i, ] - t(matrixCentres))^2)))))

  summaryClusters <- c()
  for (x in c(1:numClusters)){
    summaryClusters <- c( summaryClusters, sum(clusterLabels==x) )
  }
  results <- rbind( results, summaryClusters/sum(summaryClusters) )

  ## PCA ON SUBSET POPULATION
  ###########################
  pred <- pcaPredict(testData, testData[,c(1)], pcaModel)
  ## 2D GGPLOT OF PC1 AND PC2
  ###########################
  plotPCA( pred, subsets, size=2, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\PCA_",gsName,sep=""))
  ## PCA WITH AUTOCLUSTER LABELS
  ##############################
  autoClusterPred <- pcaPredict(testData, clusterLabels, pcaModel)
  plotPCA( autoClusterPred, subsets, size=2, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\PCAclusters_",gsName,sep=""))
  
}

cohortData <- as.data.frame(results)
rownames(cohortData) <- list.files()[testStart:testEnd]

## BUILD DATA FRAME WITH CLASS LABELS
#####################################
class <- c("P","HC","HC","P","P","HC","P","HC","P","P","HC","P","HC","HC","HC","HC")[testStart:testEnd]
labeledData <- data.frame( class, cohortData)
colnames(labeledData) <- c("class",letters[1:(dim(labeledData)[2]-1)])

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

for (cluster in colnames(labeledData)[-c(1)]){
  print(wilcox.test( eval(parse(text=cluster)) ~ class, data=labeledData)$p.value)
}

g <- ggplot(labeledData, aes(x=class, y=d)) +
  geom_dotplot(dotsize=0.75, binaxis='y', stackdir='center', stackratio=1.5)  +
  xlab('') +
  ylab('Cell Cluster %') +
  stat_summary(fun.data=data_summary, color="red", size=1.25) +
  theme_bw()
g
#ggsave('C:\\Users\\Laura\\Desktop\\Alice\\comparePopulations.png',g)

# pca <- getPCAModel(labeledData)
# prediction <- predict(pca, labeledData[,-c(1)])
# pred <- data.frame(class, prediction)
# pred <- as.data.frame(pred)
# pred$class <- as.factor(pred$class)
# 
# g <- ggplot(pred, aes(PC1, PC2)) +
#  geom_point(aes(colour=class,size=2))
# g
# 

## LOOK At CLUSTER CHARACTERISTICS
##################################
fullData <- cbind( fm@Label, trainData)
colnames(fullData) <- c("clusters",colnames(fullData)[2:length(fullData)])
fullData$clusters <- factor(fullData$clusters)

g <- ggplot(fullData, aes(x=CD24))+
  geom_density(aes(group=clusters, colour=clusters, fill=clusters), alpha=0.2)+
  theme_bw()
g
ggsave('C:\\Users\\Laura\\Desktop\\Alice\\cd24_densities.png',g)
