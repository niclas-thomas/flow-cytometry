rm(list=ls())

library(ggplot2)
library(openCyto)
library(flowPeaks)

source("C:\\Users\\Laura\\Desktop\\Alice\\code\\functions.R")

setwd("C:\\Users\\Laura\\Desktop\\Alice\\gatingSets")
folderFiles <- list.files()

myMarkers <- c("CD24","CD38","CD27","IgM","IgD")

## LOOP OVER GATING SETS IN FOLDER
##################################

trainData <- c()
for (k in c(4,5,8,10,15)){
  
  gs <- load_gs(folderFiles[k])
  metadata <- getData(gs)[[1]]
  gsName <- description(metadata)$GUID
  
  ## DEFINE WHAT SUBSETS ARE DESIRED FOR ANALYSIS
  getNodes( gs, order="bfs")
  subsets <- getNodes( gs, order="bfs", path="auto")[c(9:16)]

  oneSampleData <- get.underlying.data( subsets, gs, metadata, myMarkers)
  trainData <- rbind( trainData, oneSampleData )
  
}
## FLOW PEAK CLUSTERS
#####################
fp <- flowPeaks( trainData[,-c(1)], tol=0)
fp

results <- c()
testStart = 1
testEnd = 15
numClusters = length(fp$peaks$cid)

for (k in c(testStart:testEnd)){
  
  gs <- load_gs(folderFiles[k])
  metadata <- getData(gs)[[1]]
  gsName <- description(metadata)$GUID
  
  ## DEFINE WHAT SUBSETS ARE DESIRED FOR ANALYSIS
  getNodes( gs, order="bfs")
  subsets <- getNodes( gs, order="bfs", path="auto")[c(10,11,15,16)]
  
  testData <- get.underlying.data( subsets, gs, metadata, myMarkers)
  clusterLabels <- assign.flowPeaks( fp, testData[,-c(1)], tol=0, fc=0)
  
  summaryClusters <- c()
  for (x in c(1:numClusters)){
    summaryClusters <- c( summaryClusters, sum(clusterLabels==x) )
  }
  results <- rbind( results, summaryClusters/sum(summaryClusters) )
  
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

g <- ggplot(labeledData, aes(x=class, y=c)) +
  geom_dotplot(dotsize=0.75, binaxis='y', stackdir='center', stackratio=1.5)  +
  xlab('') +
  ylab('Cell Cluster %') +
  stat_summary(fun.data=data_summary, color="red", size=1.25) +
  theme_bw()
g

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
