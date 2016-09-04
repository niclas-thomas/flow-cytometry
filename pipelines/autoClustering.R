rm(list=ls())

library(flowMeans)

myDataDir <- "C:\\Users\\Laura\\Desktop\\bCellData\\"
myGitDir <- "C:\\Users\\Laura\\Desktop\\git\\flowCytometry\\"
dir.create(file.path(myDataDir,"clusteredDataFrames"))
dir.create(file.path(myDataDir,"figures"))
source(paste(myGitDir,"\\src\\functions.R",sep=""))

folderFiles <- list.files(paste(myDataDir,"gatedDataFrames",sep=""))
myMarkers <- c("CD19","CD20","CD27","IgM","IgD","CD24","CD38")

## LOOP OVER DATAFRAMES IN FOLDER
#################################
results <- c()

for (k in c(1)){

  trainData <- read.csv(file=paste(myDataDir,"gatedDataFrames\\",folderFiles[k],sep=""),header=TRUE)

  ## FLOW MEANS CLUSTERS
  ######################
  fm <- getKMeans( trainData[,myMarkers], myMarkers, paste(myDataDir,"figures\\kMeans_",folderFiles[k],sep=""))
  
  numClusters = length(table(fm@Label))
  
  centres <- getClusterCentres(trainData,numClusters)
  matrixCentres <- as.matrix(centres)
  
  matrixTrainData <- as.matrix(trainData[,-c(1)])
  
  clusterLabels <- unlist(lapply(seq_len(nrow(matrixTrainData)), function(i) which.min(sqrt(colSums((matrixTrainData[i, ] - t(matrixCentres))^2)))))
  
  summaryClusters <- c()
  for (x in c(1:numClusters)){
    summaryClusters <- c( summaryClusters, sum(clusterLabels==x) )
  }
  results <- rbind( results, summaryClusters/sum(summaryClusters) )
  
  appendedData <- cbind(trainData, clusterLabels)
  write.csv(appendedData,
            file = paste(myDataDir,"clusteredDataFrames\\",folderFiles[k],"_clustered",sep=""),
            row.names=FALSE,
            quote = FALSE)
}

testStart = 2
testEnd = 2

for (k in c(testStart:testEnd)){
  
  testData <- read.csv(file=paste(myDataDir,"gatedDataFrames\\",folderFiles[k],sep=""))
  matrixTestData <- as.matrix(testData[,-c(1)])
  
  clusterLabels <- unlist(lapply(seq_len(nrow(matrixTestData)), function(i) which.min(sqrt(colSums((matrixTestData[i, ] - t(matrixCentres))^2)))))

  summaryClusters <- c()
  for (x in c(1:numClusters)){
    summaryClusters <- c( summaryClusters, sum(clusterLabels==x) )
  }
  results <- rbind( results, summaryClusters/sum(summaryClusters) )
  
  appendedData <- cbind(testData, clusterLabels)
  write.csv(appendedData,
            file = paste(myDataDir,"clusteredDataFrames\\",folderFiles[k],"_clustered",sep=""),
            row.names=FALSE,
            quote = FALSE)

}
