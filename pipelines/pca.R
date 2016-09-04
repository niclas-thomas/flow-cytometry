rm(list=ls())

library(ggplot2)
myDataDir <- "C:\\Users\\Laura\\Desktop\\bCellData\\"
myGitDir <- "C:\\Users\\Laura\\Desktop\\git\\flowCytometry\\"
dir.create(file.path(myDataDir,"figures"))
source(paste(myGitDir,"\\src\\functions.R",sep=""))

myMarkers <- c("CD19","CD20","CD27","IgM","IgD","CD24","CD38")
folderFiles <- list.files(paste(myDataDir,"clusteredDataFrames\\",sep=""))

for (k in c(1:length(folderFiles))){
  df <- read.csv(file=paste(myDataDir,"\\clusteredDataFrames\\",folderFiles[k],sep=""))
  pca <- getPCAModel(df[,myMarkers])
  
  predGated <- pcaPredict(df, df[,c("subset")], pca)
  plotPCA(predGated, paste(myDataDir,"figures\\",folderFiles[k],"_pcaGates",sep=""))
  
  predClustered <- pcaPredict(df, df[,c("clusterLabels")], pca)
  plotPCA(predClustered, paste(myDataDir,"figures\\",folderFiles[k],"_pcaClusters",sep=""))
}