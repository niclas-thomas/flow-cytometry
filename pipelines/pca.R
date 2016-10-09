## SET FOLDER LOCATIONS
## AND THE FILE TO ANALYSE
##########################

rm(list=ls())

myDataDir <- "C:\\Users\\Laura\\Desktop\\bCellData\\"
myCodeDir <- "C:\\Users\\Laura\\Desktop\\git\\flowCytometry\\"

## PROJECT CELL POPULATIONS
## IN 2-D SPACE USING PCA
#########################

library(ggplot2)

options(warn=-1)
dir.create(file.path(myDataDir,"figures"))
source(paste(myCodeDir,"\\src\\functions.R",sep=""))

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