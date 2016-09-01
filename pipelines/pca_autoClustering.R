rm(list=ls())

library(ggplot2)
library(openCyto)
library(flowMeans)

source("C:\\Users\\Laura\\Desktop\\Alice\\code\\functions.R")

setwd("C:\\Users\\Laura\\Desktop\\Alice\\gatingSets")

gs_list <- lapply(list.files(),function(this_folder){
  load_gs(this_folder)
})

## DEFINE WHAT SUBSETS ARE DESIRED FOR ANALYSIS
getNodes(gs_list[[1]],order="bfs")
subsets <- getNodes(gs_list[[1]],order="bfs",path="auto")[c(10,11,14,15,16)]
myMarkers <- c("CD19","CD20","CD27","IgM")

## LOOP OVER ALL GATING SETS IN FOLDER
######################################

for (k in c(1:length(gs_list))){
  
  gs <- gs_list[[k]]
  metadata <- getData(gs)[[1]]
  gsName <- description(metadata)$GUID
  print(paste("Analysing fcs file: ",gsName,sep=""))
  
  myData <- get.underlying.data( subsets, gs, metadata, myMarkers)
  
  if (k == 1){
    
    ## PCA ON SUBSET POPULATION
    ###########################
    pcaModel <- getPCAModel(myData)
    pred <- pcaPredict(myData, myData[,c(1)], pcaModel)

    ## 2D GGPLOT OF PC1 AND PC2
    ###########################
    plotPCA( pred, subsets, size=2, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\PCA_",gsName,sep=""))
    
    ## FIND CELL CLUSTERS AND PLOT
    ##############################
    clusters <- getKMeans( myData, myMarkers, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\kMeans_",gsName,sep=""))
    
    ## PCA WITH AUTOCLUSTER LABELS
    ##############################
    autoClusterPred <- pcaPredict(myData, clusters@Label, pcaModel)
    plotPCA( autoClusterPred, subsets, size=2, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\PCAclusters_",gsName,sep=""))
    
    ## PLOT CLUSTER MEMBERSHIP
    ##########################
    #plotClusterMembership(clusters, myData, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\clusterMembership_",gsName,sep="") )
    
  }
  else {
    
    ## PCA ON SUBSET POPULATION
    ###########################
    pred <- pcaPredict(myData, myData[,c(1)], pcaModel)

    ## 2D GGPLOT OF PC1 AND PC2
    ###########################
    plotPCA( pred, subsets, size=2, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\PCA_",gsName,sep=""))
    
    ## FIND CELL CLUSTERS AND PLOT
    ##############################
    clusters <- getKMeans( myData, myMarkers, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\kMeans_",gsName,sep=""))
    
    ## PCA WITH AUTOCLUSTER LABELS
    ##############################
    autoClusterPred <- pcaPredict(myData, clusters@Label, pcaModel)
    plotPCA( autoClusterPred, subsets, size=2, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\PCAclusters_",gsName,sep=""))
    
    ## PLOT CLUSTER MEMBERSHIP
    ##########################
    #plotClusterMembership(clusters, myData, paste("C:\\Users\\Laura\\Desktop\\Alice\\figures\\clusterMembership_",gsName,sep="") )
      
  }

}
