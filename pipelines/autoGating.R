## SET FOLDER LOCATIONS
## AND THE FILE TO ANALYSE
##########################

rm(list=ls())

myDataDir <- "C:\\Users\\Laura\\Desktop\\bCellData\\"
myCodeDir <- "C:\\Users\\Laura\\Desktop\\git\\flowCytometry\\"

plot.strategy = FALSE
plot.gates = FALSE

## BEGIN AUTO GATING
####################

library(openCyto)

options(warn=-1)
dir.create(file.path(myDataDir,"gatedDataFrames"))
source(paste(myCodeDir,"\\src\\functions.R",sep=""))

folderFiles <- list.files(myDataDir,pattern=".fcs")

## READ IN T CELL GATING TEMPLATES
##################################
bcell.gating.template <- gatingTemplate(paste(myCodeDir,"gatingTemplates\\","gatingtemplate_bcell.csv",sep=""))

for (i in c(1:length(folderFiles))){
  
  fcsFile <- folderFiles[i]
  print(paste("Processing ",fcsFile,"...",sep=""))
  
  ## READ IN FCS FILES, SPLITTING FILES INTO ISOTYPE CONTROL AND OTHER FILES
  ##########################################################################
  samplesFlowSet <- read.flowSet(paste(myDataDir,fcsFile,sep=""),alter.names=TRUE)
  sampleNames(samplesFlowSet) <- strsplit(fcsFile,".fcs")[[1]]
  
  ## COMPENSATION
  ###############
  apply.compensation <- function(frame){
    colnames(keyword(frame)$`SPILL`) <- gsub("-",".",colnames(keyword(frame)$`SPILL`))
    colnames(keyword(frame)$`SPILL`) <- gsub(" ",".",colnames(keyword(frame)$`SPILL`))
    comp <- keyword(frame)$`SPILL`
    new_frame <- compensate(frame,comp)
    new_frame
  }
  
  samplesFlowSet.comp <- fsApply(samplesFlowSet,apply.compensation)
  
  ## REMOVE VARIABLES FSC, SSC .. FROM LOGICLE TRANSFORM
  ######################################################
  vars <- colnames(samplesFlowSet.comp)
  vars <- vars[-grep("FSC",vars)]
  vars <- vars[-grep("SSC",vars)]
  
  ## TRANSFORM DATA USING LOGICLE TRANSFORM
  #########################################
  lgcl <- estimateLogicle(samplesFlowSet.comp[[1]], channels = vars)
  samplesFlowSet.trans <- transform(samplesFlowSet.comp, lgcl)
  
  ## PERFORM AUTOMATED GATING
  ###########################
  samples.gating.set <- GatingSet(samplesFlowSet.trans)
  gating(bcell.gating.template, samples.gating.set, mc.cores=1, parallel_type = "multicore")
  
  ## GET POPULATION STATS
  #######################
  summary.stats.samples <- getPopStats(samples.gating.set[[1]])
  summary.stats.samples
  
  ## PLOT RESULTS
  ###############
  if (plot.strategy == TRUE){
    jpeg(paste(name[length(name)],'.jpg',sep=""))
    plot(samples.gating.set)
    def.off()
  }
  
  if (plot.gates == TRUE){
    jpeg(paste(name[length(name)],'.jpg',sep=""))
    plotGate(samples.gating.set[[1]],
             sample.ratio=0.5,
             #xbin=100,
             default.y="SSC.A",
             par.settings=list( cex=10,
                                axis.text = list(cex = 1.5),
                                gate.text = list( background = list(fill = "white"),cex = 1.5),
                                #panel.background = list(col = "white"),
                                par.xlab.text = list(cex = 1.5),
                                par.ylab.text = list(cex = 1.5))
    )
    dev.off()
  }
  
  ## SAVE GATING SET AS DATAFRAME
  ###############################
  metadata <- getData(samples.gating.set)[[1]]
  gsName <- description(metadata)$GUID
  getNodes( samples.gating.set, order="bfs")
  subsets <- getNodes( samples.gating.set, order="bfs", path="auto")
  df <- get.underlying.data( subsets[4:length(subsets)], samples.gating.set, metadata)
  write.csv(df,
            file = paste(myDataDir,"gatedDataFrames\\",gsName,".csv", sep=""),
            row.names=FALSE,
            quote = FALSE)
  
}