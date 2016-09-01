## OPENCYTO PACKAGE PROCESSING FACS FILES
#########################################

rm(list=ls())

library(openCyto)

setwd("C:\\Users\\Laura\\Desktop\\Alice\\")
fcsFile <- "Panel 1_HMB82-3_001.fcs"

## READ IN T CELL GATING TEMPLATES
##################################
bcell.gating.template <- gatingTemplate(paste(getwd(),"/gatingTemplates/","gatingtemplate_bcell.csv",sep=""))

## READ IN FCS FILES, SPLITTING FILES INTO ISOTYPE CONTROL AND OTHER FILES
##########################################################################
samplesFlowSet <- read.flowSet(paste("data\\Blood\\Healthy\\",fcsFile,sep=""),alter.names=TRUE)
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
#plot(samples.gating.set)
#plotGate(samples.gating.set[[1]],xbin=50,default.y="SSC.A")
subsets <- getNodes(samples.gating.set[[1]],order="bfs")
for (i in subsets[c(9,10,11,15,16)]){
  name = strsplit(i,split="/")[[1]]
  jpeg(paste(name[length(name)],'.jpg',sep=""))
  plotGate(samples.gating.set[[1]], i,
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

## SAVE GATING SET TO DISK
##########################
#save_gs(samples.gating.set,path=paste(getwd(),"/gatingSets/",sampleNames(samplesFlowSet),sep=""))
