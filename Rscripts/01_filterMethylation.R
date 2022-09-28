## Set path to input files:
rm(list=ls())
wd<-getwd()
baseDir <- gsub("/<PATH-TO-results>", "", wd)
scriptDir <- file.path(baseDir, "scripts")
## Load packages silently:
suppressPackageStartupMessages({
  library("data.table")
  library("plyr") #ddply
  library("dplyr") #select
  library("reshape2") #melt
  library("tidyverse") #mutate
  library("ggplot2")
  source(file.path(scriptDir, "commonFunctions.R"))
})

## process both data set
RE<-c("AseI-NsiI", "Csp6I-NsiI")
for (r in 1:length(RE)){
  designTable <- file.path(paste0(baseDir, "/<PATH-TO-rawData>/",RE[r], "_Design_withPlotInfos.txt"))
  infileName <- file.path(paste0(baseDir,"/<PATH-TO-rawData>/",RE[r],"_methylation.bed"))
  annotationFile <- file.path(paste0(baseDir, "/<PATH-TO-annotation>/",RE[r], "_mergedAnnot.csv"))
  
  rDir <- file.path(baseDir)
  ## Load and explore data
  myData <- f.load.methylation.bed(infileName) 
  colNamesForGrouping <- c("Treat")
  sampleTab <- f.read.sampleTable(designTable, colNamesForGrouping) # see commonFunctions.R
  str(myData)
  totalCols <- grep("_total$", colnames(myData), value = TRUE)
  methCov <- myData[,totalCols]
  colnames(methCov) <- gsub("_total$", "", colnames(methCov))
  
  #Number of Citosynes (C) in original data set. 
  dimA<-dim(myData)[1]
  
  ## Match to remove bad samples
  commonSamples <- sort(intersect(colnames(methCov), rownames(sampleTab)))
  #allSamples <- union(colnames(methCov), rownames(sampleTab))
  samplesToRemove <- setdiff(colnames(methCov), commonSamples) #order mathers! biggest first
  if (length(samplesToRemove) > 0) {
    f.print.message("Removing", length(samplesToRemove), "samples!")
    cat(paste0(samplesToRemove, collapse = '\n'), '\n')
  }
  ### remove bad samples
  sampleTab <- sampleTab[commonSamples,]
  dim(methCov)
  methCov <- methCov[, commonSamples]
  myData <- myData[, c("chr", "pos", "context", "samples_called", paste0(rep(commonSamples, each = 2), c("_methylated", "_total")))]

  ## Filtering data by minimum and maximum coverage
  #First start with minimum coverage
  bins<-c(0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9, 0.95,0.96,0.97,0.98, 0.99, 0.991, 0.995, 0.998,0.999,0.9991, 0.9995,0.9999, 0.99999, 1)
  quaCov<-quantile(methCov, bins, na.rm=TRUE) # to detect high coverage reads
  #transform in NA all values below 10X
  minCov<-10
  methCov[(methCov < minCov)] <- NA

  quaCovB<-quantile(methCov, bins, na.rm=TRUE) # to detect high coverage reads
  #set number of samples that are the 80% of samples
  minCountPerGroup <- round((length(commonSamples)/2)*0.8, 0)
  tabForSummary <- data.frame(sample = rownames(sampleTab), group = sampleTab$group, stringsAsFactors = FALSE)
  callsPerGroup <- f.summarize.columns(!is.na(methCov), tabForSummary, sum)
  maskToKeep <- rowSums(callsPerGroup >= minCountPerGroup) == ncol(callsPerGroup) # c's that are in at least 80% samples

  if (sum(maskToKeep) == 0) {
    f.print.message("No cytosine passed the filter, no output.")
    quit("no", 0)
  }

  #remove C's that arenot in 80% samples
  myData <- myData[maskToKeep,]
  #subset and store

  dimB<-dim(myData)[1]
  dimB

  #Second, set maximum filtering after removing minimum filtering
  totalCols <- grep("_total$", colnames(myData), value = TRUE)
  methCov <- myData[,totalCols]
  colnames(methCov) <- gsub("_total$", "", colnames(methCov))

  maxCov<-as.numeric(quaCovB[25]) #99.99%
  methCov[(methCov > maxCov)] <- NA
  minCountPerGroup <- round((length(commonSamples)/2)*0.8, 0)

  tabForSummary <- data.frame(sample = rownames(sampleTab), group = sampleTab$group, stringsAsFactors = FALSE) #summarize data for each sample
  callsPerGroup <- f.summarize.columns(!is.na(methCov), tabForSummary, sum)
  maskToKeep <- rowSums(callsPerGroup >= minCountPerGroup) == ncol(callsPerGroup) # c's that are in at least 80% samples

  if (sum(maskToKeep) == 0) {
    f.print.message("No cytosine passed the filter, no output.")
    #quit("no", 0)
  }

  #subset and store
  myData <- myData[maskToKeep,]
  write.table(myData, file = paste0(baseDir,"/<PATH-TO-results/",RE[r],"_methylation.filtered"), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  #C's in final data
  totalCols <- grep("_total$", colnames(myData), value = TRUE)
  methCov <- myData[,totalCols]
  methCov[(methCov < minCov) | (methCov > maxCov)] <- NA
  dimC<-dim(myData)[1]
  quaCovC<-quantile(methCov, bins, na.rm=TRUE)
  length(unique(myData$chr)) 

  #see quantile for each sample
  qua4samples<-c()
  for (s in 1:ncol(methCov)){
    qDis<-quantile(methCov[,1], na.rm=TRUE)
    sampleDis<-c(colnames(methCov)[s], qDis)
    qua4samples<-rbind(qua4samples, sampleDis)
  }

  # Save results
  sink(paste0(baseDir,"/<PATH-TO-tmp/", RE[r],"_Filtering_report.txt"))
  paste0("Number of Cs from bed file: ", dimA)
  paste0("Quantile coverage from bed file: ")
  quaCov
  paste0("Number of Cs after 10x filtering: ", dimB)
  paste0("Quantile coverage after 10X: ")
  quaCovB
  paste0("Number of Cs final file: ", dimC)
  paste0("Quantile coverage final file: ")
  quaCovC
  sink()
  # check coverage
  ctxt<-c("CG", "CHG", "CHH")
  tmp <- myData[,c("context",totalCols)]
  qDisTotal<-quantile(tmp[,2:dim(tmp)[2]], na.rm=TRUE)
  totalMean<-melt(tmp[,2:dim(tmp)[2]],)
  totalKeep<-mean(totalMean[,2], na.rm=TRUE)
  qDisTotal<-c(data.frame(qDisTotal)[,1], totalKeep)
  out<-c()
  for(c in ctxt){
    tmp2<-subset(tmp, context==c)
    toMean<-melt(tmp2[,2:dim(tmp2)[2]],)
    toKeep<-mean(toMean[,2], na.rm=TRUE)
    qDis<-quantile(tmp2[,2:dim(tmp2)[2]], na.rm=TRUE)
    count<-data.frame(qDis)[,1]
    tmp_out<-c(count, toKeep)
    out<-cbind(out,tmp_out)
  }
  out<-cbind(out,qDisTotal)
  rownames(out)<-c(rownames(data.frame(qDis)),"mean")
  colnames(out)<-c(ctxt,"qDisTotal")
  write.table(out, file = paste0(baseDir,"/<PATH-TO-results/",RE[r],"_qDis.txt"), sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
}