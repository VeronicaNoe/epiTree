rm(list=ls())
wd<-getwd()
baseDir <- gsub("/results", "", wd)
scriptDir <- file.path(baseDir, "scripts")

## load packages silently
suppressPackageStartupMessages({
  library("data.table") # file reading
  library("DSS")
  library("qqman")
  library("plyr")
  library("dplyr")
  library(vegan)
  library("reshape2")
  library("tidyr")
  library("tidyverse")
  source(file.path(scriptDir, "commonFunctions.R"))
})

RE<-c("AseI", "Csp6")
for (r in 1:length(RE)){
  designTable <- file.path(paste0(baseDir, "/rawData/",RE[r], "_Design_with_converage_info.txt"))
  infileName <- file.path(paste0(baseDir,"/results/",RE[r],"_methylation.filtMETH"))
  annotationFile <- file.path(paste0(baseDir, "/rawData/",RE[r], "_mergedAnnot.csv"))
  sampleTab <- f.read.sampleTable(designTable) # see commonFunctions.R
  myData <- f.load.methylation.bed(infileName) # see commonFunctions.R
  
  totalCols <- grep("_total$", colnames(myData), value = TRUE)
  methCov <- myData[,totalCols]
  colnames(methCov) <- gsub("_total$", "", colnames(methCov))
  
  bins<-c(0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.85,0.9, 0.95,0.96,0.97,0.98, 0.99, 0.991, 0.995, 0.998,0.999,0.9991, 0.9995,0.9999, 0.99999, 1)
  quaCov<-quantile(methCov, bins, na.rm=TRUE) # to detect high coverage reads
  print(quaCov)
  methMelt<-reshape2::melt(methCov, variable="Samples", value.name="Tot_Methylation")
  print(mean(methMelt$Tot_Methylation, na.rm=TRUE))
  }
  