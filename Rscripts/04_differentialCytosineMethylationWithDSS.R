## Set path to input files:
rm(list=ls())
wd<-getwd()
getwd()
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

## process both data set
RE<-c("AseI", "Csp6")
for (r in 1:length(RE)){
  designTable <- file.path(paste0(baseDir, "/rawData/",RE[r], "_Design_withPlotInfos.txt"))
  infileName <- file.path(paste0(baseDir,"/results/",RE[r],"_methylation.filtMETH"))
  annotationFile <- file.path(paste0(baseDir, "/rawData/",RE[r], "_mergedAnnot.csv"))
  sampleTab <- f.read.sampleTable(designTable) # see commonFunctions.R
  Data <- f.load.methylation.bed(infileName) # see commonFunctions.R
  Trt<-c()
  Ac<-c()
  feature <- c("gene", "transposon", "repeat", "nothing")
  ctxt <- c("CHH", "CG", "CHG")
  dataFeatCtxt<-c()
  dfDMC.treat<-NULL
  dfDMC.acc<-NULL
  for (i in 1:length(feature)){
    subAnno <- f.load.merged.annotation(annotationFile, feature[i]) # see commonFunctions.R
    if (feature[i] != "all") {
      toKeep <- gsub("chr", "", rownames(subAnno))
      before <- nrow(Data)
      commonChr <- sort(intersect(toKeep, Data$chr))
      myData <- Data[as.character(Data$chr) %in% as.character(commonChr),]
      afterFeature <- nrow(myData)
      cat(paste("Subsetting by feature: ",feature[i], dim(myData)[1], "rows \n"))
      }
    for (j in 1:length(ctxt)){
      if (ctxt[j] != "all") {
        contextFilter <- ctxt[j]
        myD <- subset(myData, context == contextFilter)
        }
      allSamples <- gsub("_total$", "", grep("_total$", colnames(myD), value = TRUE))
      sampleTab <- sampleTab[allSamples,]
      myD$chr <- as.numeric(myD$chr)
      forDSS <- list()
      for (curSample in allSamples) {
        tempTab <- data.frame(
          chr = myD$chr,
          pos = myD$pos,
          N = myD[[paste0(curSample, "_total")]],
          X = myD[[paste0(curSample, "_methylated")]],
          stringsAsFactors = FALSE
          )
        forDSS[[curSample]] <- tempTab
        }
      myBS <- makeBSseqData(forDSS, names(forDSS))
      myFit <- DMLfit.multiFactor(myBS, sampleTab, formula=~Treat+Acc)
      cat(paste(ctxt[j],dim(myD)[1], "rows \n"))
      testRes.Treat <- DMLtest.multiFactor(myFit, term="Treat")
      testRes.Acc <- DMLtest.multiFactor(myFit, term="Acc")
      write.csv(testRes.Treat, paste0(baseDir,"/tmp/",RE[r],"_","Treat_",ctxt[j],"_", feature[i], "_DMC_analysis.csv"),row.names = FALSE)
      write.csv(testRes.Acc, paste0(baseDir,"/tmp/",RE[r],"_","Acc_",ctxt[j],"_", feature[i],"_DMC_analysis.csv"),row.names = FALSE)
    }
  }
}
#### summary table
inFiles <- list.files(path=paste0(baseDir,"/tmp"), pattern = "_DMC_analysis.csv$")
toSave<-c()
outdf<-matrix(NA, nrow=length(inFiles), ncol = 7)
  for (i in 1:length(inFiles)){
    minFDR<-0.05
    input<-read.csv(paste0(baseDir,"/tmp/",inFiles[i]),header=TRUE, stringsAsFactors = FALSE, sep=",")
    fileNames<-strsplit(inFiles[i], "_" )
    subDF<-filter(input, fdrs <= 0.05)
    colnames(outdf)<-c("RE","Factor","Context","Feature","# DMC","Total_Cs","Region")
    outdf[i,1]<-fileNames[[1]][1] #which RE
    outdf[i,2]<-fileNames[[1]][2] #which factor
    outdf[i,3]<-fileNames[[1]][3] #which context
    outdf[i,4]<-fileNames[[1]][4] #which feature
    outdf[i,5]<-nrow(subDF)
    outdf[i,6]<-nrow(input)
    uniReg<-unique(input$chr)
    temp<-matrix(NA, nrow=length(uniReg), ncol = 2)  
    for(j in 1:length(uniReg)){
      hits<-sum(subDF$chr==uniReg[j])
      colnames(temp) <- c("chr","ocurrences")
      temp[j,1]<-uniReg[j]
      temp[j,2]<-hits
    }
    temp<-data.frame(temp)
    temp <- temp[order(-temp$ocurrences),]
    temp<-dplyr::filter(temp, ocurrences>=5)
    if (dim(temp)[1]!= 0){
      outdf[i,7]<-"yes"
      write.table(temp,paste0(inFiles[i], "_DMC_region.csv"),row.names = FALSE, sep="\t", col.names=T, quote=FALSE)
    } else {
      outdf[i,7]<-"no"
    }
  }
write.table(outdf,paste0("00_DMC_summary.csv"),row.names = FALSE, sep="\t", col.names=T, quote=FALSE)

#### unique table
inFiles <- list.files(path=paste0(baseDir,"/tmp"), pattern = "_DMC_analysis.csv$")
toSave<-c()
for (i in 1:length(inFiles)){
  input<-read.csv(paste0(baseDir,"/tmp/",inFiles[i]), header=TRUE, stringsAsFactors = FALSE, sep=",")
  fileNames<-strsplit(inFiles[i], "_" )
  input$RE<-rep(fileNames[[1]][1], times=nrow(input))
  input$factor<-rep(fileNames[[1]][2], times=nrow(input))
  input$context<-rep(fileNames[[1]][3], times=nrow(input))
  input$feature<-rep(fileNames[[1]][4], times=nrow(input))
  toSave<-rbind(toSave, input)
}
write.table(toSave, "00_DMC_table.csv",row.names = FALSE, sep="\t", col.names=T, quote=FALSE)
