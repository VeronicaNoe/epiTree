## Set path to input files:
rm(list=ls())
wd<-getwd()
baseDir <- gsub("/results", "", wd)
scriptDir <- file.path(baseDir, "scripts")

## Load packages silently:
suppressPackageStartupMessages({
  library(data.table)
  library(vioplot) # plotting
  library(vegan)
  source(file.path(scriptDir, "commonFunctions.R"))
})

## process both data set
RE<-c("AseI", "Csp6")
for (r in 1:length(RE)){
  designTable <- file.path(paste0(baseDir, "/rawData/",RE[r], "_Design_withPlotInfos.txt"))
  infileName <- file.path(paste0(baseDir,"/results/",RE[r],"_methylation.filtMETH"))
  annotationFile <- file.path(paste0(baseDir, "/rawData/",RE[r], "_mergedAnnot.csv"))
  ## Load data
  sampleTab <- f.read.sampleTable(designTable)
  mePerc <- f.load.methylation.bed(infileName, percentages = TRUE) # see commonFunctions.R
  infoColumns <- c("chr", "pos", "context")
  allSamples <- setdiff(colnames(mePerc), infoColumns)
  sampleTab <- sampleTab[allSamples,]

  ## Average within groups 
  #Select on what to average
  aveData <- f.summarize.columns(mePerc, data.frame(sample = rownames(sampleTab), group = sampleTab$Treat, stringsAsFactors = FALSE), function(x) mean(x, na.rm = TRUE))
  aveDataInfo <- mePerc[,infoColumns]
  rownames(aveDataInfo) <- paste0("chr", aveDataInfo$chr, "_", aveDataInfo$pos)
  rownames(aveData) <- rownames(aveDataInfo)
  ## Choose a group order for the plot and set the colors
  forPlotOrder <- c("Control", "Shade")
  aveData <- aveData[, match(forPlotOrder, colnames(aveData))]
  temp <- unique(sampleTab[,c("Treat","color")])
  plotColors <- temp$color; names(plotColors) <- temp$Treat
  ## Draw histograms
  allContexts <- c("CG", "CHG", "CHH")#
  numPlotRows <- 1 #
  numPlotCols <- length(allContexts)
  pdf(file.path(paste0(baseDir,"/plots/", RE[r], "_histo.pdf")))
  par(oma = c(2, 2, 2, 2))
  layout(matrix(1:(numPlotRows*numPlotCols), nrow = numPlotRows, byrow = TRUE))
  for (ctxt in allContexts) {
    subData <- subset(mePerc, mePerc$context==ctxt)
    toPlot<-reshape2::melt(subData, id=infoColumns)
    histo<-hist(toPlot$value, breaks=seq(0,100,10),  plot=FALSE)
    ymax<-max(histo$counts)
    plot(histo, main=ctxt, xlab="", ylim=c(0, ymax*1.5))
  }
  invisible(dev.off())

  ## Draw violin plots
  allContexts <- c("CG", "CHG", "CHH")#, "all") 
  numPlotRows <- 1 
  numPlotCols <- length(allContexts)
  allMeans <- matrix(NA, nrow = length(forPlotOrder), ncol = numPlotCols, dimnames = list(forPlotOrder, allContexts))
  aveData <- aveData[rownames(aveDataInfo),] 
  
  pdf(paste0(baseDir,"/plots/", RE[r], "_Context_methylationLevelsViolinPlot.pdf"), height = 5*numPlotRows, width = 2+length(forPlotOrder)*numPlotCols)
  par(oma = c(12, 8, 3, 0), mar = c(0, 0, 0, 0))
  layout(matrix(1:(numPlotRows*numPlotCols), nrow = numPlotRows, byrow = TRUE))
  for (ctxt in allContexts) {
    if (ctxt == "all") {
      subData <- aveData
    } else {
      subData <- aveData[aveDataInfo$context == ctxt,]
    }
    plot(NA, main = ctxt, bty = "n", xaxs = "r", yaxs = "r", xlab = "", ylab = "", las = 1, cex = 0.2, tck = 0.01, xlim = c(0.5, length(forPlotOrder)+0.5), ylim = c(0, 100), xaxt = "n", yaxt = "n")
    curPos <- 1
    for (curGroup in forPlotOrder) {
    toPlot <- subData[,curGroup]
    toPlot <- toPlot[!is.na(toPlot)]
    curCol <- plotColors[curGroup]
    if (sum(toPlot > 0) > 4) {
      vioplot(toPlot, names = c(curGroup), col = curCol, ylim = c(0,100), drawRect = TRUE, add = TRUE, at = curPos)
      }
    curMean <- mean(toPlot)
    lines(c(curPos-0.3,curPos+0.3), c(curMean, curMean), col = "black", lwd = 4, lty = 1)
    curPos <- curPos + 1
    allMeans[curGroup,ctxt] <- curMean # add the mean to the collection
    }
    if (ctxt != "all") { axis(2, at = seq(0, 100, by = 20), labels = seq(0, 100, by = 20), outer = TRUE, las = 1, line=2, lwd=2, cex.axis=2) }
    axis(1, at = 1:length(forPlotOrder), labels = forPlotOrder, outer = TRUE, las = 2, line=2, lwd=2, cex.axis=3)
  }
  invisible(dev.off())
  write.csv(round(allMeans, 3), file.path(paste0(baseDir,"/tmp/",RE[r], "_Context_methylationLevelsViolinPlot_means.csv")))

  ## Get the average methylation level per group, context, feature
  allFeatures <- c("gene", "transposon", "repeat", "nothing")
  forMask <- paste0("chr", aveDataInfo$chr)
  listForPlot <- list()
  
  for (ctxt in allContexts) {
    if (ctxt == "all") {
      contextMask <- rep(TRUE, nrow(aveDataInfo))
    } else {
      contextMask <- aveDataInfo$context == ctxt
    }
    featureMeans <- matrix(NA, nrow = length(allFeatures), ncol = length(forPlotOrder), dimnames = list(allFeatures, forPlotOrder))
    for (feature in allFeatures) {
      mergedAnno <- f.load.merged.annotation(annotationFile, feature)
      annoMask <- forMask %in% rownames(mergedAnno)
      subData <- aveData[annoMask & contextMask,]
      featureMeans[feature, colnames(aveData)] <- colMeans(subData, na.rm = TRUE)
    }
    listForPlot[[ctxt]] <- featureMeans
  }
  
  ## Do the plot
  dirOut<-paste0(baseDir,"/plots/",RE[r], "_Context_methylationLevelsPerFeature.pdf")
  imageColors <- f.blackblueyellowredpinkNICE(51) 
  pdf(dirOut, height = 5*numPlotRows, width = 2+length(forPlotOrder)*numPlotCols)
  layout(matrix(1:(numPlotRows*numPlotCols), nrow = numPlotRows, byrow = TRUE))
  for (ctxt in allContexts) {
    temp <- listForPlot[[ctxt]]
    f.image.without.text(forPlotOrder, allFeatures, t(temp), xLabel = "", yLabel = "", mainLabel = ctxt, useLog = FALSE, col = imageColors, zlim = c(0, 100))
    write.csv(round(temp, 3), paste0(baseDir,"/tmp/",RE[r],"_methylationLevelsPerFeature_",ctxt, ".csv"))
    }
  invisible(dev.off())
}