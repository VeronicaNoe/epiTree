## Set path to input files:
rm(list=ls())
wd<-getwd()
baseDir <- gsub("/<PATH-TO-results>", "", wd)
scriptDir <- file.path(baseDir, "scripts")

## load packages silently
suppressPackageStartupMessages({
  library("data.table")
  library("vegan")
  source(file.path(scriptDir, "commonFunctions.R"))
})

## process both data set
RE<-c("AseI-NsiI", "Csp6I-NsiI")
for (r in 1:length(RE)){
  designTable <- file.path(paste0(baseDir, "/<PATH-TO-rawData>/",RE[r], "_Design_withPlotInfos.txt"))
  infileName <- file.path(paste0(baseDir,"/<PATH-TO-results/",RE[r],"_methylation.filtered"))
  annotationFile <- file.path(paste0(baseDir, "/<PATH-TO-annotation>/",RE[r], "_mergedAnnot.csv"))
  
  colNamesForGrouping <- c("Treat")
  sampleTab <- f.read.sampleTable(designTable, colNamesForGrouping) # see commonFunctions.R
  Data <- f.load.methylation.bed(infileName) # see commonFunctions.R
  ctxt <- c("CHH", "CG", "CHG")#, "all")
  #loop for all context
  for (j in 1:length(ctxt)){
    myD <- subset(Data, context == ctxt[j])
    totalCols <- grep("_total$", colnames(myD), value = TRUE)
    methCols <- grep("_methylated$", colnames(myD), value = TRUE)
    totCov <- myD[,totalCols]
    methCov <- myD[,methCols]
    colnames(totCov) <- gsub("_total$", "", colnames(totCov))
    colnames(methCov) <- gsub("_methylated$", "", colnames(methCov))
    # match the samples
    commonSamples <- sort(intersect(colnames(totCov), rownames(sampleTab)))
    out<-f.meth.distances(methCov, totCov)
    #plot the non-Metric MDS
    m<-as.matrix(t(out))
    d<-as.dist(m)
    forPlot <- tryCatch(MASS::isoMDS(d)$points, error = function(e) {NA}, finally = cat("###\n"))
    pdf(paste0(baseDir,"/<PATH-TO-plots>/", RE[r], "_",ctxt[j],"_isoMDS", ".pdf"), height = 4, width = 4)
    plot(forPlot, pch=sampleTab[rownames(forPlot), "pch"], col = sampleTab[rownames(forPlot), "color"], main = ctxt[j])
    invisible(dev.off())
    #permanova in distance matrix
    eti<-sampleTab[order(rownames(sampleTab)),1:3]
    eti[, 1:3]<-lapply(eti[,1:3], as.factor)
    fit<-adonis(d ~eti$Acc+eti$Treat, data=eti, permutation=10000)
    #save results
    write.csv(fit$aov.tab, file = paste0(baseDir,"/<PATH-TO-tmp>/",RE[r],"_",ctxt[j], "_adonis.csv"), row.names = TRUE)
    write.csv(out, file = paste0(baseDir,"/<PATH-TO-tmp>/",RE[r],"_",ctxt[j], "_Distances.csv"), row.names = TRUE)
  }
}

  