rm(list=ls())
wd<-getwd()
getwd()
baseDir <- gsub("/results", "", wd)
scriptDir <- file.path(baseDir, "scripts")

## load packages silently
suppressPackageStartupMessages({
  library(data.table) # file reading
  library(qqman)
  library(plyr)
  library(dplyr)
  library(tidyr)
  source(file.path(scriptDir, "commonFunctions.R"))
})

## process both data set
RE<-c("AseI-NsiI", "Csp6I-NsiI")
infileName <- file.path(paste0(baseDir,"/rawData/","00_DMC_table.csv"))
data<-read.csv(infileName, header=TRUE, sep="\t")
head(data)
#for (r in 1:1){
for (r in 1:length(RE)){
  toUse<-RE[r]
  dfDMC<-dplyr::filter(data, RE==toUse)
  dfDMC<-dplyr::filter(dfDMC, fdrs!="NA")
  dfDMC<-dplyr::filter(dfDMC, factor=="Treat")
  dfDMC<-unite(dfDMC, chrPos, c(chr, pos), sep="_", remove=FALSE)
  dfDMC$chr<-as.numeric(dfDMC$chr)
  
  DMC<-dplyr::filter(dfDMC, fdrs<=0.05)
  uniReg<-unique(DMC$chr)
  
  outdf<-matrix(NA, nrow=length(uniReg), ncol = 2)  
  for(i in 1:length(uniReg)){
    hits<-sum(DMC$chr==uniReg[i])
    colnames(outdf) <- c("chr","ocurrences")
    outdf[i,1]<-uniReg[i]
    outdf[i,2]<-hits
  }
  outdf<-data.frame(outdf)
  outdf <- outdf[order(-outdf$ocurrences),]
  out<-dplyr::filter(outdf, ocurrences>=5)
  out<-out[order(out$chr),]
  write.csv(out, paste0(baseDir,"/results/",RE[r],"_DMRegions.csv"),row.names = FALSE)
  
  outDir<-paste0(baseDir,"/plots/", RE[r],"_")
  ctxt <- c("CG","CHG","CHH")
  for (i in 1:length(ctxt)){
    df<-dplyr::filter(dfDMC, context==ctxt[i])
    sum(is.na(df))
    toPlot<-dplyr::filter(DMC, DMC$context==ctxt[i])
    snpOfInterest<-intersect(toPlot$chr,out$chr)
    snpOfInterest<-sort(snpOfInterest)
    pdf(paste0(outDir,ctxt[i],"_manhattanPlot.pdf"))
    manhattan(df, chr = "chr", bp = "pos", p = "fdrs", highlight =snpOfInterest, 
              snp = "chr", annotatePval = 0.05, annotateTop=TRUE, col = c("gray60"), 
              chrlabs = NULL, main = paste0("DMC_for_",ctxt[i], "_context"," (",RE[r],")"), 
              suggestiveline = -log10(5e-02), xlab="epiGBS fragment", ylab=expression('-log'[10]*' (FDR)'))
    dev.off()
    }
  }
