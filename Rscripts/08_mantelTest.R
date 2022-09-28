## Set path to input files:
rm(list=ls())
wd<-getwd()
baseDir <- gsub("/<PATH-TO-results>", "", wd)
scriptDir <- file.path(baseDir, "scripts")

## load packages silently
suppressPackageStartupMessages({
  library(data.table) # file reading
  library(vegan)
  library(tidyr)
  library(adegenet)
  library(poppr)
  source(file.path(scriptDir, "commonFunctions.R"))
})

#load data
ssrData<-file.path(paste0(baseDir, "/<PATH-TO-rawData>/","SSR_Data.csv"))
RE<-c("AseI-NsiI", "Csp6I-NsiI")
for (r in 1:length(RE)){
  #genetic data
  SSR<-read.csv(ssrData, header=TRUE, row.names = 1, check.names = FALSE, sep="\t")
  ind<-rownames(SSR)
  site<-as.character(SSR[,2])
  SSR[,1:2]<-NULL
  SSRgen<-df2genind(SSR, ploidy = 3, ind.names = ind, pop=site, sep="_")
  SSRdist<-nei.dist(SSRgen)
  
  #epigenetic data
  designTable <- file.path(paste0(baseDir, "/<PATH-TO-rawData>/",RE[r], "_Design_with_converage_info.txt"))
  infileName <- file.path(paste0(baseDir,"/<PATH-TO-rawData>/",RE[r],"_methylation.filtered"))
  annotationFile <- file.path(paste0(baseDir, "/<PATH-TO-annotation>/",RE[r], "_mergedAnnot.csv"))
  
  sampleTab <- f.read.sampleTable(designTable) # see commonFunctions.R
  myData <- f.load.methylation.bed(infileName) # see commonFunctions.R
  myData<-unite(myData, chrPos, c(chr, pos), sep="_", remove=FALSE)
  rownames(myData)<-myData$chrPos
  ctxt<-c("CG","CHG","CHH")
  for (c in 1:length(ctxt)){
    toKeepCTXT<-ctxt[c]
    myD <- subset(myData, context == toKeepCTXT)
    #skeep only control samples
    sampleNames<-sampleTab[order(rownames(sampleTab)),]
    controlSamples<-subset(sampleTab, Treat=="Control")
    
    #do mantel test for each feature
    ftr <- c("all","gene", "transposon", "repeat", "nothing")
    out<-c()
    mantelTestResult<-c()
    numDMC<-c()
    for (f in 1: length(ftr)){
      cat(ftr[f], "\n")
      subAnno <- f.load.merged.annotation(annotationFile, ftr[f])
      toKeep <- gsub("chr", "", rownames(subAnno))
      commonChr <- sort(intersect(toKeep, as.character(myD$chr)))
      mD <- myD[as.character(myD$chr) %in% as.character(commonChr),]
      lenDMC<-c(ftr[f],nrow(mD))
      numDMC<-rbind(numDMC, lenDMC)  
      numDMC<-data.frame(numDMC)
      
      totalCols <- grep("_total$", colnames(mD), value = TRUE)
      methCols <- grep("_methylated$", colnames(mD), value = TRUE)
      totCov <- mD[,totalCols]
      methCov <- mD[,methCols]
      colnames(totCov) <- gsub("_total$", "", colnames(totCov))
      colnames(methCov) <- gsub("_methylated$", "", colnames(methCov))
      mePerc<-methCov/totCov
      datos<-t(mePerc)
      datos<-datos[order(rownames(datos)),]
      controlDF<-datos[rownames(controlSamples),]
      
      # distances  
      control.dis.std<-vegdist(decostand(controlDF, "norm",na.rm=TRUE), "euclidean", na.rm=TRUE)/2
      
      ####mantel_test
      controlMantelTest<-mantel(control.dis.std, SSRdist, "pearson", permutation= 1000)
      toSave<-c(ftr[f], controlMantelTest$statistic, controlMantelTest$signif)
      cat(toSave, "\n")
      out<-rbind(out,toSave)
    }
    rownames(out)<-out[,1]    
    out<-out[,-1]
    colnames(out)<-c("r2", "p-value")
    out<-as.data.frame(out)
    whichAnalisis<-rep("all", times=nrow(out))
    a<-cbind(whichAnalisis, out)
    mantelTestResult<-rbind(mantelTestResult, a)   
    ######################      
    #para filtrar por DMC
    inFile<-list.files(pattern = "_table.csv")
    dfDMC<-read.csv(inFile[1], header=TRUE, sep="\t")
    dfDMC<-unite(dfDMC, chrPos, c(chr, pos), sep="_", remove=FALSE)
    toKeepRE<-RE[r]
    toKeepCTXT<-ctxt[c]
    dfDMC<-dplyr::filter(dfDMC, RE==toKeepRE)
    dfDMC<-dplyr::filter(dfDMC, factor=="Acc")
    dfDMC<-dplyr::filter(dfDMC, context==toKeepCTXT)
    dfDMC<-dplyr::filter(dfDMC, pvals<=0.05)
    rownames(dfDMC)<-dfDMC$chrPos
    
    keep<-intersect(rownames(myD),rownames(dfDMC))
    mD<-myD[keep,]
    ftr<-c("all",unique(dfDMC$feature))
    out<-c()
    
    for (f in 1: length(ftr)){
      if(ftr[f]=="all"){
        temp<-dfDMC
        tokeep<-intersect(sort(rownames(mD)), sort(rownames(temp)))
        df<-mD[tokeep,]
        lenDMC<-c(ftr[f],length(tokeep))
        numDMC<-rbind(numDMC,lenDMC)
      }else{
        temp<-dplyr::filter(dfDMC, feature==ftr[f])
        tokeep<-intersect(sort(rownames(mD)), sort(rownames(temp)))
        cat(ftr[f],length(tokeep),"\n")
        lenDMC<-c(ftr[f],length(tokeep))
        numDMC<-rbind(numDMC,lenDMC)
        if(length(tokeep)<=1){
          next
        }else{
          df<-mD[tokeep,]
        }
      }
      
      totalCols <- grep("_total$", colnames(df), value = TRUE)
      methCols <- grep("_methylated$", colnames(df), value = TRUE)
      totCov <- df[,totalCols]
      methCov <- df[,methCols]
      colnames(totCov) <- gsub("_total$", "", colnames(totCov))
      colnames(methCov) <- gsub("_methylated$", "", colnames(methCov))
      mePerc<-methCov/totCov
      datos<-t(mePerc)
      datos<-datos[order(rownames(datos)),]
      controlDF<-datos[rownames(controlSamples),]
      
      # distances  
      control.dis.std<-vegdist(decostand(controlDF, "norm",na.rm=TRUE), "euclidean", na.rm=TRUE)
      
      ####mantel_test
      controlMantelTest<-mantel(control.dis.std, SSRdist, "pearson", permutation= 1000)
      toSave<-c(ftr[f], controlMantelTest$statistic, controlMantelTest$signif)
      cat(toSave, "\n")
      out<-rbind(out,toSave)
    }
    
    rownames(out)<-out[,1]    
    out<-out[,-1]
    colnames(out)<-c("r2", "p-value")
    out<-as.data.frame(out)
    whichAnalisis<-rep("DMC", times=nrow(out))
    b<-cbind(whichAnalisis, out)
    mantelTestResult<-rbind(mantelTestResult, b)   
    colnames(numDMC)<-c("feature","number_of_C")
    numDMC$number_of_C<-as.numeric(numDMC$number_of_C)
    numDMC<-subset(numDMC, numDMC$number_of_C >1)
    outDF<-cbind(mantelTestResult, numDMC)
    write.table(outDF, file = paste0(RE[r],"_",ctxt[c],"_MantelTestResult.csv"), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}
