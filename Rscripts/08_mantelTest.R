rm(list=ls())
wd<-getwd()
getwd()
baseDir <- gsub("/results", "", wd)
scriptDir <- file.path(baseDir, "scripts")
#ssrData<-file.path("PATH_TO_FILE/SSR_Data.csv")
ssrData<-file.path("/Users/macbook/Documents/Wageningen/toSubmit/scripts/SSR_Data.csv")


## load packages silently
suppressPackageStartupMessages({
  library("data.table") # file reading
  library("plyr")
  library("dplyr")
  library(vegan)
  library("reshape2")
  library("tidyr")
  library("tidyverse")
  library(adegenet)
  library(ape)
  library(poppr)
  #source(file.path(scriptDir, "commonFunctions.R"))
  source(file.path("../../toSubmit/scripts", "commonFunctions.R"))
})

RE<-c("AseI", "Csp6")
for (r in 1:length(RE)){
  #genetic data
  SSR<-read.csv(ssrData, header=TRUE,  check.names = FALSE, sep=",")
  ind<-rownames(SSR)
  site<-as.character(SSR[,2])
  SSR[,1:2]<-NULL
  SSRgen<-df2genind(SSR, ploidy = 3, ind.names = ind, pop=site, sep="_")
  SSRdist<-nei.dist(SSRgen)
  
  #epigenetic data
  designTable <- file.path(paste0(baseDir, "/rawData/",RE[r], "_Design_with_converage_info.txt"))
  infileName <- file.path(paste0(baseDir,"/results/",RE[r],"_methylation.filtMETH"))
  annotationFile <- file.path(paste0(baseDir, "/rawData/",RE[r], "_mergedAnnot.csv"))
  
  sampleTab <- f.read.sampleTable(designTable) # see commonFunctions.R
  myData <- f.load.methylation.bed(infileName) # see commonFunctions.R
  myData<-unite(myData, chrPos, c(chr, pos), sep="_", remove=FALSE)
  rownames(myData)<-myData$chrPos
  ctxt<-c("CG","CHG","CHH")
  for (c in 1:length(ctxt)){
    toKeepCTXT<-ctxt[c]
    myD <- subset(myData, context == toKeepCTXT)
    #split data for control and Shade samples
    sampleNames<-sampleTab[order(rownames(sampleTab)),]
    controlSamples<-subset(sampleTab, Treat=="Control")
    shadeSamples<-subset(sampleTab, Treat=="Shade")
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
      shadeDF<-datos[rownames(shadeSamples),]
      rownames(shadeDF)<-rownames(controlDF)
      
      # distances  
      control.dis.std<-vegdist(decostand(controlDF, "norm",na.rm=TRUE), "euclidean", na.rm=TRUE)/2
      shade.dis.std<-vegdist(decostand(shadeDF, "norm",na.rm=TRUE), "euclidean", na.rm=TRUE)/2
      
      ####mantel_test
      bothTreat<-(control.dis.std+shade.dis.std)/2
      controlMantelTest<-mantel(control.dis.std, SSRdist, "pearson", permutation= 1000)
      shadeMantelTest<-mantel(shade.dis.std, SSRdist, "pearson", permutation= 1000)
      bothMantelTest<-mantel(bothTreat, SSRdist, "pearson", permutation= 1000)
      #toSave<-c(ftr[f],as.numeric(c(controlMantelTest$statistic, controlMantelTest$signif,shadeMantelTest$statistic, shadeMantelTest$signif, bothMantelTest$statistic, bothMantelTest$signif)))
      toSave<-c(ftr[f], controlMantelTest$statistic, controlMantelTest$signif)
      cat(toSave, "\n")
      out<-rbind(out,toSave)
    }
    rownames(out)<-out[,1]    
    out<-out[,-1]
    #colnames(out)<-c("r2_control", "p-value_control","r2_shade", "p-value_shade","r2_both", "p-value_both")
    colnames(out)<-c("r2", "p-value")
    out<-as.data.frame(out)
    whichAnalisis<-rep("all", times=nrow(out))
    a<-cbind(whichAnalisis, out)
    mantelTestResult<-rbind(mantelTestResult, a)   
    ######################      
    #para filtrar por DMC
    inFile<-list.files(pattern = "_table.csv")
    dfDMC<-read.csv(inFile, header=TRUE, sep="\t")
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
      shadeDF<-datos[rownames(shadeSamples),]
      rownames(shadeDF)<-rownames(controlDF)
      
      # distances  
      control.dis.std<-vegdist(decostand(controlDF, "norm",na.rm=TRUE), "euclidean", na.rm=TRUE)/2
      shade.dis.std<-vegdist(decostand(shadeDF, "norm",na.rm=TRUE), "euclidean", na.rm=TRUE)/2
      
      ####mantel_test
      bothTreat<-(control.dis.std+shade.dis.std)/2
      controlMantelTest<-mantel(control.dis.std, SSRdist, "pearson", permutation= 1000)
      shadeMantelTest<-mantel(shade.dis.std, SSRdist, "pearson", permutation= 1000)
      bothMantelTest<-mantel(bothTreat, SSRdist, "pearson", permutation= 1000)
      #toSave<-c(ftr[f],as.numeric(c(controlMantelTest$statistic, controlMantelTest$signif,shadeMantelTest$statistic, shadeMantelTest$signif, bothMantelTest$statistic, bothMantelTest$signif)))
      toSave<-c(ftr[f], bothMantelTest$statistic, bothMantelTest$signif)
      cat(toSave, "\n")
      out<-rbind(out,toSave)
    }
    
    rownames(out)<-out[,1]    
    out<-out[,-1]
    #colnames(out)<-c("r2_control", "p-value_control","r2_shade", "p-value_shade","r2_both", "p-value_both")
    colnames(out)<-c("r2", "p-value")
    out<-as.data.frame(out)
    whichAnalisis<-rep("DMC", times=nrow(out))
    b<-cbind(whichAnalisis, out)
    mantelTestResult<-rbind(mantelTestResult, b)   
    colnames(numDMC)<-c("feature","number_of_C")
    numDMC<-as.data.frame(numDMC)
    numDMC<-filter(numDMC, number_of_C!=1)
    outDF<-cbind(mantelTestResult, numDMC)
    write.table(outDF, file = paste0(RE[r],"_",ctxt[c],"_MantelTestResult.csv"), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}