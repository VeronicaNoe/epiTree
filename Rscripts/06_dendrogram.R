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
# process both data set
RE<-c("AseI", "Csp6")
for (r in 1:length(RE)){
  designTable <- file.path(paste0(baseDir, "/rawData/",RE[r], "_Design_with_converage_info.txt"))
  infileName <- file.path(paste0(baseDir,"/results/",RE[r],"_methylation.filtMETH"))
  annotationFile <- file.path(paste0(baseDir, "/rawData/",RE[r], "_mergedAnnot.csv"))
  sampleTab <- f.read.sampleTable(designTable) # see commonFunctions.R
  myData <- f.load.methylation.bed(infileName) # see commonFunctions.R
  myData<-unite(myData, chrPos, c(chr, pos), sep="_", remove=FALSE)
  sort(myData$chrPos)
  rownames(myData)<-myData$chrPos

  inFile<-list.files(pattern = "_table.csv")
  temp<-read.csv(inFile, header=TRUE, sep="\t")
  temp<-unite(temp, chrPos, c(chr, pos), sep="_", remove=FALSE)
  toKeepRE<-RE[r]
  temp<-dplyr::filter(temp, RE==toKeepRE)
  temp<-dplyr::filter(temp, factor=="Acc")
  temp<-dplyr::filter(temp, pvals<=0.05)
  sort(temp$chrPos)
  rownames(temp)<-temp$chrPos

  context <- c("CG", "CHH", "CHG")
  whichAnalysis<-c("all", "DMC")
  for (j in 1:length(context)){
    if (context[j] != "all") {
      contextFilter <- context[j]
      myD <- subset(myData, context == contextFilter)
      } else if (context[j] == "all"){
      myD<-myData
      }
  # subset meth data that has DMC
    for (a in 1:length(whichAnalysis)){
      if(whichAnalysis[a]=="all"){
        mD<-myD
        feature <- c("all","gene", "transposon", "repeat", "nothing")
      } else {
        keep<-intersect(sort(rownames(temp)),sort(rownames(myD)))
        mD<-myD[keep,]
        feature<-c("all",unique(temp$feature))
      }
  ## a plot for each feature
      for (i in 1:length(feature)){
        subAnno <- f.load.merged.annotation(annotationFile, feature[i])
        toKeep <- gsub("chr", "", rownames(subAnno))
        commonChr <- sort(intersect(toKeep,as.character(mD$chr)))
        df <- mD[as.character(mD$chr) %in% as.character(commonChr),]
        if(dim(df)[1]<=4){
          next
        }else{
          ##dendrogram
          totalCols <- grep("_total$", colnames(df), value = TRUE)
          methCols <- grep("_methylated$", colnames(df), value = TRUE)
          totCov <- df[,totalCols]
          methCov <- df[,methCols]
          colnames(totCov) <- gsub("_total$", "", colnames(totCov))
          colnames(methCov) <- gsub("_methylated$", "", colnames(methCov))
          mePerc<-methCov/totCov
          datos<-t(mePerc) 
          datos<-datos[order(rownames(datos)),]
          sampleNames<-sampleTab[order(rownames(sampleTab)),]
          sampleNames$Sample_name<-rownames(sampleNames)
          sampleNames<-unite(sampleNames, roName, Treat,Acc, sep=":", remove=TRUE)
          rownames(datos)<-sampleNames$roName
          #onlyControl<-grep("Control:", rownames(datos))
          #datos<-datos[onlyControl,]
          dis <- vegdist(datos, na.rm=TRUE, "euclid")
          dist<-vegdist(decostand(dis, "norm"), "euclidean", na.rm=TRUE)
          clus<-hclust(dist, "average")
          dirOut<-paste0(baseDir,"/plots/",RE[r],"_all" ,"_Context_methylationLevelsPerFeature.pdf")
          pdf(paste0(baseDir,"/plots/",RE[r],"_","Dendrogram_","allSamples_",context[j],"_",whichAnalysis[a],"_", feature[i], ".pdf"))
          plot(clus, cex = 0.6, hang=-0.5,main=paste0("Dendrogram_",RE[r],"_",context[j],"_",whichAnalysis[a],"_", feature[i]))
          dev.off()
          }
       }
    }
  }
}

### SSR dendrogram
data<-read.csv(paste0("<PATH_TO_FILE>/SSR_Data.csv"), header=TRUE, stringsAsFactors = FALSE, row.names = NULL, sep=",")
str(data)
rownames(data)<-data$Pop
data$Pop<-NULL
data.sim<-vegdist(data,  na.rm=TRUE, "euclid")
clus<-hclust(data.sim, "average")
pdf("Dendrogram.pdf")
pdf(paste0(baseDir,"/plots/","SSR_Dendrogram.pdf"))
plot(clus, hang = -1, cex = 0.6, main="SSR_dendrogram", ylab="Euclidean distance")
dev.off()
