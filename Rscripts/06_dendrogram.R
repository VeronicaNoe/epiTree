## Set path to input files:
rm(list=ls())
wd<-getwd()
baseDir <- gsub("/<PATH-TO-results>", "", wd)
scriptDir <- file.path(baseDir, "scripts")

## Load packages silently:
suppressPackageStartupMessages({
  library(data.table) # file reading
  library(adegenet)
  library(poppr)
  library(vegan)
  library(tidyr)
  source(file.path(scriptDir, "commonFunctions.R"), local=TRUE)
})

# for genetic data
data<-read.csv(paste0(baseDir,"/<PATH-TO-rawData>/SSR_data.csv"), header=TRUE, stringsAsFactors = FALSE, row.names = NULL, sep="\t")
str(data)
rownames(data)<-data$Pop
ind<-as.character(data[,1])
site<-as.character(data[,2])
data[,1:3]<-NULL
data_gen<-df2genind(data, ploidy = 3, ind.names = ind, pop=site, sep="_")
# plot
pdf(paste0(baseDir,"/<PATH-TO-plots>/SSR_NeiDistance_UPGMA.pdf"))
data_gen %>% 
  genind2genpop() %>%
  aboot(cutoff = 60, quiet = TRUE, sample = 1000, distance = nei.dist)
dev.off()

# for epigenetic data
## process both data set
RE<-c("AseI-NsiI", "Csp6I-NsiI")
for (r in 1:length(RE)){
  designTable <- file.path(paste0(baseDir, "/<PATH-TO-rawData>/",RE[r], "_Design_with_converage_info.txt"))
  infileName <- file.path(paste0(baseDir,"/<PATH-TO-rawData>/",RE[r],"_methylation.filtered"))
  annotationFile <- file.path(paste0(baseDir, "/<PATH-TO-annotation>/",RE[r], "_mergedAnnot.csv"))
  
  sampleTab <- f.read.sampleTable(designTable) # see commonFunctions.R
  myData <- f.load.methylation.bed(infileName) # see commonFunctions.R
  myData<-unite(myData, chrPos, c(chr, pos), sep="_", remove=FALSE)
  sort(myData$chrPos)
  rownames(myData)<-myData$chrPos
  
  inFile<-list.files(pattern = "_table.csv")
  temp<-read.csv(inFile[1], header=TRUE, sep="\t")
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
          onlyControl<-grep("Control:", rownames(datos))
          datos<-datos[onlyControl,]
          dis <- vegdist(datos, na.rm=TRUE, "euclid")
          dist<-vegdist(decostand(dis, "norm"), "euclidean", na.rm=TRUE)
          clus<-hclust(dist, "average")
          dirOut<-paste0(baseDir,"/<PATH-TO-plots>/",RE[r],"_all" ,"_Context_methylationLevelsPerFeature.pdf")
          pdf(paste0(baseDir,"/<PATH-TO-plots>/",RE[r],"_","Dendrogram_","allSamples_",context[j],"_",whichAnalysis[a],"_", feature[i], "_2.pdf"))
          plot(clus, cex = 0.6, hang=-0.5,main=paste0("Dendrogram_",RE[r],"_",context[j],"_",whichAnalysis[a],"_", feature[i]))
          dev.off()
        }
      }
    }
  }
}
