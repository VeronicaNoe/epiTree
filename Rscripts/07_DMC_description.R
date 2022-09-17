rm(list=ls())
wd<-getwd()
getwd()
baseDir <- gsub("/results", "", wd)
scriptDir <- file.path(baseDir, "scripts")

suppressPackageStartupMessages({
  library("data.table") # file reading
  library("vioplot") # plotting
  library("tidyverse")
  library("plyr") #ddply
  library(emmeans)
  library(car) #anova
  source(file.path(scriptDir, "commonFunctions.R"))
})

####### process both data set
RE<-c("AseI", "Csp6")
for (r in 1:length(RE)){
#for (r in 1:1){
  designTable <- file.path(paste0(baseDir, "/rawData/",RE[r], "_Design_withPlotInfos.txt"))
  infileName <- file.path(paste0(baseDir,"/results/",RE[r],"_methylation.filtMETH"))
  annotationFile <- file.path(paste0(baseDir, "/rawData/",RE[r], "_mergedAnnot.csv"))
  mePerc <- f.load.methylation.bed(infileName, percentages = TRUE) # see commonFunctions.R
  sampleTab <- f.read.sampleTable(designTable)
  mePerc<-unite(mePerc, chrPos, c(chr, pos), sep="_", remove=FALSE)
  mePerc<-cbind(mePerc[,1:4],mePerc[,rownames(sampleTab)])
  colNam<-c(colnames(mePerc[1:4]),sampleTab$group)
  colnames(mePerc)<-colNam
  rownames(mePerc)<-mePerc$chrPos
  

# summary info of DMC only
  inFile<-list.files(pattern = "_table.csv")
  ctxt<-c("CG", "CHG", "CHH")
  # if(RE[r]=="AseI"){
  #   ctxt<-c("CHG", "CHH")
  # } else {
  #   ctxt<-c("CG", "CHG", "CHH")
  # }
  for(c in 1:length(ctxt)){
    temp<-read.csv(inFile, header=TRUE, sep="\t")
    temp<-unite(temp, chrPos, c(chr, pos), sep="_", remove=FALSE)
    toKeepRE<-RE[r]
    temp<-dplyr::filter(temp, RE==toKeepRE)
    temp<-dplyr::filter(temp, factor=="Treat")
    temp<-dplyr::filter(temp, context==ctxt[c])
    #temp<-dplyr::filter(temp, feature!="gene")
    temp<-dplyr::filter(temp, fdrs<=0.05)
    sort(temp$chrPos)
    rownames(temp)<-temp$chrPos

# subset meth data that has DMC
    dfDMC<-mePerc[rownames(temp),]
    #if(rownames(dfDMC)[10]==rownames(temp)[10]){
    if(sum(rownames(dfDMC)==rownames(temp))==nrow(temp)){
     cat("Same order in temp and DMC data!\n")
      cat("Adding feature column to data.\n")
      dfDMC$ftre<-temp$feature
    }
#reorder table
    DMCmelt<-reshape2::melt(dfDMC, id=c("ftre", "chrPos","chr", "pos","context"), variable="Samples")
    DMCmelt$Samples<-gsub("[[:digit:]]","", DMCmelt$Samples)
    write.table(DMCmelt, paste0(baseDir,"/results/",RE[r],"_methyltation_",ctxt[c],"_DMC.csv"),row.names = FALSE, sep="\t", col.names=T, quote=FALSE)
    ## cal dif
    genFeature<-unique(DMCmelt$ftre)
    for(f in 1:length(genFeature)){
      tempDMC<-dplyr::filter(DMCmelt, ftre==genFeature[f])
      glm.1<-glm(value ~ Samples, family = poisson, data=tempDMC)
      #Type III Analysis of Variance Table with Satterthwaite's method
      fit<-Anova(glm.1) #from car package
      fit  # revisar, si no hay interaccion hacer el analisis para cada efecto
      #Degrees-of-freedom method: kenward-roger
      #P value adjustment: tukey method for comparing a family of 9 estimates
      #another option: adjust="bonferroni"
      tukey.treat<-emmeans(glm.1, list(pairwise ~ Samples), adjust = "tukey")
      sink(paste0(baseDir,"/results/",RE[r],"_DMC_",ctxt[c],"_",genFeature[f],"_Report.txt"))
      cat("#########   ANOVA analysis for each variable  #######\n")
      print(fit)
      cat("#########   Variable Names  #######\n")
      print(tukey.treat)
      sink()
    }
    #calculate mean and sd for each factor and feature
    summDMC<-plyr::ddply(DMCmelt, c("ftre", "Samples"), summarise,
                     mean = mean(value,na.rm=TRUE),
                     events=length(value),
                     #sd = sd(value,na.rm=TRUE))
                     sem = sd(value,na.rm=TRUE)/sqrt(length(value)))
#save summary
   write.table(summDMC, paste0(baseDir,"/results/",RE[r],"_methLeveles4_",ctxt[c],"_DMC.csv"),row.names = FALSE, sep="\t", col.names=T, quote=FALSE)
#order x label
    summDMC$ftre <- factor(summDMC$ftre,levels = c("gene","transposon", "repeat", "nothing"))
# Plot from subset
   ggplot(summDMC,aes(x=ftre,y=mean, fill=Samples))+ geom_bar(stat="identity", position=position_dodge()) +
   geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2, position=position_dodge(.9)) +
   scale_y_continuous(limits = c(0, 100)) +
   scale_fill_manual(values=c('cyan4','hotpink3')) +
   labs(title=ctxt[c], x="", y = "Methylation (%)")+ theme_bw()
   ggsave(paste0(baseDir,"/plots/",RE[r],"_DMC4",ctxt[c],"_context.pdf"))
  }
}