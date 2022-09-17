rm(list=ls())
wd<-getwd()
getwd()
baseDir <- gsub("/results", "", wd)
scriptDir <- file.path(baseDir, "scripts")

## load packages silently
library(ggplot2)
library(webr)
library(dplyr)

# Building a table with the data for the plot
tbl<-paste0(baseDir,"/rawData/featurebycontext.csv")
data<-read.csv(tbl, header=TRUE, check.names = FALSE, sep="\t")
re<-unique(data$RE)
ctxt<-unique(data$Context)
head(data)
for(r in 1:length(re)){
#for(r in 1:1){
  df<-subset(data,RE==re[r])
  for(c in 1:length(ctxt)){
  #for(c in 1:1){
    subdata<-subset(df,Context==ctxt[c])
    PD<-subdata %>% group_by(Feature) %>% summarise(n = sum(Total_Cs))
    print(PD)
    p<-PieDonut(PD, aes(Feature, count=n), 
             ratioByGroup = FALSE, showPieName = FALSE,
             explode = 2, addPieLabel = FALSE)
    pdf(paste0(baseDir,"/plots/",re[r],"_","FeaturesSequenced_",ctxt[c], ".pdf"))
    plot(p, cex = 0.6, hang=-0.5,main=paste0(re[r],"_",ctxt[c]))
    dev.off()
  }
}

tbl<-paste0(baseDir,"/rawData/featurebycontext_DMC.csv")
data<-read.csv(tbl, header=TRUE, check.names = FALSE, sep="\t")
re<-unique(data$RE)
ctxt<-unique(data$Context)
head(data)
for(r in 1:length(re)){
#for(r in 1:1){
  df<-subset(data,RE==re[r])
  for(c in 1:length(ctxt)){
  #for(c in 1:1){
    subdata<-subset(df,Context==ctxt[c])
    PD<-subdata %>% group_by(Feature, Factor) %>% summarise(n = sum(Cs))
    print(PD)
    pdf(paste0(baseDir,"/plots/",re[r],"_","FeaturesSequenced_DMC_",ctxt[c], ".pdf"))
    plot(PieDonut(PD, aes(Feature,Factor, count=n), 
                  ratioByGroup = FALSE, showPieName = FALSE,
                  addPieLabel = FALSE))
    dev.off()
  }
}


