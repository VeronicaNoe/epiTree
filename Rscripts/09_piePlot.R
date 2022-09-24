## Set path to input files:
rm(list=ls())
wd<-getwd()
baseDir <- gsub("/results", "", wd)
scriptDir <- file.path(baseDir, "scripts")

## Load packages silently:
suppressPackageStartupMessages({
  library(ggplot2)
  library(webr)
  library(dplyr)
})

# Load data
tbl<-paste0(baseDir,"/rawData/featurebycontext.csv")
# process both datasets
data<-read.csv(tbl, header=TRUE, check.names = FALSE, sep="\t")
re<-unique(data$RE)
ctxt<-unique(data$Context)
head(data)
for(r in 1:length(re)){
  df<-subset(data,RE==re[r])
  for(c in 1:length(ctxt)){
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
