#'@title a function for distances in methylation file.
#'@author Verónica Noé Ibañez \email{veronicanoeibanez@@gmail.com}.
#'@export
f.meth.distances <- function (methCov, totCov) {
  percMe <- methCov/totCov
  binTotCov <- !is.na(totCov)
  diffSum <- matrix(NA, nrow = length(commonSamples), ncol = length(commonSamples), dimnames = list(commonSamples, commonSamples))
  matchCount <- matrix(NA, nrow = length(commonSamples), ncol = length(commonSamples), dimnames = list(commonSamples, commonSamples))
  for (i in 1:(length(commonSamples)-1)) {
    sampleA <- commonSamples[i]
    for (j in (i+1):length(commonSamples)) {
      sampleB <- commonSamples[j]
      commonMask <- rowSums(binTotCov[,c(sampleA, sampleB)]) == 2
      if (sum(commonMask) == 0) { next }
      forDist <- percMe[commonMask, c(sampleA, sampleB)]
      if (is.vector(forDist)) {
        matchCount[sampleA, sampleB] <- 1
        diffSum[sampleA, sampleB] <- abs(forDist[sampleA]-forDist[sampleB])
      } else {
        matchCount[sampleA, sampleB] <- nrow(forDist)
        curDiffs <- forDist[,sampleA]-forDist[,sampleB]
        diffSum[sampleA, sampleB] <- sum(abs(curDiffs))
      }
    }
  }
  methDist <- diffSum/matchCount
  return(methDist)
}

#'@title a function for printing
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

#' a color gradient that I like
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.pinkredyellowblueblackNICE <- function (x) {
  #r <- approx(c(0, 0.25, 0.5, 0.75, 1), c(154/255, 205/255, 238/255, 24/255, 0), n = x)$y
  r <- approx(c(0, 0.25, 0.5, 0.75, 1), c(117/255, 239/255, 238/255, 24/255, 0), n = x)$y
  #g <- approx(c(0, 0.25, 0.5, 0.75, 1), c(50/255, 38/255, 201/255, 116/255, 0), n = x)$y
  g <- approx(c(0, 0.25, 0.5, 0.75, 1), c(50/255, 38/255, 208/255, 118/255, 0), n = x)$y
  #b <- approx(c(0, 0.25, 0.5, 0.75, 1), c(205/255, 38/255, 0/255, 205/255, 0), n = x)$y
  b <- approx(c(0, 0.25, 0.5, 0.75, 1), c(156/255, 38/255, 0/255, 156/255, 0), n = x)$y
  return(rgb(r, g, b))
}

#' a color gradient that I like
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.blackblueyellowredpinkNICE <- function(x) {
  return(rev(f.pinkredyellowblueblackNICE(x)))
}

#' add transparancy to an existing color
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
f.add.alpha <- function(col, alpha=1){
  if(missing(col)){stop("Please provide a vector of colours.")}
  apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))
}

#'@title summarize a table according to groups and names
#'@param x <dataframe, matrix>: numeric, colnames as in byTab
#'@param byTab <dataframe>: at least two columns: sample and group
#'@param summaryFunction: a function like mean, median, sum
#'@return the summarized table (a matrix)
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.summarize.columns <- function(x, byTab, summaryFunction) {
  byTab$sample <- as.character(byTab$sample)
  byTab$group <- as.character(byTab$group)
  groupnames <- unique(byTab$group)
  out <- matrix(0, nrow = nrow(x), ncol = length(groupnames), dimnames = list(rownames(x), groupnames))
  for (gn in groupnames) {
    toSummarize <- byTab$sample[byTab$group == gn]
    if (length(toSummarize) > 1) {out[,gn] <- apply(x[,toSummarize], 1, summaryFunction)}
    else { out[,gn] <- x[,toSummarize] }
  }
  return(out)
}

#'@title summarize a table according to groups and names
#'@param x <dataframe, matrix>: numeric, colnames as in byTab
#'@return the summarized table (a matrix)
#'@author Verónica Noé Ibañez \email{veronicanoeibanez@@gmail.com}.
#'@export
f.summarize.badSamples <- function(x) {
  methCovStacked<-reshape2::melt(x, variable.name="Samples",value.name="ReadCov")
  summMethCov<-ddply(methCovStacked, c("Samples"), summarise, Sum=sum(ReadCov, na.rm=TRUE)) ## sum all reads per sample
  summMethCov$Samples <- factor(summMethCov$Samples, levels=unique(summMethCov$Samples))
  return(summMethCov)
}


#'@title simplify detailed repeatmasker annotations
#'@param seqAnno annotation data (see \code{\link{f.load.seq.annotation}})
#'@return the simplified annotation data
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.simplify.repeat.masker <- function(seqAnno) {
  featureSimplification <- list(
    DNA_Ginger = c("DNA_Ginger"),
    DNA_hAT = c("DNA_hAT", "DNA_hAT_Ac", "DNA_hAT_Charlie", "DNA_hAT_Tag1", "DNA_hAT_Tip100"),
    DNA_TcMar = c("DNA_TcMar_Pogo", "DNA_TcMar_Mariner", "DNA_TcMar_Tc1", "DNA_TcMar_Stowaway"),
    DNA_MuLE = c("DNA_MULE_MuDR", "DNA_MULE"),
    DNA_PIF_Harbinger = c("DNA_PIF_Harbinger"),
    DNA_CMC = c("DNA_CMC_EnSpm"),
    DNA_Sola = c("DNA_Sola_1"),
    DNA = c("DNA"),
    LTR_Cassandra = c("LTR_Cassandra"),
    LTR_Caulimovirus = c("LTR_Caulimovirus"),
    LTR_Copia = c("LTR_Copia"),
    LTR_Gypsy = c("LTR_Gypsy"),
    LTR = c("LTR"),
    RC_Helitron = c("RC_Helitron"),
    LINE_Penelope = c("LINE_Penelope"),
    LINE_CRE = c("LINE_CRE"),
    LINE_RTE = c("LINE_RTE_BovB"),
    LINE_L1 = c("LINE_L1", "LINE_L1_Tx1"),
    LINE_I = c("LINE_I"),
    SINE_tRNA = c("SINE_tRNA", "SINE_tRNA_RTE"),
    SINE = c("SINE"),
    retroposon = c("Retroposon"),
    lowComplexity = c("Low_complexity"),
    uknRepeats = c("unClassRep", "Unknown"),
    simple_repeat = c("Simple_repeat"),
    satellite_centromeric = c("Satellite_centr"),
    satellite = c("Satellite"),
    tRNA = c("tRNA"),
    rRNA = c("rRNA"),
    snRNA = c("snRNA"),
    others = c("Other_Composite", "ARTEFACT", "nothing", "unAnnotated")
  )
  extendedToSimple <- list()
  for (simple in names(featureSimplification)) {
    for (entry in featureSimplification[[simple]]){
      extendedToSimple[[entry]] <- simple
    }
  }
  extendedToSimple <- unlist(extendedToSimple)
  # check what's missing (should be added manually)
  missing <- unique(seqAnno$repMas[!(seqAnno$repMas %in% names(extendedToSimple))])
  if (length(missing) > 0) {
    f.print.message("Missing repeatmasker features (add in f.simplify.repeat.masker):")
    cat(paste0(missing, collapse = '\n'), '\n')
  }
  seqAnno$repMas <- extendedToSimple[seqAnno$repMas]
  return(seqAnno)
}

#'@title read gene, transposon, and repeat annotations
#'@param geneFile path to the file with the genes
#'@param repMasFile path to the file with the repeats/transposons
#'@param allSeqs (optional) a vector with all sequence numbers
#'@param simplifyRepMasker attempt to simplify the detailed repeatmasker annotation (see \code{\link{f.simplify.repeat.masker}})
#'@return a merged annotation
#'@note
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.load.seq.annotation <- function(geneFile, repMasFile, allSeqs = c(), simplifyRepMasker = FALSE) {
  if (grepl("\\.gz$", geneFile)) {
    geneAnno <- read.table(gzfile(geneFile), sep = '\t', stringsAsFactors = FALSE, header = FALSE, col.names = c("seqID", "refID", "score"))
  } else {
    geneAnno <- read.table(geneFile, sep = '\t', stringsAsFactors = FALSE, header = FALSE, col.names = c("seqID", "refID", "score"))
  }
  if (grepl("\\.gz$", repMasFile)) {
    repMasAnno <- read.table(gzfile(repMasFile), sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  } else {
    repMasAnno <- read.table(repMasFile, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
  }
  geneAnno$seqID <- paste0("chr", geneAnno$seqID)
  repMasAnno$queryID <- paste0("chr", repMasAnno$queryID)
  repMasAnno$repeatClassOrFamilySimple <- sapply(repMasAnno$repeatClassOrFamily, function(x) unlist(strsplit(x, "\\?|\\/"))[1])
  repMasAnno$teOrRep <- "transposon"
  repMasAnno$teOrRep[repMasAnno$repeatClassOrFamilySimple %in% c("Low_complexity", "Other", "rRNA", "Satellite", "Simple_repeat", "snRNA", "tRNA", "Unknown")] <- "repeat"
  if (length(allSeqs) == 0) { allSeqs <- union(geneAnno$seqID, repMasAnno$queryID) }
  out <- data.frame(gene = rep("no", length(allSeqs)),
                    repMas = rep("nothing", length(allSeqs)),
                    repMasSimple = rep("nothing", length(allSeqs)),
                    repMasVerySimple = rep("nothing", length(allSeqs)),
                    geneScore = rep(0, length(allSeqs)),
                    repMasScore = rep(0, length(allSeqs)),
                    stringsAsFactors = FALSE, row.names = allSeqs)
  out[geneAnno$seqID, "gene"] <- "yes"
  out[repMasAnno$queryID, "repMas"] <- repMasAnno$repeatClassOrFamily
  out[repMasAnno$queryID, "repMasSimple"] <- repMasAnno$repeatClassOrFamilySimple
  out[repMasAnno$queryID, "repMasVerySimple"] <- repMasAnno$teOrRep
  out[geneAnno$seqID, "geneScore"] <- geneAnno$score
  out[repMasAnno$queryID, "repMasScore"] <- repMasAnno$swScore
  # replace some odd characters
  out$repMas <- gsub("-", "_", out$repMas, fixed = TRUE)
  out$repMas <- gsub("/", "_", out$repMas, fixed = TRUE)
  out$repMas <- gsub("?", "", out$repMas, fixed = TRUE)
  if (simplifyRepMasker) {
    out <- f.simplify.repeat.masker(out)
  }
  return(out)
}

#'@title extract all IDs from a fasta file
#'@param geneFile path to the fasta file (can be .gz)
#'@return vector with IDs
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.extraxt.fasta.IDs <- function(infileName) {
  if (grepl("\\.gz$", infileName)) {
    infile <- gzcon(file(infileName, "rb"))
  }
  else {
    infile <- file(infileName, open = "rt")
  }
  out <- rep(NA, 1e6)
  counter <- 0
  while (length(oneLine <- readLines(infile, n = 1, warn = FALSE)) > 0) {
    if (substr(oneLine, 1, 1) != ">") { next }
    counter <- counter + 1
    seqID <- gsub(">|[[:space:]]", "", oneLine)
    if (counter > length(out)) {
      out <- c(out, rep(NA, 1e6))
    }
    out[counter] <- seqID
  }
  close(infile)
  out <- out[!is.na(out)]
  return(out)
}

#'@title read the design table and check if the sample names are ok
#'@param designTablePath path to the design table file
#'@param colNamesForGrouping a vector with the column names used for grouping, if given, will check and add a column "group" with all of them pasted together
#'@return sampleTab
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.read.sampleTable <- function(designTablePath, colNamesForGrouping = c()) {
  sampleTab <- read.table(designTablePath, sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  # avoid using anything else than underscores, letters, and numbers in sample names
  rownames(sampleTab) <- gsub("-|\\.", "_", rownames(sampleTab))
  # and don't start sample names with a number
  if (sum(grepl("^[[:digit:]]", rownames(sampleTab))) > 0) {
    rownames(sampleTab) <- paste0("sample_", rownames(sampleTab))
  }
  # check if some samples need to be removed
  if ("sampleRemovalInfo" %in% colnames(sampleTab)) {
    f.print.message("Removing", sum(sampleTab$sampleRemovalInfo == "ERASE"), "samples due to the sampleRemovalInfo column")
    sampleTab <- subset(sampleTab, sampleRemovalInfo != "ERASE")
  }
  # check for the grouping variables
  if (length(colNamesForGrouping) > 0) {
    if (sum(colNamesForGrouping %in% colnames(sampleTab)) != length(colNamesForGrouping)) {
      f.print.message("The column names for grouping don't exist!")
      cat(paste0(setdiff(colNamesForGrouping, colnames(sampleTab)), collapse = '\n'), '\n')
      f.print.message("Columns in the design table:")
      cat(paste0(colnames(designTable), collapse = '\n'), '\n')
      stop("STOPPING!")
    }
    if (length(colNamesForGrouping) > 1) {
      sampleTab$group <- apply(sampleTab[,colNamesForGrouping], 1, function(x) paste(x, collapse = '_'))
    } else {
      sampleTab$group <- sampleTab[[colNamesForGrouping]]
    }
  }
  return(sampleTab)
}

#'@title read a methylation.bed(.gz)
#'@param infileName path to the bed file
#'@param removeMulticontextPositions if TRUE, positions with more than one context are removed
#'@param percentages if TRUE, the methylated and total columns are converted into percentages
#'@return table with the methylation data
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@author edited by Verónica Noé Ibañez \email{veronicanoeibanez@@gmail.com}.
#'@export
f.load.methylation.bed <- function(infileName, removeMulticontextPositions = TRUE, percentages = FALSE) {
  require("data.table")
  if (grepl("\\.filtMETH$", infileName)) {
    myData <- fread(infileName, sep = '\t', data.table = FALSE, fill = TRUE)
  } else {
    myData <- fread(infileName, sep = '\t', data.table = FALSE, fill = TRUE, na.string=c("None")) #na.string="None"  nrows = 1000
  }
  # columns with sampleInfo
  sampleInfoColNumbers <- grep("_total$|_methylated$", colnames(myData), value = FALSE)
  # avoid using anything else than underscores, letters, and numbers in sample names
  colnames(myData)[sampleInfoColNumbers] <- gsub("-|\\.", "_", colnames(myData)[sampleInfoColNumbers])
  # and don't start sample names with a number
  if (sum(grepl("^[[:digit:]]", colnames(myData)[sampleInfoColNumbers])) > 0) {
    colnames(myData)[sampleInfoColNumbers] <- paste0("sample_", colnames(myData)[sampleInfoColNumbers])
  }
  # remove positions with multiple contexts if requested
  if (removeMulticontextPositions) {
    temp <- unique(myData[,c("chr", "pos", "context")])
    ctxtCount <- table(paste0(temp$chr, '_', temp$pos))
    multipleContexts <- names(ctxtCount)[ctxtCount>1]
    if (length(multipleContexts) > 0) {
      f.print.message("removing", length(multipleContexts), "positions with multiple contexts.")
      myData <- subset(myData, !(paste0(chr, '_', pos) %in% multipleContexts))
    }
  }
  # convert to percentages if requested
  if (!percentages) {
    out <- myData
  } else {
    out <- myData[,c("chr", "pos", "context")]
    methCols <- grep("_methylated$", colnames(myData), value = TRUE)
    totCols <- grep("_total$", colnames(myData), value = TRUE)
    mePerc <- round(myData[,methCols]/myData[,totCols]*100, 3)
    colnames(mePerc) <- gsub("_methylated$", "", colnames(mePerc))
    for (toAdd in colnames(mePerc)) {
      out[[toAdd]] <- mePerc[[toAdd]]
    }
  }
  return(out)
}


#'@title read a merged annotation file and subset for a feature if requested
#'@param infileName path to the merged annotation (see the script mergeAnnotation.R)
#'@param selectedFeature all, gene, transposon, repeat, nothing
#'@param simplifyRepMasker attempt to simplify the detailed repeatmasker annotation (see \code{\link{f.simplify.repeat.masker}})
#'@return table with the snp data
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.load.merged.annotation <- function(infileName, selectedFeature = "all", simplifyRepMasker = FALSE) {
  mergedAnno <- read.csv(infileName, row.names = 1, stringsAsFactors = FALSE, header = TRUE)
  # subset feature
  if (selectedFeature != "all") {
    if (selectedFeature == "gene") {
      out <- subset(mergedAnno, (gene == "yes" ))
    } else if (selectedFeature == "repeat") {
      out <- subset(mergedAnno, repMasVerySimple == "repeat")
    } else if (selectedFeature == "transposon") {
      out <- subset(mergedAnno, repMasVerySimple == "transposon")
    } else if (selectedFeature == "nothing") {
      out <- subset(mergedAnno, (repMasVerySimple == "nothing") & (gene == "no"))
    } else {
      cat("ERROR, feature", feature, "does not exist.\n")
    }
  } else {
    out <- mergedAnno
  }
  if (simplifyRepMasker) {
    out <- f.simplify.repeat.masker(out)
  }
  return(out)
}

#'@title remove col/rows with NA from a matrix
#'@return a distance matrix without NA
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.remove.NA.from.distance.matrix <- function(x) {
  toRemove <- which(is.na(x), arr.ind = TRUE)
  toKeep <- setdiff(rownames(x), rownames(toRemove))
  x <- x[toKeep,toKeep]
  return(x)
}


#'@param allDistanceFiles a vector with path to distance files
#'@param imputeNA TRUE/FALSE if NAs should be imputed
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.load.many.distance.files <- function(allDistanceFiles, imputeNA) {
  require("ape")
  allDistances <- list()
  for (curFile in allDistanceFiles) {
    curFileName <- basename(curFile)
    curName <- gsub("\\.csv$", "", curFileName)
    temp <- read.csv(curFile, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
    rownames(temp) <- gsub("-|\\.", "_", rownames(temp)) # avoid using anything else than underscores, letters, and numbers in sample names
    colnames(temp) <- gsub("-|\\.", "_", colnames(temp)) # avoid using anything else than underscores, letters, and numbers in sample names
    if (sum(grepl("^[[:digit:]]", colnames(temp))) > 0) {# and don't start with a number
      rownames(temp) <- paste0("sample_", rownames(temp))# and don't start with a number
      colnames(temp) <- paste0("sample_", colnames(temp))# and don't start with a number
    }# and don't start with a number
    temp <- as.matrix(temp)
    diag(temp) <- 0
    temp[lower.tri(temp)] <- t(temp)[lower.tri(t(temp))]
    if (imputeNA) {
      if (mean(is.na(temp)) > 0.2) {
        cat("Skipping data set because more than 20% are NA.\n")
        next
      }
      if (sum(is.na(temp)) > 0) {
        prevNames <- dimnames(temp)
        temp <- additive(temp)
        dimnames(temp) <- prevNames
      }
      if (sum(temp < 0) > 0) {
        cat("Detected negative distances. Removing some samples, trying again\n")
        percNeg <- rowMeans(temp<0)
        toKeep <- setdiff(rownames(temp), names(percNeg)[percNeg>0.05]) # remove the samples with more than 5 % negative distances, then set neg distances to NA and remove them
        temp <- temp[toKeep, toKeep]
        temp[temp<0] <- NA
        temp <- f.remove.NA.from.distance.matrix(temp)
        if (nrow(temp) < 4) {
          cat("Skipping table with less than 4 samples after NA removal.\n")
          next
        }
      }
    } 
    }
    allDistances[[curName]] <- temp
  return(allDistances)
}

#'@title plot an image without text
#'@param x x-coords
#'@param y y-coords
#'@param z the image value
#'@param xLabel label for x-axis
#'@param yLabel label for y-axis
#'@param mainLabel a title
#'@param useLog TRUE/FALSE if data should be log2(x+1) transformed
#'@param ... other parameters forwarded to \code{\link{image}}
#'@author Marc W. Schmid \email{contact@@mwschmid.ch}.
#'@export
f.image.without.text <- function(x, y, z, xLabel, yLabel, mainLabel, useLog = FALSE, ...) {
  xChars <- as.character(x)
  yChars <- as.character(y)
  xPos <- 1:length(xChars)
  yPos <- 1:length(yChars)
  if (useLog) {
    toPlot <- log2(z+1)
  } else {
    toPlot <- z
  }
  image(xPos, yPos, toPlot, xlab = xLabel, ylab = yLabel, main = mainLabel, yaxt = "n", xaxt = "n", ...)
  axis(1, at = xPos, labels = xChars, outer = FALSE, las = 2)
  axis(2, at = yPos, labels = yChars, outer = FALSE, las = 1)
  return(NULL)
}
