#!/usr/bin/env Rscript
rm(list=ls())
sequenceFile <- "<PATH_TO_FILE>/consensus_cluster.renamed.fa"
geneFile <- "<PATH_TO_FILE>/diamond_output.gz"
repeatmaskerfile <- "<PATH_TO_FILE>/repeatMasker_output.gz"
outfileName <- "<PATH_TO_FILE>/mergedAnnot.csv"
source(file.path("<PATH_TO_FILE>/", "commonFunctions.R"))

################################################################################################
### load data, first get all possible sequence IDs from the fasta
allPossibleIDs <- f.extraxt.fasta.IDs(sequenceFile) # see commonFunctions
allPossibleIDs <- paste0("chr", allPossibleIDs)
seqAnno <- f.load.seq.annotation(geneFile, repeatmaskerfile, allPossibleIDs) # see commonFunctions
write.csv(seqAnno, outfileName)

sink("Annotation_summary.txt")
cat("Loaded", length(allPossibleIDs), "sequences.\n")
cat("Annotated as gene:", sum(seqAnno$gene == "yes"), "\n")
cat("Annotated as repeat:", sum(seqAnno$repMasVerySimple == "repeat"), "\n")
cat("Annotated as transposon:", sum(seqAnno$repMasVerySimple == "transposon"), "\n")
cat("Unannotated:", sum((seqAnno$gene == "no")&(seqAnno$repMasVerySimple=="nothing")), "\n")
sink()