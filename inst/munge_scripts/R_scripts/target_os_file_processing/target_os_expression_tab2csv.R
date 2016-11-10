

#   1.  Merge all RNAseq gene expression files 
#       as one big file in csv format
#   ====================================================================
rm(list=ls(all=TRUE))

input.files <- list.files(path="target_os_expression", 
                pattern="gene", full.names=TRUE)
length(input.files) #   [1] 86

system("wc -l target_os_expression/*gene.quantification.txt")
#   3911022 total

aFile <- 1;
message(aFile, ": ", input.files[aFile]);
bigData <- read.table(input.files[aFile], header=TRUE, sep="\t", quote="")
dim(bigData)    #   [1] 45476     5

colnames(bigData)

#   [1] "ensembl_gene_id" "mean_length"     "mean_eff_length" 
#   [4] "est_counts"      "tpm"


#   for data check
#   ===============================================
#
totalRows <- nrow(bigData)
totalCols <- ncol(bigData)

geneID        <- as.character(bigData$ensembl_gene_id)
meanLength    <- as.numeric(bigData$mean_length)

sampleID <- sub(".{1,}/", "", input.files[aFile]);
sampleID <- sub(".gene.quantification.txt", "", sampleID);
bigData <- cbind(sampleID, bigData);
head(bigData)

for(aFile in seq_len(length(input.files))[-1])
{
    message(aFile, ": ", input.files[aFile]);
    moreData <- read.table(input.files[aFile], header=TRUE, sep="\t", quote="")

    if(nrow(moreData) != totalRows || ncol(moreData) != totalCols)
        stop("Bad file dimensions!\n")

    #   New sampleID 
    #
    moreSampleID <- sub(".{1,}/", "", input.files[aFile]);
    moreSampleID <- sub(".gene.quantification.txt", "", moreSampleID);
    if(moreSampleID == sampleID) stop("Sample ID error!")

    #   geneID in all files should be a same set
    #
    moreGeneID  <- as.character(moreData$ensembl_gene_id)
    if(sum(moreGeneID == geneID) != totalRows) 
        stop("ensembl_gene_id error!")

    #   mean_length in all files should be a same set
    #
    moreMeanLength <- as.character(moreData$mean_length)
    if(sum(moreMeanLength == meanLength) != totalRows)
        stop("mean_length error!")

    moreData <- cbind(sampleID=moreSampleID, moreData);
    print(head(moreData));
    
    bigData <- rbind(bigData, moreData);
}

#   data check
#   =======================================================
dim(bigData)                        #   [1] 3910936       6
totalRows*length(input.files)       #   [1] 3910936
length(unique(bigData$sampleID))    #   [1] 86
head(bigData)
tail(bigData)


#   Save output to files
#   ===========================================================================
save(bigData, file="target_os_RNAseq_gene_expression_level_3_all_samples.RData")
write.table(bigData, sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE,
    file="target_os_RNAseq_gene_expression_level_3_all_samples.csv" )


#   data check
#   ===============================================================
system("head -n 10 target_os_RNAseq_gene_expression_level_3_all_samples.csv")
system("tail -n 10 target_os_RNAseq_gene_expression_level_3_all_samples.csv")



#   2.  Merge all RNAseq isoform expression files 
#       as one big file in csv format
#   ====================================================================

rm(list=ls(all=TRUE))

input.files <- list.files(path="target_os_expression", 
                pattern="isoform", full.names=TRUE)
length(input.files) #   [1] 86

system("wc -l target_os_expression/*isoform.quantification.txt")
#   16509248 total

aFile <- 1;
message(aFile, ": ", input.files[aFile]);
bigData <- read.table(input.files[aFile], header=TRUE, sep="\t", quote="")
dim(bigData)    #   [1] 191967      5

colnames(bigData)

#   [1] "target_id"  "length"     "eff_length" "est_counts" "tpm"



#   for data check
#   ===============================================
#
totalRows <- nrow(bigData)
totalCols <- ncol(bigData)

targetID     <- as.character(bigData$target_id)
targetLength <- as.numeric(bigData$length)

sampleID <- sub(".{1,}/", "", input.files[aFile]);
sampleID <- sub(".isoform.quantification.txt", "", sampleID);

bigData <- cbind(sampleID, bigData);
head(bigData)

for(aFile in seq_len(length(input.files))[-1])
{
    message(aFile, ": ", input.files[aFile]);
    moreData <- read.table(input.files[aFile], header=TRUE, sep="\t", quote="")

    if(nrow(moreData) != totalRows || ncol(moreData) != totalCols)
        stop("Bad file dimensions!\n")

    #   New sampleID 
    #
    moreSampleID <- sub(".{1,}/", "", input.files[aFile]);
    moreSampleID <- sub(".isoform.quantification.txt", "", moreSampleID);
    if(moreSampleID == sampleID) stop("Sample ID error!")

    #   geneID in all files should be a same set
    #
    moreTargetID  <- as.character(moreData$target_id)
    if(sum(moreTargetID == targetID) != totalRows) 
        stop("target_id error!")

    #   mean_length in all files should be a same set
    #
    moreLength <- as.character(moreData$length)
    if(sum(moreLength == targetLength) != totalRows)
        stop("length error!")

    moreData <- cbind(sampleID=moreSampleID, moreData);
    print(head(moreData));
    
    bigData <- rbind(bigData, moreData);
}


#   data check
#   ========================================================
dim(bigData)                        #   [1] 16509162       6
totalRows*length(input.files)       #   [1] 16509162
length(unique(bigData$sampleID))    #   [1] 86
head(bigData)



#   Save output to files
#   ===========================================================================
outFile <- "target_os_RNAseq_isoform_expression_level_3_all_samples.RData";
save(bigData, file=outFile);
write.table(bigData, sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE,
    file="target_os_RNAseq_isoform_expression_level_3_all_samples.csv");


#   data check
#   ===============================================================
system("head -n 10 target_os_RNAseq_isoform_expression_level_3_all_samples.csv")
system("tail -n 10 target_os_RNAseq_isoform_expression_level_3_all_samples.csv")














