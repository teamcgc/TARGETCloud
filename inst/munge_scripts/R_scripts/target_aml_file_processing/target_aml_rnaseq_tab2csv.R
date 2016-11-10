#
#   Merge all AML RNASeq data into one table for google cloud storage.
#   
#   __________________________________________________________________
#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


rm(list=ls(all=TRUE))

setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_aml/", 
        "target_aml_gene_rnaseq"))

sampleInfo <- read.table("TARGET_AML_mRNA-seq_20160826.sdrf.txt", 
        header=TRUE, sep="\t", quote="")
sampleInfo <- sampleInfo[, c(1,25)];
sampleInfo <- sampleInfo[-which(sampleInfo[,2] == ""),]
sampleInfo <- sampleInfo[-which(duplicated(sampleInfo[,2])),]

dim(sampleInfo)                                 #   [1] 66  2

file.pattern <- ".gene.quantification.txt"
rnaseq.files <- list.files(pattern=file.pattern)
length(rnaseq.files)                            #   [1] 66
system("wc -l AML*gene.quantification.txt")     #   3424608 total

aFile <- 1;
message(aFile, ": ", rnaseq.files[aFile]);
aml_rnaseq <- read.table(rnaseq.files[aFile], header=TRUE, sep="\t", quote="")
dim(aml_rnaseq)    #   [1]  51887     3

colnames(aml_rnaseq)    #   [1] "gene"      "raw_count" "rpkm" 
head(aml_rnaseq)
#        gene raw_count      rpkm
#   1     7SK       119 1.1883301
#   2    A1BG       215 0.0578736
#   3    A1CF         0 0.0000000
#   4   A2LD1         8 0.0770183
#   5     A2M        21 0.2397502
#   6 A2M-AS1         0 0.0000000

sampleName <- sub(file.pattern, "", rnaseq.files[aFile])
sampleName <- sub("\\.", "_", sampleName)
nameRow <- which(sampleInfo[,2] == sampleName)
sampleID <- as.character(sampleInfo[nameRow,1])
aml_rnaseq <- data.frame(sampleID=rep(sampleID, nrow(aml_rnaseq)),
        alternateID=rep(sampleName, nrow(aml_rnaseq)), aml_rnaseq)

for(aFile in 2:length(rnaseq.files))
{
    message(aFile, ": ", rnaseq.files[aFile]);
    moreData <- read.table(rnaseq.files[aFile], header=TRUE, sep="\t", quote="")

    sampleName <- sub(file.pattern, "", rnaseq.files[aFile])
    sampleName <- sub("\\.", "_", sampleName)
    nameRow <- which(sampleInfo[,2] == sampleName)
    sampleID <- as.character(sampleInfo[nameRow,1])

    moreData <- data.frame(sampleID=rep(sampleID, nrow(moreData)),
        alternateID=rep(sampleName, nrow(moreData)), moreData)
    aml_rnaseq <- rbind(aml_rnaseq, moreData);
}


#   Save output to files
#   ===========================================================================
save(aml_rnaseq, file="target_aml_RNAseq_gene_level_3_all_samples.RData")
write.table(aml_rnaseq, sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE,
    file="target_aml_RNAseq_gene_level_3_all_samples.csv" )
