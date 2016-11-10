#
#   Covert gene expression data from Affymetrix array platform to a
#   csv format for Google cloud storage.
#
#   Input data may be a table with rows for genes and columns for 
#   samples or in separated files.
#
#   Columns for BigQuery are:
#   sampleID        bar code of taget project samples  
#   probeSetID      Affymetrix probe set ID 
#   geneSymbols     HUGO gene symbols 
#   targetID        Ensembl transcript ID (ENSTxxxxxxx) 
#   expValues       log2 of signal intensity
#
#   ________________________________________________________________
#   <BigQuery table><BigQuery table><BigQuery table><BigQuery table>


rm(list=ls(all=TRUE))

#   1. AML samples. probeID in aml data is transcript_cluster_id
#   ================================================================

setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_aml/", 
        "target_aml_gene_array"))

geneExp <- read.table("TARGET_AML_GE_level3.txt", header=TRUE, 
                    sep="\t", quote="")
NL.BM.columns <- grep("NL.BM", colnames(geneExp))
geneExp <- geneExp[,-NL.BM.columns]

   totalRows <- nrow(geneExp)           #   [1] 33297
totalSamples <- ncol(geneExp) - 2;      #   [1] 261

geneExp <- geneExp[order(as.character(geneExp[,1])), ]


load(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_aml//",
    "target_aml_gene_array/HuGene.1.1.st.v1.na36.hg19.gene.RData"))
    
redundents <- which(duplicated(HuGene.1.1.st.v1.na36.hg19.gene$transcriptID))
annotation <- HuGene.1.1.st.v1.na36.hg19.gene[-redundents,]
annotation <- annotation[order(as.character(annotation[,2])), ]

dim(annotation)     #   [1] 33297     4
dim(geneExp)        #   [1] 33297   263
sum(as.character(annotation[,2]) == as.character(geneExp[,1])) #    [1] 33297

#   Reformat the sample IDs
#
sampleNames <- colnames(geneExp)[-c(1,2)]
sampleNames <- gsub("\\.", "-", sampleNames);
sampleID <- rep(sampleNames, each=totalRows)

probeSetID  <- rep(as.character(geneExp[,1]), times=totalSamples)
geneSymbols <- rep(as.character(geneExp[,2]), times=totalSamples)
targetID    <- rep(as.character(annotation$targetID), times=totalSamples)
exprValues  <- as.numeric(as.matrix(geneExp[, -c(1:2)]))

target_aml_affy_exp <- data.frame(sampleID, probeSetID, geneSymbols,
                                targetID, exprValues)
dim(target_aml_affy_exp)    #   [1] 8690517       5

save(target_aml_affy_exp, file="target_aml_affy_exp.RData")
write.table(target_aml_affy_exp, file="target_aml_affy_exp.csv",
        sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)


#   data check
#   =================================
length(sampleID)
length(probeSetID)
length(targetID)
length(geneSymbols)
length(exprValues)

unique(sampleID[1:length(sampleID)])

for(x in 1:totalSamples)
{
    first <- (x-1) * totalRows + 1;
    last  <- x*totalRows;
    columns <- x + 2;

    print(unique(sampleID[first:last]));
    same <- sum(probeSetID[first:last] == geneExp[,1]);
    message("    probe set ID same: ", same)

    same <- sum(geneSymbols[first:last] == geneExp[,2]);
    message("    geneSymbols same: ", same)

    same <- sum(exprValues[first:last] == geneExp[,columns])
    message("    Expr values same: ", same)
}



#   2.  CCSK samples. Input data is saved as gct format
#   ==============================================================

rm(list=ls(all=TRUE))

setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_ccsk/", 
        "target_ccsk_gene_expression"))

geneExp <- read.table("TARGET_CCSK_GEX_L3_NotCollapsed_20130522.gct", 
                    skip=2, header=TRUE, sep="\t", quote="")

geneExp[,2] <- gsub(".{1,}, ", "", geneExp[,2]);
geneExp[,2] <- gsub("\"", "", geneExp[,2])
geneExp <- geneExp[order(as.character(geneExp[,1])), ]

   totalRows <- nrow(geneExp)
totalSamples <- ncol(geneExp) - 2;

annotation.file <- paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
    "annotations/HG-U133_Plus_2.na36.annot.csv")
annotation <- read.csv(annotation.file, header=TRUE, sep=",", quote="\"",
                            comment.char="#")

targetID <- rep("", nrow(annotation))
Ensembl <- as.character(annotation$Ensembl);
for(aID in 1:length(Ensembl))
{
    ensembl.gene <- gsub("/", "", Ensembl[aID]);
    ensembl.gene <- unique(unlist(strsplit(ensembl.gene, "  ")));
    ENSG <- grep("ENSG", ensembl.gene); 
    if(length(ENSG) == 0) next;
    
    if(length(ensembl.gene) == 1 ) {
        targetID[aID] <- ensembl.gene;
    } else { 
        if(length(ENSG) == 1) {
            targetID[aID] <- ensembl.gene[ENSG];
        } else {
            targetID[aID] <- paste(ensembl.gene[ENSG], collapse=";");
        }
    }
} 

#   Reformat the sample IDs
#
sampleNames <- colnames(geneExp)[-c(1,2)]
sampleNames <- gsub("\\.", "-", sampleNames);
sampleID <- rep(sampleNames, each=totalRows)

probeSetID  <- rep(as.character(geneExp[,1]), times=totalSamples)
geneSymbols <- rep(as.character(geneExp[,2]), times=totalSamples)
   targetID <- rep(targetID, times=totalSamples)
exprValues  <- as.numeric(as.matrix(geneExp[,-c(1:2)]))

target_ccsk_affy_exp <- data.frame(sampleID, probeSetID, geneSymbols,
                                targetID, exprValues)
dim(target_ccsk_affy_exp)    #   [1] 710775       5

save(target_ccsk_affy_exp, file="target_ccsk_affy_exp.RData")
write.table(target_ccsk_affy_exp, file="target_ccsk_affy_exp.csv",
        sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
rm(list=ls(all=TRUE))

#


#   3.  NBL samples. Data are saved in separated files by sample
#       at two genomic levels (gene and transcript)
#   ============================================================


#   Gene expression in "Full" scope
#   ===========================================================
#
setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_nbl/", 
        "target_nbl_gene_array/gene/Full"))
expr.files <- list.files(pattern="TARGET")

aFile <- 1
geneExp <- read.table(expr.files[aFile], header=TRUE, sep="\t", quote="");
sampleID <- gsub("\\.", "-", colnames(geneExp)[2]);
target_nbl_affy_exp <- data.frame(
        sampleID=rep(sampleID, nrow(geneExp)),
        probeSetID=rep("", nrow(geneExp)),
        geneSymbols=as.character(geneExp[,1]),
        targetID=rep("", nrow(geneExp)),
        exprValues=as.numeric(geneExp[,2]));

for(aFile in 2:length(expr.files))
{
    print(expr.files[aFile]);
    geneExp <- read.table(expr.files[aFile], header=TRUE, sep="\t", quote="");
    sampleID <- gsub("\\.", "-", colnames(geneExp)[2]);
    more.data <- data.frame(
        sampleID=rep(sampleID, nrow(geneExp)),
        probeSetID=rep("", nrow(geneExp)),
        geneSymbols=as.character(geneExp[,1]),
        targetID=rep("", nrow(geneExp)),
        exprValues=as.numeric(geneExp[,2]));

    target_nbl_affy_exp <- rbind(target_nbl_affy_exp, more.data)
}

save(target_nbl_affy_exp, file="target_nbl_gene_exp_affy_full.RData")
write.table(target_nbl_affy_exp, sep=",", quote=FALSE, 
        row.names=FALSE, col.names=TRUE,
        file="target_nbl_gene_exp_affy_full_level_3.csv")


#   Transcription expression in "Full" scope
#   ===========================================================
#
rm(list=ls(all=TRUE))
setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_nbl/", 
        "target_nbl_gene_array/transcript/Full"))
expr.files <- list.files(pattern="TARGET")

#   Due to the large number of files it is bettwe to generate
#   csv file for each sample then concatenate them with shell
#   scripts

for(aFile in 1:length(expr.files))
{
    message(aFile, ": ", expr.files[aFile]);
    geneExp <- read.table(expr.files[aFile], header=TRUE, sep="\t", quote="");
    sampleID <- gsub("\\.", "-", colnames(geneExp)[2]);
    more.data <- data.frame(
        sampleID=rep(sampleID, nrow(geneExp)),
        probeSetID=as.character(geneExp[,1]),
        geneSymbols=rep("", nrow(geneExp)),
        targetID=rep("", nrow(geneExp)),
        exprValues=as.numeric(geneExp[,2]));

    out.file <- sub("txt", "csv", expr.files[aFile]);
    write.table(more.data, file=out.file, sep=",", quote=FALSE,
                row.names=FALSE, col.names=TRUE)
}



aFile <- 1
geneExp <- read.table(expr.files[aFile], header=TRUE, sep="\t", quote="");
sampleID <- gsub("\\.", "-", colnames(geneExp)[2]);
target_nbl_affy_exp <- data.frame(
        sampleID=rep(sampleID, nrow(geneExp)),
        probeSetID=as.character(geneExp[,1]),
        geneSymbols=rep("", nrow(geneExp)),
        targetID=rep("", nrow(geneExp)),
        exprValues=as.numeric(geneExp[,2]));

for(aFile in 2:length(expr.files))
{
    message(aFile, ": ", expr.files[aFile]);
    geneExp <- read.table(expr.files[aFile], header=TRUE, sep="\t", quote="");
    sampleID <- gsub("\\.", "-", colnames(geneExp)[2]);
    more.data <- data.frame(
        sampleID=rep(sampleID, nrow(geneExp)),
        probeSetID=as.character(geneExp[,1]),
        geneSymbols=rep("", nrow(geneExp)),
        targetID=rep("", nrow(geneExp)),
        exprValues=as.numeric(geneExp[,2]));

    target_nbl_affy_exp <- rbind(target_nbl_affy_exp, more.data)
}

save(target_nbl_affy_exp, file="target_nbl_transcript_exp_affy_full.RData")
write.table(target_nbl_affy_exp, sep=",", quote=FALSE, 
        row.names=FALSE, col.names=TRUE,
        file="target_nbl_transcript_exp_affy_full_level_3.csv")
rm(list=ls(all=TRUE))



#   4.  WT samples. Input data is saved as gct format
#   ==============================================================
rm(list=ls(all=TRUE))
setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_wt/", 
        "target_wt_gene_array"))

geneExp <- read.table("TARGET_WT_GEX_L3_NotCollapsed_20130912.gct", 
                    skip=2, header=TRUE, sep="\t", quote="")

geneExp[,2] <- gsub(".{1,}, ", "", geneExp[,2]);
geneExp[,2] <- gsub("\"", "", geneExp[,2])

totalRows <- nrow(geneExp)
totalSamples <- ncol(geneExp) - 2;

#   Reformat the sample IDs
#
sampleNames <- colnames(geneExp)[-c(1,2)]
sampleNames <- gsub("\\.", "-", sampleNames);
sampleID <- rep(sampleNames, each=totalRows)

probeSetID <- as.character(geneExp[,1]);
probeSetID <- rep(probeSetID, times=totalSamples)

geneSymbols <- as.character(geneExp[,2])
geneSymbols <- rep(geneSymbols, times=totalSamples)

targetID <- rep("", times=totalRows*totalSamples)
exprValues <- as.numeric(as.matrix(geneExp[,-c(1:2)]))

target_wt_affy_exp <- data.frame(sampleID, probeSetID, geneSymbols,
                                targetID, exprValues)
dim(target_wt_affy_exp)    #   [1]  6998400       5

save(target_wt_affy_exp, file="target_wt_affy_exp.RData")
write.table(target_wt_affy_exp, file="target_wt_affy_exp.csv",
        sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)
rm(list=ls(all=TRUE))



