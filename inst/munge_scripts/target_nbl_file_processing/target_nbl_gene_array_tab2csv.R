#
#   Covert NBL gene expression array data to csv format. 
#   Input data are saved in separated files by sample at two genomic 
#   levels (gene and transcript) and platform is HuEx-1_0-st-v2_01
#   ________________________________________________________________
#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



#   1.  Annotation data (HuEx-1_0-st-v2.na36.hg19). Use transcript
#       annoation file to extract ensembl gene inforamtion since 
#       they are not included in probe set annotation file.
#   ===========================================================
#
source(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/",
        "R_scripts/affymetrix_huex_probe2gene.R"));

annotation.file <- paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
        "annotations/source.text.files/", 
        "HuEx-1_0-st-v2.na36.hg19.transcript.csv");

HuEx.1.0.st.v2.na36.hg19 <- affymatrix.csv.probe2gene(
        annotation.file=annotation.file, ensembl="ENST");
        
dim(HuEx.1.0.st.v2.na36.hg19)   #   [1] 337842      4

save(HuEx.1.0.st.v2.na36.hg19, file="HuEx.1.0.st.v2.na36.hg19.RData");
write.table(HuEx.1.0.st.v2.na36.hg19, sep="\t", quote=FALSE,
    col.names=FALSE, row.names=TRUE, file="HuEx.1.0.st.v2.na36.hg19.txt");


#   2.  Gene expression in "Full" scope. Input data has only 
#       gene symbols and expression values (log2)
#   ===========================================================
#
setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_nbl/", 
                                "target_nbl_gene_array/gene/Full"))
expr.files <- list.files(pattern="TARGET", path="original.text.files", full.names=TRUE)
length(expr.files)      #   [1] 249

aFile <- 1
geneExp <- read.table(expr.files[aFile], header=TRUE, sep="\t", quote="");
sampleID <- gsub("\\.", "-", colnames(geneExp)[2]);

dim(geneExp)    #   [1] 22985     2
head(geneExp)

#     probeset_id TARGET.30.PAAPFA
#   1     AADACL3          5.52132
#   2     AADACL4          5.17251
#   3       ACADM          7.47834
#   4       ACAP3          7.80589
#   5      ACOT11          6.75915
#   6       ACOT7          8.11807


#   2.1     Get transcriptID and targetID for each gene
#   ==============================================================

geneSymbols <- as.character(geneExp[,1]);
geneList <- as.character(HuEx.1.0.st.v2.na36.hg19$geneSymbols);
targetID <- rep("", length(geneSymbols));
probeID  <- rep("", length(geneSymbols));

targetList <- as.character(HuEx.1.0.st.v2.na36.hg19$targetID);
probeList  <- as.character(HuEx.1.0.st.v2.na36.hg19$transcriptID);

for(aGene in seq_along(geneSymbols))
{
        rows <- grep(geneSymbols[aGene], geneList);
        if(length(rows) == 0) next;

        print(paste0(aGene, ": ", geneSymbols[aGene]));
        targetID[aGene] <- targetList[rows[1]];
        probeID[aGene] <- probeList[rows[1]];
}


#       2.2     Merge all samples into one table
#       ================================================================

#       The first sample
#       =================================================
target_nbl_affy_exp <- data.frame(
        sampleID=rep(sampleID, nrow(geneExp)),
        probeSetID=probeID,
        geneSymbols=geneSymbols,
        targetID=targetID,
        exprValues=as.numeric(geneExp[,2]));
out.file <- paste0("target_nbl_", sampleID, "_affy_gene_expr_full.csv")
write.table(target_nbl_affy_exp, file=out.file, sep=",", quote=FALSE,
                row.names=FALSE, col.names=TRUE);

#       All other samples
#       =================================================
for(aFile in 2:length(expr.files))
{
        print(paste0(aFile, ": ", expr.files[aFile]));
        geneExp <- read.table(expr.files[aFile], header=TRUE, sep="\t", quote="");
        if(sum(as.character(geneExp[,1]) == geneSymbols) != length(geneSymbols))
                stop("Gene symbols do not match.")

        sampleID <- gsub("\\.", "-", colnames(geneExp)[2]);
        more.data <- data.frame(
                sampleID=rep(sampleID, nrow(geneExp)),
                probeSetID=probeID,
                geneSymbols=as.character(geneExp[,1]),
                targetID=targetID,
                exprValues=as.numeric(geneExp[,2]));

        out.file <- paste0("target_nbl_", sampleID, "_affy_gene_expr_full.csv")
        write.table(more.data, file=out.file, sep=",", quote=FALSE,
                                row.names=FALSE, col.names=TRUE);
        target_nbl_affy_exp <- rbind(target_nbl_affy_exp, more.data);
}

dim(target_nbl_affy_exp)    #   [1] 5723265       5

save(target_nbl_affy_exp, file="target_nbl_affy_exp_gene_full_all_samples.RData")
write.table(target_nbl_affy_exp, sep=",", quote=FALSE, 
                        row.names=FALSE, col.names=TRUE,
                        file="target_nbl_affy_exp_gene_full_level_3_all_samples.csv")


#       3.  Transcription expression in "Full" scope
#       ===========================================================
#

setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_nbl/", 
                        "target_nbl_gene_array/transcript/Full"))
expr.files <- list.files(path="target_nbl_gene_array_transcript_full_txt_format",
                                pattern="TARGET",  full.names=TRUE)
length(expr.files)      #       [1] 249

aFile <- 1
transcriptExp <- read.table(expr.files[aFile], header=TRUE, sep="\t", quote="");
transcriptExp <- transcriptExp[order(as.character(transcriptExp$probeset_id)),]
sampleID <- gsub("\\.", "-", colnames(transcriptExp)[2]);

dim(transcriptExp)              #   [1] 273551     2
head(transcriptExp)

#     probeset_id TARGET.30.PAAPFA
#   1     2315100          5.81806
#   2     2315106          6.99158
#   3     2315109          7.26259
#   4     2315113          4.57804
#   5     2315115          6.04346
#   6     2315117          8.29391


#       3.1     Get geneSymbol and targetID for each transcript
#       ==============================================================

probeID <- as.character(transcriptExp$probeset_id)
annotID <- as.character(HuEx.1.0.st.v2.na36.hg19$transcriptID)

length(which(annotID %in% probeID))   #   [1] 273551

rows <- which(annotID %in% probeID)
transcriptAnnot <- HuEx.1.0.st.v2.na36.hg19[rows,];
the.order <- order(as.character(transcriptAnnot$transcriptID))
transcriptAnnot <- transcriptAnnot[the.order, ]

sum(probeID == as.character(transcriptAnnot$transcriptID))      #   [1] 273551
colnames(transcriptAnnot)
#      [1] "probeSetID"   "transcriptID" "geneSymbols"  "targetID"



#       3.2     Merge all samples into one table
#       ================================================================

#       The first sample
#       =================================================
target_nbl_affy_exp <- data.frame(
        sampleID=rep(sampleID, nrow(transcriptExp)),
        probeSetID=probeID,
        geneSymbols=transcriptAnnot$geneSymbols,
        targetID=transcriptAnnot$targetID,
        exprValues=as.numeric(transcriptExp[,2]));

out.file <- paste0("target_nbl_", sampleID, "_affy_exon_expr_full.csv")
write.table(target_nbl_affy_exp, file=out.file, sep=",", quote=FALSE,
                row.names=FALSE, col.names=TRUE);

#       All other samples
#       =================================================
for(aFile in 2:length(expr.files))
{
        print(paste0(aFile, ": ", expr.files[aFile]));
        transcriptExp <- read.table(expr.files[aFile], header=TRUE, 
                                                sep="\t", quote="");
        the.order <- order(as.character(transcriptExp$probeset_id))
        transcriptExp <- transcriptExp[the.order,];
    
        if(sum(as.character(transcriptExp[,1]) == probeID) != length(probeID))
                stop("probeSetID do not match.")

        sampleID <- gsub("\\.", "-", colnames(transcriptExp)[2]);
        more.data <- data.frame(
                sampleID=rep(sampleID, nrow(transcriptExp)),
                probeSetID=as.character(transcriptExp[,1]),
                geneSymbols=transcriptAnnot$geneSymbols,
                targetID=transcriptAnnot$targetID,
                exprValues=as.numeric(transcriptExp[,2]));

        out.file <- paste0("target_nbl_", sampleID, "_affy_exon_expr_full.csv")
        write.table(more.data, file=out.file, sep=",", quote=FALSE,
                                row.names=FALSE, col.names=TRUE);
        target_nbl_affy_exp <- rbind(target_nbl_affy_exp, more.data);
}

dim(target_nbl_affy_exp)    #           [1] 68114199        5

save(target_nbl_affy_exp, file="target_nbl_affy_exp_exon_full_all_samples.RData")
write.table(target_nbl_affy_exp, sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE,
                        na="", file="target_nbl_affy_exp_exon_full_level_3_all_samples.csv")


#       Merge all files as one if there are too many samples
#       =====================================================================
#   
#   #   check out header line for each file
#   #   ====================================
#   #
#   for i in *.csv; do head -n 1 $i; done
#
#   #   Write header line to a new file
#   #   =====================================
#   #
#   head -n 1 target_nbl_TARGET-30-PASLGS_affy_exon_expr_full.csv 
#   >target_nbl_expr_array_exon
#   
#   #   Concatenate all data lines from each file to the new file above
#   #   ===============================================================
#
#   [hzhang@helix Full]$ for i in *.csv; 
#   do 
#   echo $i; 
#   sed '1d' $i >>target_nbl_expr_array_exon;  
#   done
#
#   #   Data check
#   #   ===================================================
#   [hzhang@helix Full]$ wc -l *.csv
#    .
#    .
#    .
#    273552 target_nbl_TARGET-30-PASLGS_affy_exon_expr_full.csv
#    68114448 total
#
#   [hzhang@helix Full]$ wc -l target_nbl_expr_array_exon
#    68114200 target_nbl_expr_array_exon
# 
#   [hzhang@helix Full]$ ls -l *.csv | wc -l     
#   249
#
#   # 249*273552-248 = 68114200
#
#
#   #   Rename the new file above 
#   #   ==================================================
#
#   mv target_nbl_expr_array_exon /
#   target_nbl_affy_exp_exon_full_level_3_all_samples.csv
#
#   ==================================================================
#