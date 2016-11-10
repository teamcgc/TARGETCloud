#   Generate RNAseq data file from NBL RNAseq data for google
#   cloud storage. 
#
#   There are two data sets: BCCA and NCI_Khan.
#
#   NCI_Khan data has three types of RNAseq data: gene, exon, 
#   and isoform level counts.
#
#   BCCA data has three types of RNASeq data: gene, exon, and 
#   splijxn lvele counts. BCCA data will not be used as it has 
#   only 10 samples.
#   ____________________________________________________________
#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


rm(list=ls(all=TRUE))



#   1.  Sample information. There are replicates for some samples
#   ==============================================================

setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_nbl/", 
        "target_nbl_rnaseq"))

sampleInfo <- read.table("TARGET_NBL_mRNA-seq_20160826.sdrf.txt", 
    header=TRUE, sep="\t", quote="")
dim(sampleInfo)     #   [1] 732  55

sampleInfo <- sampleInfo[, c(1,54)]
sampleInfo[,2] <- gsub("\\..{1,}", "", sampleInfo[,2])
sampleInfo <- sampleInfo[-which(duplicated(sampleInfo[,2])),]
dim(sampleInfo)     #   [1] 181   2

head(sampleInfo)
#           Source.Name Derived.Array.Data.File.1
#   1  TARGET-30-PAIFXV                    NB2240
#   5  TARGET-30-PAIPGU                    NB2059
#   9  TARGET-30-PAISNS                    NB2325
#   13 TARGET-30-PAITCI                    NB2026
#   17 TARGET-30-PAITEG                    NB2024
#   21 TARGET-30-PAIVHE                    NB2023



#   2.  Gene level FPKM counts. Files are in two directories and
#       one file defines counts as FPKM and one as RPKM. Here 
#       they are treated as same.
#   ===============================================================

#   NCI_Khan samples
#
setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_nbl/", 
        "target_nbl_rnaseq/NCI_Khan"))
system("wc -l *gene.fpkm.txt")

#      .
#      .
#      56146 NB2331.gene.fpkm.txt
#      9039506 total

file.pattern <- "gene.fpkm.txt"
gene.expr.files <- list.files(pattern=file.pattern);

for(aFile in 1:length(gene.expr.files))
{
    sampleName <- sub("\\..{1,}", "", gene.expr.files[aFile]);
    sampleID <- as.character(sampleInfo[which(sampleInfo[,2]==sampleName), 1]);
    print(paste(aFile, sampleName, sampleID));

    gene <- read.table(gene.expr.files[aFile], 
        header=TRUE, sep="\t", quote="");
    inData <- data.frame(sampleID=rep(sampleID, times=nrow(gene)),
        alternative_id=rep(sampleName, times=nrow(gene)),
        ensembl_gene_id=as.character(gene$gene_id),
        gene_symbols=as.character(gene$gene_short_name),
        fpkm=as.numeric(gene$FPKM));

    #write.table(inData, file=paste0(sampleName, ".rnaseq.gene.expr.csv"), 
    #    sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE);
    
    if(aFile == 1) {
        target_nbl_rnaseq_gene <- inData;
    } else {
        target_nbl_rnaseq_gene <- rbind(target_nbl_rnaseq_gene, inData);
    }
}
dim(target_nbl_rnaseq_gene)             #   [1] 9039345       5
system("wc -l NB2300.gene.fpkm.txt")    #   56146 NB2300.gene.fpkm.txt
56146*161-161                           #   [1] 9039345

save(target_nbl_rnaseq_gene, file="target_nbl_NCI_Khan_rnaseq_gene.RData")
write.table(target_nbl_rnaseq_gene, sep=",", quote=FALSE, 
    row.names=FALSE, col.names=TRUE, 
    file="target_nbl_NCI_Khan_rnaseq_gene_level_3_all_samples.csv")


#   BCCA samples. Columns are different with NCI_Khan samples
#
setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_nbl/", 
        "target_nbl_rnaseq/BCCA"))
system("wc -l *gene.quantification.txt")

#      .
#      .
#      .
#      51738 HS1803.gene.quantification.txt
#      517380 total

file.pattern <- "gene.quantification.txt"
gene.expr.files <- list.files(pattern=file.pattern);

for(aFile in 1:length(gene.expr.files))
{
    sampleName <- sub("\\..{1,}", "", gene.expr.files[aFile]);
    sampleID <- as.character(sampleInfo[which(sampleInfo[,2]==sampleName), 1]);
    print(paste(sampleName, sampleID));

    gene <- read.table(gene.expr.files[aFile], 
        header=TRUE, sep="\t", quote="");

    inData <- data.frame(sampleID=rep(sampleID, times=nrow(gene)),
        alternative_id=rep(sampleName, times=nrow(gene)),
        ensembl_gene_id=gsub(".{1,}[|]", "", gene$gene),
        gene_symbols=gsub("[|].{1,}", "", gene$gene),
        fpkm=as.numeric(gene$RPKM));

    target_nbl_rnaseq_gene <- rbind(target_nbl_rnaseq_gene, inData);
}

dim(target_nbl_rnaseq_gene)     #   [1] 9556715       5
9039506 - 161 + 517380 - 10     #   [1] 9556715

save(target_nbl_rnaseq_gene, file="target_nbl_rnaseq_gene.RData")
write.table(target_nbl_rnaseq_gene, sep=",", quote=FALSE, row.names=FALSE,
    col.names=TRUE, file="target_nbl_rnaseq_gene_level_3_all_samples.csv")



#   3.  isoform level RPKM counts. The data is contained in NCI_Khan
#       samples only
#   ====================================================================
setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_nbl/", 
        "target_nbl_rnaseq/NCI_Khan"))
        
file.pattern <- ".isoform.fpkm.txt"
isoform.expr.files <- list.files(pattern=file.pattern);

for(aFile in 1:length(isoform.expr.files))
{
    sampleName <- sub(file.pattern, "", isoform.expr.files[aFile]);
    sampleID <- as.character(sampleInfo[which(sampleInfo[,2]==sampleName), 1]);
    print(paste(aFile, sampleName, sampleID));

    isoform <- read.table(isoform.expr.files[aFile], 
        header=TRUE, sep="\t", quote="");
    inData <- data.frame(sampleID=rep(sampleID, times=nrow(isoform)),
        alternative_id=rep(sampleName, times=nrow(isoform)),
        target_id=as.character(isoform$tracking_id),
        length=as.numeric(isoform$length),
        geneSymbol=as.character(isoform$gene_short_name),
        fpkm=as.numeric(isoform$FPKM));
        
    if(aFile==1) {
        target_nbl_rnaseq_isoform_exp <- inData;
    } else {
        target_nbl_rnaseq_isoform_exp <- rbind(target_nbl_rnaseq_isoform_exp, 
            inData);
    }
}

dim(target_nbl_rnaseq_isoform_exp)      #   [1] 30317427        6
length(isoform.expr.files)              #   [1] 161
system("wc -l NB2331.isoform.fpkm.txt") #   188308 NB2331.isoform.fpkm.txt
188308*161-161                          #   [1] 30317427

#       Replace scientific notation with real numbers
#
index <- which(target_nbl_rnaseq_isoform_exp$fpkm == 10000000)
target_nbl_rnaseq_isoform_exp[1753443,]
#                       sampleID alternative_id       target_id length geneSymbol  fpkm
#       1753443 TARGET-30-PAMEZH         NB2013 ENST00000461719     46      IGHJ4 1e+07

index <- which(target_nbl_rnaseq_isoform_exp$fpkm < 0.001)
length(index)           #       [1] 14764935
target_nbl_rnaseq_isoform_exp$fpkm[index] <- 0.001
target_nbl_rnaseq_isoform_exp[index[1:5],]
#                  sampleID alternative_id       target_id length     geneSymbol  fpkm
#       1  TARGET-30-PALBFW         NB2001 ENST00000492842    940        OR4G11P 0.001
#       3  TARGET-30-PALBFW         NB2001 ENST00000461467    590        FAM138A 0.001
#       4  TARGET-30-PALBFW         NB2001 ENST00000335137    918          OR4F5 0.001
#       12 TARGET-30-PALBFW         NB2001 ENST00000496488    457 RP11-34P13.9.1 0.001
#       13 TARGET-30-PALBFW         NB2001 ENST00000410691    104             U6 0.001

save(target_nbl_rnaseq_isoform_exp, 
    file="target_nbl_NCI_Khan_rnaseq_isoform_exp.RData");
write.table(target_nbl_rnaseq_isoform_exp, sep=",", quote=FALSE, 
    row.names=FALSE, col.names=TRUE, na="",
    file="target_nbl_NCI_Khan_RNAseq_isoform_expression_level_3_all_samples.csv")





#   3.  Exon level RPKM counts. 
#
#   There is no header in NCI_Khan data files. Columns in data files
#   are ensembl gene ID with exon number and raw counts. Last five
#   lines in each file are summary data
#
#   BCCA data has genomic interval info, raw counts, median length and 
#   RPKM values
#   ========================================================


setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_nbl/", 
        "target_nbl_rnaseq/NCI_Khan"))
        
file.pattern <- ".exon.count.txt"
exon.expr.files <- list.files(pattern=file.pattern);

aFile <- 1;
sampleName <- sub(file.pattern, "", exon.expr.files[aFile]);
sampleID <- as.character(sampleInfo[which(sampleInfo[,2]==sampleName), 1]);
print(paste(aFile, sampleName, sampleID));

exon <- read.table(exon.expr.files[aFile], 
        header=FALSE, sep="\t", quote="");
exon <- exon[1:(nrow(exon)-5),];
exonID <- as.character(exon[,1]);
ensembl_gene_id <- gsub(":.{1,}", "", exon[,1]);
exon_num <- as.numeric(gsub(".{1,}:", "", exon[,1]));

target_nbl_rnaseq_exon_exp <- data.frame(
        sampleID=rep(sampleID, times=nrow(exon)),
        alternative_id=rep(sampleName, times=nrow(exon)),
        ensembl_gene_id=ensembl_gene_id,
        exon_num=exon_num,
        exon_count=as.numeric(exon[,2]));

for(aFile in 2:length(exon.expr.files))
{
    sampleName <- sub(file.pattern, "", exon.expr.files[aFile]);
    sampleID <- as.character(sampleInfo[which(sampleInfo[,2]==sampleName), 1]);
    print(paste(aFile, sampleName, sampleID));

    exon <- read.table(exon.expr.files[aFile], 
        header=FALSE, sep="\t", quote="");
    exon <- exon[1:(nrow(exon)-5),];
    
    if(sum(exonID==as.character(exon[,1])) != length(exonID))
        stop("Row headers do not match!")
    
    inData <- data.frame(sampleID=rep(sampleID, times=nrow(exon)),
        alternative_id=rep(sampleName, times=nrow(exon)),
        ensembl_gene_id=ensembl_gene_id,
        exon_num=exon_num,
        exon_count=as.numeric(exon[,2]));

    target_nbl_rnaseq_exon_exp <- rbind(target_nbl_rnaseq_exon_exp, inData);
}

save(target_nbl_rnaseq_exon_exp, file="target_nbl_NCI_Khan_rnaseq_exon_exp.RData")
write.table(target_nbl_rnaseq_exon_exp, sep=",", quote=FALSE, 
    row.names=FALSE, col.names=TRUE,
    file="target_nbl_NCI_Khan_RNAseq_exon_expression_level_3_all_samples.csv")

