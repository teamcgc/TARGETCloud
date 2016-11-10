#   Generate RNAseq data file from ccsk RNAseq data for google
#   cloud storage. There are three types of RNAseq data: gene,
#   exon, and isoform level counts.
#   ____________________________________________________________
#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


rm(list=ls(all=TRUE))
setwd("/data/CCRBioinfo/Henry/BigQuery.Data/target_ccsk/target_ccsk_rnaseq_expression")


#   Sample information
#   ==============================================================

sampleInfo <- read.table("TARGET_CCSK_mRNA-seq_20160826.sdrf.txt", 
    header=TRUE, sep="\t", quote="")
sampleInfo <- sampleInfo[-which(duplicated(sampleInfo[,1])), c(1,45)]
sampleInfo[,2] <- gsub("\\..{1,}", "", sampleInfo[,2])

dim(sampleInfo)     #   [1] 13  2
sampleInfo
#           Source.Name Derived.Array.Data.File
#   1  TARGET-51-PAEALX                 CCSK002
#   5  TARGET-51-PAJLIV                 CCSK003
#   9  TARGET-51-PAJLWU                 CCSK004
#   13 TARGET-51-PAJMFS                 CCSK005
#   17 TARGET-51-PAJMNM                 CCSK006
#   21 TARGET-51-PAJNCV                 CCSK007
#   25 TARGET-51-PAJPFB                 CCSK008
#   29 TARGET-51-PAKWMM                 CCSK009
#   33 TARGET-51-PALEIR                 CCSK010
#   37 TARGET-51-PALFEF                 CCSK011
#   41 TARGET-51-PALFYG                 CCSK012
#   45 TARGET-51-PALKEI                 CCSK013
#   49 TARGET-51-PALLXV                 CCSK014


#   Gene level FPKM counts
#   ==================================================

file.pattern <- ".gene.fpkm.txt"
gene.expr.files <- list.files(pattern=file.pattern);

for(aFile in 1:length(gene.expr.files))
{
    sampleName <- sub(file.pattern, "", gene.expr.files[aFile]);
    sampleID <- as.character(sampleInfo[which(sampleInfo[,2]==sampleName), 1]);
    print(paste(sampleName, sampleID));

    gene <- read.table(gene.expr.files[aFile], 
        header=TRUE, sep="\t", quote="");
    inData <- data.frame(sampleID=rep(sampleID, times=nrow(gene)),
        alternative_id=rep(sampleName, times=nrow(gene)),
        ensembl_gene_id=as.character(gene$gene_id),
        gene_symbols=as.character(gene$gene_short_name),
        fpkm=as.numeric(gene$FPKM));
        
    if(aFile==1) {
        target_ccsk_rnaseq_gene_exp <- inData;
    } else {
        target_ccsk_rnaseq_gene_exp <- rbind(target_ccsk_rnaseq_gene_exp, 
            inData);
    }
}

save(target_ccsk_rnaseq_gene_exp, file="target_ccsk_rnaseq_gene_exp.RData")
write.table(target_ccsk_rnaseq_gene_exp, sep=",", quote=FALSE, 
    row.names=FALSE, col.names=TRUE,
    file="target_ccsk_RNAseq_gene_expression_level_3_all_samples.csv")

rm(aFile); rm(file.pattern); rm(gene); rm(gene.expr.files); rm(inData); 
rm(sampleID); rm(sampleName); rm(target_ccsk_rnaseq_gene_exp)
ls()    #   [1] "sampleInfo"


#       Exon level RPKM counts. There is no header in data files. Columns in data files are:
#
#       chromosome (1:22, X, and Y)
#       start position
#       end position
#       gene symbols
#       Locus ID
#       AccessionNum
#       RPKM value
#
#       It seems that the most of 5th and 6th columns are same 
#       nrow(exon)                                                                              #       [1] 230470
#       sum(as.character(exon[,5]) == as.character(exon[,6]))   #       [1] 227436
#
#       Also, and 5th and 6th columns contains "," and should be replaced with ";"
#   ======================================================================

file.pattern <- ".exon.rpkm.txt"
exon.expr.files <- list.files(pattern=file.pattern);

for(aFile in 1:length(exon.expr.files))
{
        sampleName <- sub(file.pattern, "", exon.expr.files[aFile]);
        sampleID <- as.character(sampleInfo[which(sampleInfo[,2]==sampleName), 1]);
        print(paste(sampleName, sampleID));

        exon <- read.table(exon.expr.files[aFile], header=FALSE, sep="\t", quote="");
        exon[, 5] <- gsub(",$", "", exon[, 5]);
        exon[,5]  <- gsub(",", ";", exon[,5]);
        
        inData <- data.frame(sampleID=rep(sampleID, times=nrow(exon)),
                                                alternative_id=rep(sampleName, times=nrow(exon)),
                                                gene_symbols=as.character(exon[,4]),
                                                locus_id=as.character(exon[, 5]),
                                                rpkm=as.numeric(exon[,7])
                                                );
        
        if(aFile==1) {
                target_ccsk_rnaseq_exon_exp <- inData;
        } else {        
                target_ccsk_rnaseq_exon_exp <- rbind(target_ccsk_rnaseq_exon_exp, inData);
        }
}

save(target_ccsk_rnaseq_exon_exp, file="target_ccsk_rnaseq_exon_exp.RData")
write.table(target_ccsk_rnaseq_exon_exp, sep=",", quote=FALSE, 
                        row.names=FALSE, col.names=TRUE, na="",
                        file="target_ccsk_RNAseq_exon_expression_level_3_all_samples.csv")

rm(aFile); rm(exon); rm(exon.expr.files); rm(file.pattern); rm(inData); 
rm(sampleID); rm(sampleName); rm(target_ccsk_rnaseq_exon_exp)
ls()    #   [1] "sampleInfo"



#   isoform level RPKM counts
#   ==================================================

file.pattern <- ".isoform.fpkm.txt"
isoform.expr.files <- list.files(pattern=file.pattern);

for(aFile in 1:length(isoform.expr.files))
{
    sampleName <- sub(file.pattern, "", isoform.expr.files[aFile]);
    sampleID <- as.character(sampleInfo[which(sampleInfo[,2]==sampleName), 1]);
    print(paste(sampleName, sampleID));

    isoform <- read.table(isoform.expr.files[aFile], 
        header=TRUE, sep="\t", quote="");
    inData <- data.frame(sampleID=rep(sampleID, times=nrow(isoform)),
        alternative_id=rep(sampleName, times=nrow(isoform)),
        gene_symbols=as.character(isoform$gene_short_name),
        ensembl_gene_id=as.character(isoform$gene_id),
        fpkm=as.numeric(isoform$FPKM));
        
    if(aFile==1) {
        target_ccsk_rnaseq_isoform_exp <- inData;
    } else {
        target_ccsk_rnaseq_isoform_exp <- rbind(target_ccsk_rnaseq_isoform_exp, 
            inData);
    }
}

save(target_ccsk_rnaseq_isoform_exp, file="target_ccsk_rnaseq_isoform_exp.RData")
write.table(target_ccsk_rnaseq_isoform_exp, sep=",", quote=FALSE, 
    row.names=FALSE, col.names=TRUE,
    file="target_ccsk_RNAseq_isoform_expression_level_3_all_samples.csv")




