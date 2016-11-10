#
#   Convert target ccsk gene array data to csv format
#   Original data is in gcl format and plat form is HG-U133_Plus_2.
#   Second column in expression data is gene description and gene
#   symbols, e.g., "discoidin domain receptor tyrosine kinase 1, DDR1"
#   
#   ____________________________________________________________________
#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_ccsk/", 
        "target_ccsk_gene_expression"))

geneExp <- read.table("TARGET_CCSK_GEX_L3_NotCollapsed_20130522.gct", 
                    skip=2, header=TRUE, sep="\t", quote="")

geneExp[,2] <- gsub(".{1,}, ", "", geneExp[,2]);
geneExp[,2] <- gsub("\"", "", geneExp[,2])
geneExp <- geneExp[order(as.character(geneExp[,1])), ]

head(geneExp)[,1:3]

#          Name Description TARGET.51.PAEALX.01A.01R
#   1 1007_s_at        DDR1               131.145688
#   2   1053_at        RFC2               302.388214
#   3    117_at       HSPA6                23.271172
#   4    121_at        PAX8               175.657423
#   5 1255_g_at      GUCA1A                 5.528474
#   6   1294_at        UBA7                90.329546

   totalRows <- nrow(geneExp)
totalSamples <- ncol(geneExp) - 2;

annotation.file <- paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
    "annotations/source.text.files/HG-U133_Plus_2.na36.annot.csv")
annotation <- read.csv(annotation.file, header=TRUE, sep=",", quote="\"",
                            comment.char="#")
annotation <- annotation[order(as.character(annotation[,1])), ];

dim(geneExp)                            #   [1] 54675    15
dim(annotation)                         #   [1] 54675    41
sum(geneExp[,1] == annotation[,1])      #   [1] 54675


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
#   ============================================================
#
sampleNames <- colnames(geneExp)[-c(1,2)]
sampleNames <- gsub("\\.", "-", sampleNames);
sampleID <- rep(sampleNames, each=totalRows)

probeSetID  <- rep(as.character(geneExp[,1]), times=totalSamples)
geneSymbols <- rep(as.character(geneExp[,2]), times=totalSamples)
   targetID <- rep(targetID, times=totalSamples)
exprValues  <- as.numeric(as.matrix(geneExp[,-c(1:2)]))


#   Data check
#   ================================
#
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

#   Sample of screen output
#   ==============================
#   [1] "TARGET-51-PALLXV-01A-01R"
#       probe set ID same: 54675
#       geneSymbols same: 54675
#       Expr values same: 54675


target_ccsk_affy_exp <- data.frame(sampleID, probeSetID, geneSymbols,
                            targetID, log2(exprValues))
                            
dim(target_ccsk_affy_exp)    #   [1] 710775       5

save(target_ccsk_affy_exp, file="target_ccsk_gene_exp_affy_level_3.RData")
write.table(target_ccsk_affy_exp, sep=",", quote=FALSE, 
        col.names=TRUE, row.names=FALSE,
        file="target_ccsk_gene_exp_affy_level_3_all_samples.csv")

rm(list=ls(all=TRUE))
