#
#   array type defined is dismatch between cel files 
#   and TARGET_AML_GeneExpressionArray_20160812.sdrf.txt
#
#   sdrf file:  HuGene-1_0-st-v1.
#   cel files:  HuGene-1_1-st-v1
#   _________________________________________________________________
#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
                "target_aml/target_aml_gene_array"));


#   Annotation data 
#   ==============================================================
#
source(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/",
        "R_scripts/affymetrix_huex_probe2gene.R"));

annotation.file <- paste0("/data/CCRBioinfo/Henry/BigQuery.Data/",
        "annotations/source.text.files/", 
        "HuGene-1_1-st-v1.na36.hg19.probeset.csv");
HuGene.1.1.st.v1.na36.hg19 <- affymatrix.csv.probe2gene(
        annotation.file=annotation.file, ensembl="ENST");

save(HuGene.1.1.st.v1.na36.hg19, file="HuGene.1.1.st.v1.na36.hg19.RData");
write.table(HuGene.1.1.st.v1.na36.hg19, sep="\t", quote=FALSE,
        col.names=FALSE, row.names=TRUE, 
        file="HuGene.1.1.st.v1.na36.hg19.txt");

colnames(HuGene.1.1.st.v1.na36.hg19);

#   [1] "probeSetID"   "transcriptID" "geneSymbols"  "targetID"


#   gene array expression data
#   ==============================================================
#
geneExp <- read.table("TARGET_AML_GE_level3.txt", header=TRUE, 
                    sep="\t", quote="");
NL.BM.columns <- grep("NL.BM", colnames(geneExp));
geneExp <- geneExp[,-NL.BM.columns];
geneExp <- geneExp[order(as.character(geneExp[,1])), ];

   totalRows <- nrow(geneExp);          #   [1] 33297
totalSamples <- ncol(geneExp) - 2;      #   [1] 261


#   new data frame for google storage
#   ===============================================================
#    
redundents <- which(duplicated(HuGene.1.1.st.v1.na36.hg19$transcriptID));
annotation <- HuGene.1.1.st.v1.na36.hg19[-redundents,];
annotation <- annotation[order(as.character(annotation$transcriptID)), ]

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

#   An example of screen output
#   ============================
[1] "TARGET-20-PARMZF-09A-03R"
    probe set ID same: 33297
    geneSymbols same: 33297
    Expr values same: 33297


#   Save data in csv format
#   ===============================================================
#
target_aml_affy_exp <- data.frame(sampleID, probeSetID, geneSymbols,
                                targetID, exprValues);
dim(target_aml_affy_exp)    #   [1] 8690517       5

save(target_aml_affy_exp, file="target_aml_gene_exp_affy_level_3.RData");
write.table(target_aml_affy_exp, sep=",", quote=FALSE, 
        col.names=TRUE, row.names=FALSE,
        file="target_aml_gene_exp_affy_level_3_all_samples.csv");





























