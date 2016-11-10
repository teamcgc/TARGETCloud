#   File processing for methylation data in target project for google cloud
#   Storage and BigQuery
#
#   1.  Transform tab-delimited files to csv format and remove redundant 
#       gene symbols in each row such as NIPA2;NIPA2;NIPA2;NIPA2 -> NIPA2.
#       The output files are for local temporary storage only.
#
#   2.  Concatenate all samples into one table (sample by sample) for cloud
#       storage and BigQuery table generation
#
#   ==========================================================================
#

rm(list=ls(all=T))
data.dir <- "/data/CCRBioinfo/Henry/BigQuery.Data/";

#   sample.dir <- "target_os/target_os_methylation"
#         path <- "target_os_methyl_txt_format"
#   sample.dir <- "target_aml/target_aml_methyl_array";
#         path <- "target_aml_methylation_L3_txt_format"

#   sample.dir <- "target_nbl/target_nbl_methyl_array"
#         path <- "target_nbl_methyl_txt_format"

#   sample.dir <- "target_wt/target_wt_methyl_array"
#         path <- "target_wt_methyl_txt_format"


setwd(paste0(data.dir, sample.dir));
methyl.files <- list.files(pattern="TARGET", full.names=TRUE, path=path)

length(methyl.files)   
length(unique(methyl.files))


#   Manually check the file format:
#
#   There are two different platforms for methylation arrays:
#   27k and 450k which can be told with line numbers
#   ======================================================================
#
for(i in seq_along(methyl.files)) system(paste0("wc -l ", methyl.files[i]))


#   For 450k platform (315 files). There are only 5 columns in each file
#   Also need check how many columns in the file
#   ====================================================================

#   temp <- read.table(methyl.files[1], header=TRUE, sep="\t", 
#                           quote="", skip=1)
#   colnames(temp) <- c("ReporterID", "Signal", "GeneSymbols", 
#                        "Chromosome", "Position");                        
                        
temp <- read.table(methyl.files[1], header=TRUE, sep="\t", quote="")
colnames(temp) <- c("SampleID", "ReporterID", "Signal", "GeneSymbols", 
                        "Chromosome", "Position");

geneSymbols <- as.character(temp$GeneSymbols);
ReporterID <- as.character(temp$ReporterID);

for(aGene in seq_len(length(geneSymbols)))
{
    numOfGenes <- grep(";", geneSymbols[aGene]);
    if(length(numOfGenes)>0) 
    {
        genes <- strsplit(geneSymbols[aGene], split=";")
        genes <- unique(unlist(genes));

        if(length(genes) > 1) {
            genes <- paste0(genes, ";", collapse="");
            genes <- sub(";$", "", genes);
        }
        message(aGene, ": ", geneSymbols[aGene], " -> ", genes)           
        geneSymbols[aGene] <- genes;

    } else {message(aGene, ": ", geneSymbols[aGene]);}
}

#   file.name.prefix <- paste0("target_nbl_methyl_txt_format/", 
#           "usc.edu_NBL.HumanMethylation450.Level-3.")

file.name.prefix <- "target_wt_methyl_txt_format/"


for(aFile in seq_along(methyl.files))
{
    message(aFile, ": ", methyl.files[aFile]);
    #temp <- read.table(methyl.files[aFile], header=TRUE, sep="\t", 
    #                    quote="", skip=1)
    #colnames(temp) <- c("ReporterID", "Signal", "GeneSymbols", 
    #                    "Chromosome", "Position");

    temp <- read.table(methyl.files[aFile], header=TRUE, sep="\t", quote="")
    colnames(temp) <- c("SampleID", "ReporterID", "Signal", "GeneSymbols", 
                        "Chromosome", "Position");

    if(sum(ReporterID == temp$ReporterID) != length(ReporterID))
        stop("ReportID do not match!\n")
        
    temp$GeneSymbols <- geneSymbols;
    
    file.name <- sub("_L3.txt$", ".csv", methyl.files[aFile]);
    file.name <- sub(file.name.prefix, "", file.name);
    #   SampleID  <- sub(".csv", "", file.name);
    #   out <- data.frame(SampleID, temp)
    out <- temp;
    write.table(out, file=file.name, sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE);
}


#   R script for combining all files with one head lines
#   ==================================================================
#
#   sample.dir <- "target_os_methyl_csv_format"
#   sample.dir <- "target_ccsk_methyl_csv_format"
#   sample.dir <- "target_nbl_methyl_csv_format"
sample.dir <- "target_wt_methyl_csv_format"


input.files <- list.files(path=sample.dir, pattern="csv", full.names=TRUE)
length(input.files)


aFile <- 1;
message(aFile, ": ", input.files[aFile]);
bigData <- read.table(input.files[1], header=TRUE, sep=",", quote="")
dim(bigData)    

for(aFile in seq_len(length(input.files))[-1])
{
    message(aFile, ": ", input.files[aFile]);
    moreData <- read.table(input.files[aFile], header=TRUE, sep=",", quote="")
    bigData <- rbind(bigData, moreData);
}

#   data check
#   =====================================================
dim(bigData)                    
nrow(moreData)*length(input.files)   
head(bigData)


#   save(bigData, file="target_os_methylation_level_3_all_samples.RData")
#   write.table(bigData, file="target_os_methylation_level_3_all_samples.csv",
#        sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
#

file.name <- paste0(sub("_methyl_csv_format", "", sample.dir), 
                "_methylation_level_3_all_samples")
save(bigData, file=paste0(file.name, ".RData"))
write.table(bigData, file=paste0(file.name, ".csv"), sep=",",
        quote=FALSE, row.names=FALSE, col.names=TRUE)

#





#   For 27k platform (482 files). There are 6 columns in each file
#   ====================================================================

temp <- read.table(methyl.files[1], header=TRUE, sep="\t", quote="")
colnames(temp) <- c("SampleID", "ReporterID", "Signal", "GeneSymbols", 
                        "Chromosome", "Position");
                        
geneSymbols <- as.character(temp$GeneSymbols);
ReporterID <- as.character(temp$ReporterID) 

#   generate an unique set of gene symbols for each row
#
for(aGene in seq_len(length(geneSymbols)))
{
    numOfGenes <- grep(";", geneSymbols[aGene]);
    if(length(numOfGenes)>0) {
        genes <- strsplit(geneSymbols[aGene], split=";")
        genes <- unique(unlist(genes));

        if(length(genes) > 1) {
            genes <- paste0(genes, ";", collapse="");
            genes <- sub(";$", "", genes);
        }
        message(aGene, ": ", geneSymbols[aGene], " -> ", genes)           
        geneSymbols[aGene] <- genes;

    } else {message(aGene, ": ", geneSymbols[aGene]);}
}
    
#   Save each file in csv format comma
#
for(aFile in seq_along(methyl.files))
{
    message(aFile, ": ", methyl.files[aFile]);
    temp <- read.table(methyl.files[aFile], header=TRUE, sep="\t", quote="")
    colnames(temp) <- c("SampleID", "ReporterID", "Signal", "GeneSymbols", 
                        "Chromosome", "Position");

    if(sum(ReporterID == temp$ReporterID) != length(ReporterID))
        stop("ReportID do not match!\n")
        
    temp$GeneSymbols <- geneSymbols;    
    file.name <- sub("txt$", "csv", methyl.files[aFile]);
    file.name <- sub(".{1,}/", "", file.name);
    write.table(temp, file=file.name, sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE);
}

 
#   Merge all files as one if there are too many samples
#   =====================================================================
#
#   for i in *.csv; do head -n 1 $i; done
#
#   head -n 1 TARGET-21-PATKWH-15A-01D_450k_L3.csv >target_aml_methyl_array
#   for i in *.csv; 
#   do 
#       echo $i; 
#       sed '1d' $i >>target_wt_methyl_array;  
#   done
#   mv target_aml_methyl_array target_aml_methyl_array.csv

#   wc -l *27k_L3.csv                        13293078 total
#   wc -l *450k_L3.csv                      152936595 total
#   wc -l target_aml_methyl_array.csv       166228877 target_aml_methyl_array.csv

#   Total of 797 files 
#   $ ls -l *27k_L3.csv | wc -l     482
#   $ ls -l *450k_L3.csv | wc -l    315
#
#   450k file: 485513 lines in each file
#    27k file:  27579 lines in each file
#   152936595 + 13293078 - 79 - 796 = 166228877
#   ==================================================================
#