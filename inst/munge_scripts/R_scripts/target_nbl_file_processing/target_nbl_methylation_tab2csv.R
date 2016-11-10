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
sample.dir <- "target_nbl/target_nbl_methyl_array"
path <- "target_nbl_methyl_txt_format"

setwd(paste0(data.dir, sample.dir));
methyl_files <- list.files(pattern="TARGET", full.names=TRUE, path=path)

length(methyl_files)                            #       [1] 235
length(unique(methyl_files))            #       [1] 235


#   Manually check the file format:
#
#   There are two different platforms for methylation arrays:
#   27k and 450k which can be told with line numbers
#   ======================================================================
#
for(i in seq_along(methyl.files)) system(paste0("wc -l ", methyl.files[i]))


#       For 450k platform (315 files). There are only 5 columns in each file
#       Also need check how many columns in the file
#   ====================================================================

aFile <- 1;
nbl_methyl_data <- read.table(methyl_files[aFile], skip=1, header=TRUE, sep="\t", quote="")
colnames(nbl_methyl_data) <- c("ReporterID", "Signal", "GeneSymbols",  "Chromosome", "Position");

#       remove redundant gene symbols in each row
#       ==============================================
geneSymbols <- as.character(nbl_methyl_data$GeneSymbols);
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
nbl_methyl_data$GeneSymbols <- geneSymbols;
nbl_methyl_data$Chromosome <- paste0("chr", nbl_methyl_data$Chromosome);

sampleID <- sub(".{1,}TARGET", "TARGET", methyl_files[aFile]);
sampleID <- sub(".txt", "", sampleID);
sampleID <- rep(sampleID, nrow(nbl_methyl_data));
nbl_methyl_data <- data.frame(sampleID,  nbl_methyl_data);

file_name <- sub(".{1,}TARGET", "TARGET", methyl_files[aFile]);
file_name <- sub("txt", "csv", file_name);

write.table(nbl_methyl_data, file=file_name, sep=",", quote=FALSE,
                        row.names=FALSE, col.names=TRUE, na="");

#       Convert all files to csv foramt and merge all object as one
#
ReporterID <- as.character(nbl_methyl_data$ReporterID);

for(aFile in 2:length(methyl_files))
{
        message(aFile, ": ", methyl_files[aFile]);
        methyl_data <- read.table(methyl_files[aFile], skip=1, header=TRUE, sep="\t", quote="")
        colnames(methyl_data) <- c("ReporterID", "Signal", "GeneSymbols", "Chromosome", "Position");

        if(sum(ReporterID == methyl_data$ReporterID) != length(ReporterID))
                stop("ReportID do not match!\n");

        sampleID <- sub(".{1,}TARGET", "TARGET", methyl_files[aFile]);
        sampleID <- sub(".txt", "", sampleID);
        sampleID <- rep(sampleID, nrow(methyl_data));

        methyl_data$GeneSymbols <- geneSymbols;
        methyl_data$Chromosome <- paste0("chr", methyl_data$Chromosome);
        methyl_data <- data.frame(sampleID, methyl_data);
 
        file_name <- sub(".{1,}TARGET", "TARGET", methyl_files[aFile]);
        file_name <- sub("txt", "csv", file_name);

        write.table(methyl_data, file=file_name, sep=",", quote=FALSE,
                        row.names=FALSE, col.names=TRUE, na="");
         nbl_methyl_data <- rbind(nbl_methyl_data, methyl_data);
}

dim(nbl_methyl_data);           #       
save(nbl_methyl_data, file="target_nbl_methyl_array_level_3_all_samples.RData");
write.table(nbl_methyl_data, file="target_nbl_methyl_array_level_3_all_samples.csv",
                sep=",", quote=FALSE, na="", row.names=FALSE, col.names=TRUE)





