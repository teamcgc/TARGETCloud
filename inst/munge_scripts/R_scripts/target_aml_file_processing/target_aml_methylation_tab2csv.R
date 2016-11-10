#       Convert methylation data to csv format:
#       There are two different platforms for methylation array data: 27k and 450k
#       _________________________________________________________________________
#       <aml methylation array data> <aml methylation array data> <aml methylation array data>


setwd("/data/CCRBioinfo/Henry/BigQuery.Data/target_aml/target_aml_methyl_array/");


#       1.      27k array. Data is clean and need to modify the column headers and chromosome
#               names (add prefix of "chr")
#       ====================================================================
#
methyl_array_27K_files <- list.files(path="target_aml_methylation_L3_txt_format", 
                pattern="27k_L3.txt", full.names=TRUE);
length(methyl_array_27K_files)          #       [1] 482

a_file <- 1;
methyl_data <- read.table(methyl_array_27K_files[a_file],  header=TRUE, sep="\t", quote="");
colnames(methyl_data) <- gsub("\\.", "", colnames(methyl_data));
methyl_data$Chromosome <- paste0("chr", methyl_data$Chromosome);
methyl_data$Position <- as.character(methyl_data$Position)
head(methyl_data);

for(a_file in 2:length(methyl_array_27K_files))
{
        print(paste0(a_file, ": ", methyl_array_27K_files[a_file]));
        more_data <- read.table(methyl_array_27K_files[a_file],  header=TRUE, sep="\t", quote="");

        colnames( more_data) <- gsub("\\.", "", colnames( more_data));
         more_data$Chromosome <- paste0("chr",  more_data$Chromosome);
         more_data$Position <- as.character( more_data$Position)
 
        methyl_data <- rbind(methyl_data, more_data)
}

dim(methyl_data)                #       [1] 13292596        6
head(methyl_data)

#                         SampleID  ProbeName AVG_Beta GeneSymbols Chromosome Position
#       1 TARGET-20-PADZKD-09A-01D cg15749858  0.13781         7A5       chr7 20223521
#       2 TARGET-20-PADZKD-09A-01D cg22568540  0.95763        A1BG      chr19 63556658
#       3 TARGET-20-PADZKD-09A-01D cg03586879  0.05394       A2BP1      chr16  6008832
#       4 TARGET-20-PADZKD-09A-01D cg19378133  0.06182       A2BP1      chr16  6008836
#       5 TARGET-20-PADZKD-09A-01D cg12058490  0.71548         A2M      chr12  9159646
#       6 TARGET-20-PADZKD-09A-01D cg03490200  0.89509       A2ML1      chr12  8866463

write.table(methyl_data, file="target_aml_methyl_array_27k_level_3_all_samples.csv",
                        sep=",", na="",  quote=FALSE, row.names=FALSE, col.names=TRUE);
rm(list=ls(all=TRUE));



#       2.      450k array. Data has no samplID column and chromosome names have prefix already.
#               Gene symbols have lot of redundant separated by ";"
#       =========================================================================
#
methyl_array_450K_files <- list.files(path="target_aml_methylation_L3_txt_format", 
                pattern="450k_L3.txt", full.names=TRUE);
length(methyl_array_450K_files)          #       [1]    315

a_file <- 1;
sampleName <- sub(".{1,}/", "",  methyl_array_450K_files[a_file]);
sampleName <- sub("_450k_L3.txt", "", sampleName);

methyl_data <- read.table(methyl_array_450K_files[a_file],  header=FALSE,  skip=2,  sep="\t", quote="");
colnames(methyl_data) <- c("ProbeName", "AVG_Beta", "GeneSymbols", "Chromosome", "Position");
methyl_data$Position <- as.character(methyl_data$Position);
dim(methyl_data)        #       485512

methyl_data <- methyl_data[order(methyl_data$ProbeName), ];
probeName <- as.character(methyl_data$ProbeName);
geneSymbols <- as.character(methyl_data$GeneSymbols);

for(a_symbol in 1:length(geneSymbols))
{
        if(length(grep(";", geneSymbols[a_symbol])) > 0)
       {
                gene_names <- unique(strsplit(geneSymbols[a_symbol], ";")[[1]]);
                if(length(gene_names) >1)  gene_names <- paste(gene_names, collapse=";")
           
                geneSymbols[a_symbol] <- gene_names;
                print(geneSymbols[a_symbol]);
       }
}
methyl_data$GeneSymbols <- geneSymbols;

methyl_data <- data.frame(sampleID=rep(sampleName, nrow(methyl_data)), methyl_data);
head(methyl_data);

for(a_file in 2:length(methyl_array_450K_files))
{
        print(paste0(a_file, ": ", methyl_array_450K_files[a_file]));
        sampleName <- sub(".{1,}/", "",  methyl_array_450K_files[a_file]);
        sampleName <- sub("_450k_L3.txt", "", sampleName);

        more_data <- read.table(methyl_array_450K_files[a_file],  header=FALSE,  skip=2, sep="\t", quote="");
        colnames(more_data) <- c("ProbeName", "AVG_Beta", "GeneSymbols", "Chromosome", "Position");
        more_data <- more_data[order(more_data$ProbeName), ];
         if(sum(as.character(more_data$ProbeName) == probeName)  != length(probeName))
                stop("Probe name rows do not match.")

        more_data$Position <- as.character(more_data$Position);
        more_data$GeneSymbols <- geneSymbols;
        more_data <- data.frame(sampleID=rep(sampleName, nrow(more_data)), more_data);
 
        methyl_data <- rbind(methyl_data, more_data)
}

dim(methyl_data)                #       [1] 
head(methyl_data)

save(methyl_data, file="target_aml_methyl_array_450k_level_3_all_samples.RData")
write.table(methyl_data, file="target_aml_methyl_array_450k_level_3_all_samples.csv",
                        sep=",", na="", quote=FALSE, row.names=FALSE, col.names=TRUE);



#       The data is processed on three steps and saved separately
#
#       dim(part_one)           #       [1] 51949784        6
#       dim(part_two)           #       [1] 51949784        6
#       dim(part_three)         #       [1] 49036712        6
#       samples <- c(unique(part_one$sampleID), unique(part_two$sampleID), unique(part_three$sampleID))
#       length(samples)         #       [1] 315
#       aml_methy_array_450k <- rbind(part_one, part_two, part_three)
#       dim(aml_methy_array_450k)       #       [1] 152936280         6
#
#       save(aml_methy_array_450k, file="aml_methyl_array_450k_all_samples.RData")
#       write.table(aml_methy_array_450k, file="target_aml_methyl_array_450k_level_3_all_samples.csv",
#                               sep=",", na="", quote=FALSE, row.names=FALSE, col.names=TRUE);
#
#       There are probes with position value in scientific notation and should be replaced with digitals only
#       https://www.radc.rush.edu/samples/methylation_normalization/data/cpg_info/final_cpg_list.txt
#       cg21011121      2       165000000
#       cg04665351      12      3000000
#
#       sed s/3e+06/3000000/  <target_aml_methylation_array_450k_all_samples.csv >temp.csv
#       mv target_aml_methylation_array_450k_all_samples.csv target_aml_methylation_array_450k_temp
#       sed s/1.65e+08/165000000/  <temp.csv >target_aml_methylation_array_450k_all_samples.csv
#
