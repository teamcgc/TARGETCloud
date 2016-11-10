#
#   Convert miRNAseq data from txt format to csv format
#
#


setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
        "target_wt/target_wt_mirnaseq"))
rm(list=ls(all=TRUE))
 
#   1.  convert each file to csv format and add sample ID column
#   ============================================================
target_wt.files <- list.files(pattern="TARGET");
for(aFile in 1:length(target_wt.files))
{
    sampleID <- substr(target_wt.files[aFile], 1, 16);
    alternativeID <- substr(target_wt.files[aFile], 1, 20);
    rnaseq.data <- read.table(target_wt.files[aFile], 
                    header=TRUE, sep="\t", quote="");
    out.file <- sub("txt", "csv", target_wt.files[aFile]);

    print(aFile)
    print(target_wt.files[aFile])
    print(out.file)
    print("    ")
    
    out.data <- data.frame(
            sampleID=rep(sampleID, times=nrow(rnaseq.data)),
            alternativeID=rep(alternativeID, times=nrow(rnaseq.data)),
            rnaseq.data);
    write.table(out.data, file=out.file, sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE);
    
    rm(sampleID); rm(alternativeID); rm(rnaseq.data); 
    rm(out.file); rm(out.data);
}


#   2.  Merge all csv files of same count type as one table 
#   =======================================================

length(list.files(pattern="mirna.quantification.csv"))      #   [1] 138
length(list.files(pattern="isoform.quantification.csv"))    #   [1] 138

#   expression at miRNAseq level
#   =========================================================

out.file <- "target_wt_miRNAseq_mirna_expression_level_3_all_samples.csv"
mirna.files <- list.files(pattern="mirna.quantification.csv")
write.header <- paste0("head -n 1 ", mirna.files[1], " >", out.file);
merge.file <- paste0("for i in *.mirna.quantification.csv; do echo $i;",
    " tail -n +2 $i >>", out.file, "; done")

system(write.header);
system(merge.file)

#   data check
#   ==========
system("ls -l *.mirna.quantification.csv | wc -l")  #   138
system("wc -l *.mirna.quantification.csv")          #   (1871 for each)
system("wc -l target_wt_miRNAseq_mirna_expression_level_3_all_samples.csv")
#   258061 target_wt_miRNAseq_mirna_expression_level_3_all_samples.csv
1871*138 - 137   #   [1] 258061

sampleID <- substr(target_wt.files[1], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   1870

sampleID <- substr(target_wt.files[length(target_wt.files)], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   1870

rm(out.file); rm(mirna.files); rm(write.header); rm(merge.file);



#   expression at isoform level (isoform files have different line numbers)
#   =======================================================================
out.file <- "target_wt_miRNAseq_isoform_expression_level_3_all_samples.csv"
isoform.files <- list.files(pattern="isoform.quantification.csv")
write.header <- paste0("head -n 1 ", isoform.files[1], " >", out.file);
merge.file <- paste0("for i in *.isoform.quantification.csv; do echo $i;",
    " tail -n +2 $i >>", out.file, "; done")

system(write.header);
system(merge.file)

system("ls -l *.isoform.quantification.csv | wc -l")   #   138
system("wc -l *.isoform.quantification.csv")           
#   4444 TARGET-50-CAAAAA-01A-01R.isoform.quantification.csv
#   .
#   .
#   6339 TARGET-50-PALLFB-01A-01R.isoform.quantification.csv
#   615488 total
system("wc -l *.isoform.quantification.txt")           
#   4444 TARGET-50-CAAAAA-01A-01R.isoform.quantification.txt
#   .
#   .
#   6339 TARGET-50-PALLFB-01A-01R.isoform.quantification.txt
#   615488 total
  
system("wc -l target_wt_miRNAseq_isoform_expression_level_3_all_samples.csv")
#   615351 target_wt_miRNAseq_isoform_expression_level_3_all_samples.csv
615488 - 137  #   [1] 615351

sampleID <- substr(target_wt.files[1], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   4443

sampleID <- substr(target_wt.files[length(target_wt.files)], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   6338






















