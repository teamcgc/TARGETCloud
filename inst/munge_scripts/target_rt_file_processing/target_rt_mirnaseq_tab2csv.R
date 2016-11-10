#
#   Convert miRNAseq data from txt format to csv format
#
#


setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
        "target_rt/target_rt_mirnaseq"))
rm(list=ls(all=TRUE))
 
#   1.  convert each file to csv format and add sample ID column
#   ============================================================
target_rt.files <- list.files(pattern="TARGET");
for(aFile in 1:length(target_rt.files))
{
    sampleID <- substr(target_rt.files[aFile], 1, 16);
    alternativeID <- substr(target_rt.files[aFile], 1, 19);
    rnaseq.data <- read.table(target_rt.files[aFile], 
                    header=TRUE, sep="\t", quote="");
    out.file <- sub("txt", "csv", target_rt.files[aFile]);

    print(aFile)
    print(target_rt.files[aFile])
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

length(list.files(pattern="mirna.quantification.csv"))      #   [1] 79
length(list.files(pattern="isoform.quantification.csv"))    #   [1] 79

#   expression at miRNAseq level
#   =========================================================

out.file <- "target_rt_miRNAseq_mirna_expression_level_3_all_samples.csv"
mirna.files <- list.files(pattern="mirna.quantification.csv")
write.header <- paste0("head -n 1 ", mirna.files[1], " >", out.file);
merge.file <- paste0("for i in *.mirna.quantification.csv; do echo $i;",
    " tail -n +2 $i >>", out.file, "; done")

system(write.header);
system(merge.file)

#   data check
#   ==========
system("ls -l *.mirna.quantification.csv | wc -l")  #   79
system("wc -l *.mirna.quantification.csv")          #   (1595 for each)
system("wc -l target_rt_miRNAseq_mirna_expression_level_3_all_samples.csv")
#   125927 target_rt_miRNAseq_mirna_expression_level_3_all_samples.csv
1595*79 - 78   #   [1] 125927

sampleID <- substr(target_rt.files[1], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   1594

sampleID <- substr(target_rt.files[length(target_rt.files)], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   1594

rm(out.file); rm(mirna.files); rm(write.header); rm(merge.file);



#   expression at isoform level (isoform files have different line numbers)
#   =======================================================================
out.file <- "target_rt_miRNAseq_isoform_expression_level_3_all_samples.csv"
isoform.files <- list.files(pattern="isoform.quantification.csv")
write.header <- paste0("head -n 1 ", isoform.files[1], " >", out.file);
merge.file <- paste0("for i in *.isoform.quantification.csv; do echo $i;",
    " tail -n +2 $i >>", out.file, "; done")

system(write.header);
system(merge.file)

system("ls -l *.isoform.quantification.csv | wc -l")   #   79
system("wc -l *.isoform.quantification.csv")           
#   9285 TARGET-00-NAAEMA-20.2A-01R.isoform.quantification.csv
#   .
#   .
#   6529 TARGET-52-PAWFWK-01A-01R.isoform.quantification.csv
#   555554 total
system("wc -l *.isoform.quantification.txt")           
#   9285 TARGET-00-NAAEMA-20.2A-01R.isoform.quantification.txt
#   .
#   .
#   6529 TARGET-52-PAWFWK-01A-01R.isoform.quantification.txt
#   555554 total
  
system("wc -l target_rt_miRNAseq_isoform_expression_level_3_all_samples.csv")
#   555476 target_rt_miRNAseq_isoform_expression_level_3_all_samples.csv
555554 - 79  #   [1] 555475

sampleID <- substr(target_rt.files[1], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   9284

sampleID <- substr(target_rt.files[length(target_rt.files)], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   6528






















