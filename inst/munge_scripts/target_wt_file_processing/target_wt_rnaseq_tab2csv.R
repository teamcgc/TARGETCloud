#
#       Converte RNAseq data of WT samples to csv format
#       Sample ID will be extracted from file name and  
#       added to new table as one column.
#
#       All csv files will be concatenated together latter
#       from shell
#
#   =================================================================


setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
        "target_wt/target_wt_rnaseq"))

length(list.files(pattern=".gene.quantification.txt"))      #   [1] 136
length(list.files(pattern=".exon.quantification.txt"))      #   [1] 136
length(list.files(pattern=".isoform.quantification.txt"))   #   [1] 136

options(scipen = 999)

#   1.  target_wt RNAseq files: text to csv format
#   =================================================================

convertTxt2CSV("exon.quantification.txt")
convertTxt2CSV("gene.quantification.txt")
convertTxt2CSV("isoform.quantification.txt")


convertTxt2CSV <- function(filePattern) {

	target_wt_files <- list.files(pattern=filePattern);
	table_name <- paste0("target_wt_RNAseq_", 
		sub("\\..{1,}", "", filePattern),
		"_expr_level_3_all_samples.csv")
	
	for(aFile in 1:length(target_wt_files))
	{
		message(aFile, ": ", target_wt_files[aFile]);
		out_file <- sub("txt", "csv", target_wt_files[aFile]);		
		sampleID <- substr(target_wt_files[aFile], 1, 16);
		alternativeID <- sub("R\\..{1,}", "R", target_wt_files[aFile]);
		
		rnaseq_data <- read.table(target_wt_files[aFile], 
                    header=TRUE, sep="\t", quote="");
		out_data <- data.frame(
            sampleID=rep(sampleID, times=nrow(rnaseq_data)),
            alternativeID=rep(alternativeID, times=nrow(rnaseq_data)),
            rnaseq_data);
		write.table(out_data, file=out_file, sep=",", quote=FALSE,
            na="", row.names=FALSE, col.names=TRUE);
    
		if(aFile == 1) {
			target_wt_RNAseq <- out_data;
		} else { 
			target_wt_RNAseq <- rbind(target_wt_RNAseq, out_data); 
		}

		rm(sampleID); rm(alternativeID); rm(rnaseq_data); 
		rm(out_file); rm(out_data);
	}

	save(target_wt_RNAseq, file=sub("csv", "RData", table_name));
	write.table(target_wt_RNAseq, file=table_name, sep=",", quote=FALSE,
		na="", row.names=FALSE, col.names=TRUE)
}




#	target_wt.files <- list.files(pattern="TARGET");
#	for(aFile in 1:length(target_wt.files))
#	{
#	    sampleID <- substr(target_wt.files[aFile], 1, 16);
#	    alternativeID <- substr(target_wt.files[aFile], 1, 20);
#	    rnaseq.data <- read.table(target_wt.files[aFile], 
#	                    header=TRUE, sep="\t", quote="");
#	    out.file <- sub("txt", "csv", target_wt.files[aFile]);
#	
#	    print(aFile)
#	    print(target_wt.files[aFile])
#	    print(out.file)
#	    print("    ")
#	    
#	    out.data <- data.frame(
#	           sampleID=rep(sampleID, times=nrow(rnaseq.data)),
#	            alternativeID=rep(alternativeID, times=nrow(rnaseq.data)),
#				rnaseq.data);
#	    write.table(out.data, file=out.file, sep=",", quote=FALSE,
#	            row.names=FALSE, col.names=TRUE);
#	    
#	    rm(sampleID); rm(alternativeID); rm(rnaseq.data); 
#	    rm(out.file); rm(out.data);
#	}
#	


#   2.  Shell script to concatenate all csv files as one when necessary
#   =====================================================================

length(list.files(pattern="gene.quantification.csv"))       #   [1] 136
length(list.files(pattern="exon.quantification.csv"))       #   [1] 136
length(list.files(pattern="isoform.quantification.csv"))    #   [1] 136

#   RNAseq expression at gene level
#   =========================================================

out.file <- "target_wt_RNAseq_gene_expression_level_3_all_samples.csv"
gene.files <- list.files(pattern="gene.quantification.csv")
write.header <- paste0("head -n 1 ", gene.files[1], " >", out.file);
merge.file <- paste0("for i in *.gene.quantification.csv; do echo $i;",
    " tail -n +2 $i >>", out.file, "; done")

system(write.header);
system(merge.file)

#   data check
#   ==========
system("ls -l *.gene.quantification.csv | wc -l")   #   136
system("wc -l *.gene.quantification.csv")           #   (58451 for each)
system("wc -l target_wt_RNAseq_gene_expression_level_3_all_samples.csv")
#   7949201 target_wt_RNAseq_gene_expression_level_3_all_samples.csv
58451*136 - 71   #   [1] 7949201

sampleID <- substr(target_wt.files[1], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   58450

sampleID <- substr(target_wt.files[length(target_wt.files)], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   58450

rm(out.file); rm(gene.files); rm(write.header); rm(merge.file);



#   RNAseq expression at exon level
#   =========================================================
out.file <- "target_wt_RNAseq_exon_expression_level_3_all_samples.csv"
exon.files <- list.files(pattern="exon.quantification.csv")
write.header <- paste0("head -n 1 ", exon.files[1], " >", out.file);
merge.file <- paste0("for i in *.exon.quantification.csv; do echo $i;",
    " tail -n +2 $i >>", out.file, "; done")

system(write.header);
system(merge.file)

system("ls -l *.exon.quantification.csv | wc -l")   #   136
system("wc -l *.exon.quantification.csv")           #   (332262 for each)
system("wc -l target_wt_RNAseq_exon_expression_level_3_all_samples.csv")
#   45187497 target_wt_RNAseq_exon_expression_level_3_all_samples.csv
332262*136 - 135  #   [1] 45187497

sampleID <- substr(target_wt.files[1], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   332261

sampleID <- substr(target_wt.files[length(target_wt.files)], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   332261


rm(out.file); rm(exon.files); rm(write.header); rm(merge.file);


#   RNAseq expression at isoform level. There are few normal samples
#   ====================================================================
out.file <- "target_wt_RNAseq_isoform_expression_level_3_all_samples.csv"
isoform.files <- list.files(pattern="isoform.quantification.csv")
write.header <- paste0("head -n 1 ", isoform.files[1], " >", out.file);
merge.file <- paste0("for i in *.isoform.quantification.csv; do echo $i;",
    " tail -n +2 $i >>", out.file, "; done")

system(write.header);
system(merge.file)

system("ls -l *.isoform.quantification.csv | wc -l")    #   136
system("wc -l *.isoform.quantification.csv")            #   (183986 for each)
system("wc -l target_wt_RNAseq_isoform_expression_level_3_all_samples.csv")
#   25021961 target_wt_RNAseq_isoform_expression_level_3_all_samples.csv
183986*136 - 135      #   [[1] 25021961

sampleID <- substr(target_wt.files[1], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   183985

sampleID <- substr(target_wt.files[length(target_wt.files)], 1, 16);
check.sample <- paste0("grep ", sampleID, " ", out.file, "| wc -l")
system(check.sample)    #   183985



rm(list=ls(all=TRUE))