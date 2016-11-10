#
#       Converte RNAseq data of RT samples to csv format
#       Sample ID will be extracted from file name and added to new 
#       table as one column.
#
#       All csv files will be concatenated together latter from shell
#
#   =================================================================


setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
        "target_rt/target_rt_rnaseq"))

length(list.files(pattern=".gene.quantification.txt"))      #   [1] 72
length(list.files(pattern=".exon.quantification.txt"))      #   [1] 72
length(list.files(pattern=".isoform.quantification.txt"))   #   [1] 78

options(scipen = 999)


#   1.  target_rt RNAseq exon level expression
#   =================================================================


#	filePattern <- "exon.quantification.txt";

convertTxt2CSV("exon.quantification.txt")
convertTxt2CSV("gene.quantification.txt")
convertTxt2CSV("isoform.quantification.txt")

convertTxt2CSV <- function(filePattern) {

	target_rt_files <- list.files(pattern=filePattern);
	table_name <- paste0("target_rt_RNAseq_", 
		sub("\\..{1,}", "", filePattern),
		"_expr_level_3_all_samples.csv")
	
	for(aFile in 1:length(target_rt_files))
	{
		message(aFile, ": ", target_rt_files[aFile]);
		out_file <- sub("txt", "csv", target_rt_files[aFile]);		
		sampleID <- substr(target_rt_files[aFile], 1, 16);
		alternativeID <- sub("R\\..{1,}", "R", target_rt_files[aFile]);
		
		rnaseq_data <- read.table(target_rt_files[aFile], 
                    header=TRUE, sep="\t", quote="");
		out_data <- data.frame(
            sampleID=rep(sampleID, times=nrow(rnaseq_data)),
            alternativeID=rep(alternativeID, times=nrow(rnaseq_data)),
            rnaseq_data);
		write.table(out_data, file=out_file, sep=",", quote=FALSE,
            na="", row.names=FALSE, col.names=TRUE);
    
		if(aFile == 1) {
			target_rt_RNAseq <- out_data;
		} else { 
			target_rt_RNAseq <- rbind(target_rt_RNAseq, out_data); 
		}

		rm(sampleID); rm(alternativeID); rm(rnaseq_data); 
		rm(out_file); rm(out_data);
	}

	save(target_rt_RNAseq, file=sub("csv", "RData", table_name));
	write.table(target_rt_RNAseq, file=table_name, sep=",", quote=FALSE,
		na="", row.names=FALSE, col.names=TRUE)
}

#	exon:		23922792    (332262*72 - 71)   
#	gene:		 4208400	( 58541*72 - 71)
#	isoform:	14350831	(183986*78 - 77)

 
 
#   2.  Shell script to concatenate all csv files as one
#   =====================================================================

#   RNAseq expression at gene level
#   =========================================================

out.file <- "target_rt_RNAseq_gene_expression_level_3_all_samples.csv"
gene.files <- list.files(pattern="gene.quantification.csv")
write.header <- paste0("head -n 1 ", gene.files[1], " >", out.file);
system(write.header);

merge.file <- paste0("for i in *.gene.quantification.csv; do echo $i;",
    " tail -n +2 $i >>", out.file, "; done")
system(merge.file)


system("ls -l *.gene.quantification.csv | wc -l")   #   72
system("wc -l *.gene.quantification.csv")           #   (58451 for each)
system("wc -l target_rt_RNAseq_gene_expression_level_3_all_samples.csv")
#   4208401 target_rt_RNAseq_gene_expression_level_3_all_samples.csv
58451*72 - 71   #   [1] 4208401


rm(out.file); rm(gene.files); rm(write.header); rm(merge.file);



#   RNAseq expression at exon level
#   =========================================================
out.file <- "target_rt_RNAseq_exon_expression_level_3_all_samples.csv"
exon.files <- list.files(pattern="exon.quantification.csv")
write.header <- paste0("head -n 1 ", exon.files[1], " >", out.file);
merge.file <- paste0("for i in *.exon.quantification.csv; do echo $i;",
    " tail -n +2 $i >>", out.file, "; done")

system(write.header);
system(merge.file)

system("ls -l *.exon.quantification.csv | wc -l")   #   72
system("wc -l *.exon.quantification.csv")           #   (332262 for each)
system("wc -l target_rt_RNAseq_exon_expression_level_3_all_samples.csv")
#   23922793 target_rt_RNAseq_exon_expression_level_3_all_samples.csv
332262*72 - 71  #   [1] 23922793

rm(out.file); rm(exon.files); rm(write.header); rm(merge.file);


#   RNAseq expression at isoform level. There are few normal samples
#   ====================================================================
out.file <- "target_rt_RNAseq_isoform_expression_level_3_all_samples.csv"
isoform.files <- list.files(pattern="isoform.quantification.csv")
write.header <- paste0("head -n 1 ", isoform.files[1], " >", out.file);
merge.file <- paste0("for i in *.isoform.quantification.csv; do echo $i;",
    " tail -n +2 $i >>", out.file, "; done")

system(write.header);
system(merge.file)

system("ls -l *.isoform.quantification.csv | wc -l")    #   78
system("wc -l *.isoform.quantification.csv")            #   (183986 for each)
system("wc -l target_rt_RNAseq_isoform_expression_level_3_all_samples.csv")
#   14350831 target_rt_RNAseq_isoform_expression_level_3_all_samples.csv
183986*78 - 77      #   [1] 14350831

rm(list=ls(all=TRUE))