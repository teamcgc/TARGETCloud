#
#	Convert txt format file to csv format and merge all files as one
#	for targt_wt_cn files
#
#



setwd("/data/CCRBioinfo/Henry/BigQuery.Data/target_wt/target_wt_cn")
target_wt_cn_files <- list.files(pattern="TARGET");

for(aFile in 1:length(target_wt_cn_files))
{
	message(aFile, ": ", target_wt_cn_files[aFile]);
	cn_data <- read.table(target_wt_cn_files[aFile], header=TRUE,
						sep="\t", quote="");
	colnames(cn_data) <- c("sampleID", "Chromosome", "locStart", "locEnd", 
							"numMark", "segMean");
	cn_data$sampleID <- gsub("\\.", "-", cn_data$sampleID)

	if(aFile == 1) {
		target_wt_cn_data <- cn_data;
	} else {
		target_wt_cn_data <- rbind(target_wt_cn_data, cn_data);
	}
}
target_wt_cn_data$Chromosome <- paste0("chr", target_wt_cn_data$Chromosome)

options(scipen = 999)
save(target_wt_cn_data, file="target_wt_cn_level_3_all_samples.RData");
write.table(target_wt_cn_data, sep=",", quote=FALSE, na="", 
		file="target_wt_segmentation_level_3_all_samples.csv",
		row.names=FALSE, col.names=TRUE);



