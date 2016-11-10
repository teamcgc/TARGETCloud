#       Convert methylation data to csv format:
#       There are two different platforms for methylation array data: 27k and 450k
#       _________________________________________________________________________
#       <aml methylation array data> <aml methylation array data> <aml methylation array data>


setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_wt/",
			"target_wt_methyl_array"));
options(scipen = 999)



#	450k array. Data has samplID column and chromosome names have no 
#	prefix. Gene symbols are fine.
#       =========================================================================
#
wt_methyl_files <- list.files(path="target_wt_methyl_txt_format", 
                pattern="_L3.txt", full.names=TRUE);
length(wt_methyl_files)          #       [1]    131

for(a_file in 1:length(wt_methyl_files))
{
	message(a_file, ": ", wt_methyl_files[a_file]);
	methyl_data <- read.table(wt_methyl_files[a_file],  
					header=TRUE,  sep="\t", quote="");
	colnames(methyl_data) <- c("SampelID", "ProbeName", 
		"AVG_Beta", "GeneSymbols", "Chromosome", "Position");
	methyl_data$Chromosome <- paste0("chr", methyl_data$Chromosome);
		
	if(a_file == 1) {
		target_wt_methyl_data <- methyl_data;
	} else { 
		target_wt_methyl_data <- rbind(target_wt_methyl_data, methyl_data);
	}
}

dim(target_wt_methyl_data)	#       [1] 

save(target_wt_methyl_data, 
	file="target_wt_methyl_array_level_3_all_samples.RData")
write.table(target_wt_methyl_data, sep=",", na="", quote=FALSE, 
	file="target_wt_methyl_array_level_3_all_samples.csv",
	row.names=FALSE, col.names=TRUE);



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
