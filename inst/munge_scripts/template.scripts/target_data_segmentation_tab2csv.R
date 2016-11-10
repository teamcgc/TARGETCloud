#
#   Convert tab-delimited files to comma separated values format for copy
#   number change data. Columns retained are SampleName, Chromosome, Start,
#   End, NumProbes, and SementationMean.
#
#   If the data is already has the exact columns tab-delimited file will 
#   work for google storage and there is no need to convert to csv format
#   ________________________________________________________________________
#   <BigQuery CN Data><BigQuery CN Data><BigQuery CN Data><BigQuery CN Data>


setwd("/data/CCRBioinfo/Henry/BigQuery.Data/");

#   1.  AML samples. All samples are included in one table. Header raw   
#       contains a '#' character and  need to be skipped when read-in 
#   ===================================================================

data.file <- paste0("target_aml/target_aml_cn_level3/", 
    "TARGET_AML_CN_level3_filtered_Relapse.txt");
aml_cn <- read.table(data.file, header=FALSE, sep="\t", quote="",skip=1)

head(aml_cn)
out.data <- data.frame(SampleName=as.character(aml_cn[,6]), 
                        Chromosome=as.character(aml_cn[,1]), 
                        Start=as.numeric(aml_cn[,2]),
                        End=as.numeric(aml_cn[,3]),
                        NumProbes=as.numeric(aml_cn[,9]),
                        SegmentMean=as.numeric(aml_cn[,11]))
out.file <- sub("txt", "csv", data.file);
write.table(out.data, file=out.file, sep=",", quote=F,
            col.names=TRUE, row.names=FALSE)
rm(list=ls(all=TRUE))


data.file <- paste0("target_aml/target_aml_cn_level3/", 
    "TARGET_AML_CN_level3_filtered_Diagnostic.txt");
aml_cn <- read.table(data.file, header=FALSE, sep="\t", quote="",skip=1)

head(aml_cn)
out.data <- data.frame(SampleName=as.character(aml_cn[,6]), 
                        Chromosome=as.character(aml_cn[,1]), 
                        Start=as.numeric(aml_cn[,2]),
                        End=as.numeric(aml_cn[,3]),
                        NumProbes=as.numeric(aml_cn[,10]),
                        SegmentMean=as.numeric(aml_cn[,11]))
out.file <- sub("txt", "csv", data.file);
write.table(out.data, file=out.file, sep=",", quote=F,
            col.names=TRUE, row.names=FALSE)
rm(list=ls(all=TRUE))



#   2.  CCSK samples. Each sample has one table and need to be concatenated
#       as one table. Only six columns are in the tables
#   ========================================================================


setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
        "target_ccsk/target_ccsk_cn_array"));
in.files <- list.files(pattern="TARGET");

aFile <- 1;
ccsk_cn <- read.table(in.files[aFile], header=TRUE, sep="\t", quote="");
colnames(ccsk_cn) <- c("SampleName", "Chromosome", "Start",
                    "End", "NumProbes", "SementationMean");
ccsk_cn$SampleName <- gsub("\\.", "-", ccsk_cn$SampleName);

for(aFile in 2:length(in.files))
{
    more.data <- read.table(in.files[aFile], header=TRUE, sep="\t", quote="");
    colnames(more.data) <- c("SampleName", "Chromosome", "Start",
                    "End", "NumProbes", "SementationMean");
    more.data$SampleName <- gsub("\\.", "-", more.data$SampleName);

    ccsk_cn <- rbind(ccsk_cn, more.data)
}

write.table(ccsk_cn, file="target_ccsk_segmentation_level_3_all_samples.csv", 
        sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
rm(list=ls(all=TRUE))


#   3.  WT samples. Each sample has one table and need to be concatenated
#       as one table. Only six columns are in the tables
#   ========================================================================

setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
        "target_wt/target_wt_cn"));
in.files <- list.files(pattern="TARGET");

aFile <- 1;
wt_cn <- read.table(in.files[aFile], header=TRUE, sep="\t", quote="");
colnames(wt_cn) <- c("SampleName", "Chromosome", "Start",
                    "End", "NumProbes", "SementationMean");
wt_cn$SampleName <- gsub("\\.", "-", wt_cn$SampleName);

for(aFile in 2:length(in.files))
{
    more.data <- read.table(in.files[aFile], header=TRUE, sep="\t", quote="");
    colnames(more.data) <- c("SampleName", "Chromosome", "Start",
                    "End", "NumProbes", "SementationMean");
    more.data$SampleName <- gsub("\\.", "-", more.data$SampleName);

    wt_cn <- rbind(wt_cn, more.data)
}

write.table(wt_cn, file="target_wt_segmentation_level_3_all_samples.csv", 
        sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
rm(list=ls(all=TRUE))


#

#   4.  NBL samples. Data are saved in folders by sex and samples
#       (226 females and 273 males) and the file names have a pattern
#       of "purity_20__smapleName.bam.cov.cnvs.xls"
#   ==================================================================


setwd("/data/CCRBioinfo/Henry/BigQuery.Data/target_nbl/targeted_nbl_cn");

#   Female samples
#
in.files <- list.files(pattern="cnvs.xls", path="VisCap_Female_Somatic", 
                        recursive=TRUE, full.name=TRUE);

aFile <- 1;
nbl_cn <- read.table(in.files[aFile], header=TRUE, sep="\t", quote="")
colnames(nbl_cn)

 [1] "Sample"                "CNV_id"                "CNV"
 [4] "Genes"                 "Genome_start_interval" "Genome_end_interval"
 [7] "Start_interval"        "End_interval"          "Interval_count"
[10] "Min_log2ratio"         "Median_log2ratio"      "Max_log2ratio"
[13] "Loss_threshold"        "Gain_threshold"        "Batch_size"

nbl_cn <- data.frame(
        SampleName=gsub(".bam.cov", "", as.character(nbl_cn[,1])), 
        Chromosome=gsub(":.{1,}", "", as.character(nbl_cn[,5])), 
        Start=gsub(".{1,}:", "", as.character(nbl_cn[,5])),
        End=gsub(".{1,}-", "", as.character(nbl_cn[,6])),
        NumProbes=as.numeric(nbl_cn[,9]),
        SegmentMedian=as.numeric(nbl_cn[,11]))

nbl_cn$Start <- gsub("-.{1,}", "", nbl_cn$Start)
for(aFile in 2:length(in.files))
{
    more.data <- read.table(in.files[aFile], header=TRUE, sep="\t", quote="");
    more_cn <- data.frame(
        SampleName=gsub(".bam.cov", "", as.character(more.data[,1])), 
        Chromosome=gsub(":.{1,}", "", as.character(more.data[,5])), 
        Start=gsub(".{1,}:", "", as.character(more.data[,5])),
        End=gsub(".{1,}-", "", as.character(more.data[,6])),
        NumProbes=as.numeric(more.data[,9]),
        SegmentMedian=as.numeric(more.data[,11]))

    more_cn$Start <- gsub("-.{1,}", "", more_cn$Start)
    nbl_cn <- rbind(nbl_cn, more_cn);
}

save(nbl_cn, file="target_nbl_cn_female_somatic.RData");
write.table(nbl_cn, file="target_nbl_cn_female_somatic.csv", sep=",", 
        quote=FALSE, row.names=FALSE, col.names=TRUE);
rm(list=ls(all=TRUE))


#   Male samples
#
in.files <- list.files(pattern="cnvs.xls", path="VisCap_Male_Somatic", 
                        recursive=TRUE, full.name=TRUE);

aFile <- 1;
nbl_cn <- read.table(in.files[aFile], header=TRUE, sep="\t", quote="")
colnames(nbl_cn)

nbl_cn <- data.frame(
        SampleName=gsub(".bam.cov", "", as.character(nbl_cn[,1])), 
        Chromosome=gsub(":.{1,}", "", as.character(nbl_cn[,5])), 
        Start=gsub(".{1,}:", "", as.character(nbl_cn[,5])),
        End=gsub(".{1,}-", "", as.character(nbl_cn[,6])),
        NumProbes=as.numeric(nbl_cn[,9]),
        SegmentMedian=as.numeric(nbl_cn[,11]))

nbl_cn$Start <- gsub("-.{1,}", "", nbl_cn$Start)

for(aFile in 2:length(in.files))
{
    more.data <- read.table(in.files[aFile], header=TRUE, sep="\t", quote="");
    more_cn <- data.frame(
        SampleName=gsub(".bam.cov", "", as.character(more.data[,1])), 
        Chromosome=gsub(":.{1,}", "", as.character(more.data[,5])), 
        Start=gsub(".{1,}:", "", as.character(more.data[,5])),
        End=gsub(".{1,}-", "", as.character(more.data[,6])),
        NumProbes=as.numeric(more.data[,9]),
        SegmentMedian=as.numeric(more.data[,11]))

    more_cn$Start <- gsub("-.{1,}", "", more_cn$Start)
    nbl_cn <- rbind(nbl_cn, more_cn);
}

save(nbl_cn, file="target_nbl_cn_male_somatic.RData")
write.table(nbl_cn, file="target_nbl_cn_male_somatic.csv", sep=",", 
        quote=FALSE, row.names=FALSE, col.names=TRUE);
rm(list=ls(all=TRUE))

