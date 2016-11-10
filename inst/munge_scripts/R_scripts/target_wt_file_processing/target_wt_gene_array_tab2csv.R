#   4.  WT samples. Input data is saved as gct format
#
#   Covert WT gene expression array data to csv format. Input data
#   is saved as gct format with probeSetID and gene descriptions.
#   Platform is HG-U133_Plus_2
#   ________________________________________________________________
#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm(list=ls(all=TRUE))

annotation.file <- paste0("/data/CCRBioinfo/Henry/BigQuery.Data/", 
    "annotations/source.text.files/HG-U133_Plus_2.na36.annot.csv");

HG_U133_Plus_2.na36.hg19 <- read.csv(annotation.file, header=TRUE,
    sep=",", quote="\"", comment.char="#");

the.order <- order(as.character(HG_U133_Plus_2.na36.hg19$Probe.Set.ID));
HG_U133_Plus_2.na36.hg19 <- HG_U133_Plus_2.na36.hg19[the.order,];

colnames(HG_U133_Plus_2.na36.hg19);

#    [1] "Probe.Set.ID"                     "GeneChip.Array"
#    [3] "Species.Scientific.Name"          "Annotation.Date"
#    [5] "Sequence.Type"                    "Sequence.Source"
#    [7] "Transcript.ID.Array.Design."      "Target.Description"
#    [9] "Representative.Public.ID"         "Archival.UniGene.Cluster"
#   [11] "UniGene.ID"                       "Genome.Version"
#   [13] "Alignments"                       "Gene.Title"
#   [15] "Gene.Symbol"                      "Chromosomal.Location"
#   [17] "Unigene.Cluster.Type"             "Ensembl"
#   [19] "Entrez.Gene"                      "SwissProt"
#   [21] "EC"                               "OMIM"
#   [23] "RefSeq.Protein.ID"                "RefSeq.Transcript.ID"
#   [25] "FlyBase"                          "AGI"
#   [27] "WormBase"                         "MGI.Name"
#   [29] "RGD.Name"                         "SGD.accession.number"
#   [31] "Gene.Ontology.Biological.Process" "Gene.Ontology.Cellular.Component"
#   [33] "Gene.Ontology.Molecular.Function" "Pathway"
#   [35] "InterPro"                         "Trans.Membrane"
#   [37] "QTL"                              "Annotation.Description"
#   [39] "Annotation.Transcript.Cluster"    "Transcript.Assignments"
#   [41] "Annotation.Notes"

probeSetID <- as.character(HG_U133_Plus_2.na36.hg19$Probe.Set.ID);
targetID <- as.character(HG_U133_Plus_2.na36.hg19$Ensembl);
targetID <- gsub(" [/]{1,} ", ";", targetID);
for(aID in 1:length(targetID))
{
	gene_list <- unlist(strsplit(targetID[aID], ";"))
	ensembl_gene <- grep("ENSG", gene_list);
	
	if(length(targetID) > 0) {
		targetID[aID] <- paste(gene_list[ensembl_gene], collapse=";");
	} else { stop("Missing Ensembl ID.") }
}


#   Convert geneBySample table to contentsBySample table
#   =====================================================

setwd(paste0("/data/CCRBioinfo/Henry/BigQuery.Data/target_wt/", 
        "target_wt_gene_array"))

geneExp <- read.table("TARGET_WT_GEX_L3_NotCollapsed_20130912.gct", 
                    skip=2, header=TRUE, sep="\t", quote="")
geneExp <- geneExp[order(as.character(geneExp[,1])),];

totalRows <- nrow(geneExp)
totalSamples <- ncol(geneExp) - 2;

dim(HG_U133_Plus_2.na36.hg19)                   #   [1] 54675    41
dim(geneExp)                                    #   [1] 54675   130
sum(as.character(geneExp[,1]) == probeSetID)    #   [1] 54675


#   Reformat the sample IDs
#
sampleNames <- colnames(geneExp)[-c(1,2)]
sampleNames <- gsub("\\.", "-", sampleNames);
sampleID <- rep(sampleNames, each=totalRows)

probeSetID <- as.character(geneExp[,1]);
probeSetID <- rep(probeSetID, times=totalSamples)

geneSymbols <- as.character(geneExp[,2])
geneSymbols <- gsub("\"", "", geneSymbols)
geneSymbols <- gsub(".{1,}, ", "", geneSymbols)
geneSymbols <- rep(geneSymbols, times=totalSamples)

targetID <- rep(targetID, times=totalSamples)
exprValues <- log2(as.numeric(as.matrix(geneExp[,-c(1:2)])))

target_wt_affy_exp <- data.frame(sampleID, probeSetID, geneSymbols,
                                targetID, exprValues)

dim(target_wt_affy_exp) #   [1] 6998400       5
totalRows*totalSamples  #   [1] 6998400

save(target_wt_affy_exp, 
    file="target_wt_affy_exp_gene_level_3_all_samples.RData")
write.table(target_wt_affy_exp, sep=",", quote=FALSE, na="",
    file="target_wt_affy_exp_gene_level_3_all_samples.csv",
    col.names=TRUE, row.names=FALSE);
#


#   Data check
#   ==========================================================

length(sampleID)
length(probeSetID)
length(targetID)
length(geneSymbols)
length(exprValues)

unique(sampleID[1:length(sampleID)])

geneNames <- as.character(geneExp[,2])
geneNames <- gsub("\"", "", geneNames)
geneNames <- gsub(".{1,}, ", "", geneNames)

for(x in 1:totalSamples)
{
    first <- (x-1) * totalRows + 1;
    last  <- x*totalRows;
    columns <- x + 2;

    print(unique(sampleID[first:last]));
    same <- sum(probeSetID[first:last] == as.character(geneExp[,1]));
    message("    probe set ID same: ", same)

    same <- sum(geneSymbols[first:last] == geneNames);
    message("    geneSymbols same: ", same)

    same <- sum(exprValues[first:last] == log2(geneExp[,columns]))
    message("    Expr values same: ", same)
}

#   [1] "TARGET-50-PALLFB-01A-01R"
#       probe set ID same: 54675
#       geneSymbols same: 54675
#       Expr values same: 54675
#   ________________________________________________________________________
