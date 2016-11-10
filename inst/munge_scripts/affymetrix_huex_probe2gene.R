#   Extract gene symbols and Ensembl ID from affymetrix annotation file
#
#   Input file is in csv format downloaded from Affymetrix
#
#   annot.file <- paste0("/data/CCRBioinfo/Henry/BigQuery.Data/annotations/",
#                       "HuGene-1_1-st-v1.na36.hg19.probeset.csv")
#   hugene.1.1.st.v1 <- read.csv(annot.file, header=TRUE, sep=",", quote="\"",
#                        comment.char="#")
#   ___________________________________________________________________________
#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



affymatrix.csv.probe2gene <- function(annotation.file=NULL,ensembl=NULL)
{
        annotation <- read.csv(annotation.file, header=TRUE, sep=",", 
                                quote="\"", comment.char="#");
        gene.assignment <- as.character(annotation$gene_assignment);

        geneSymbols <- rep("", nrow(annotation));
        targetID    <- rep("", nrow(annotation));

        for(aRow in seq_along(gene.assignment))
        {
                if(gene.assignment[aRow] == "---") next;

                #       gene symbols are in the second field of gene_assigment
                #
                geneList <- gsub("/", "", gene.assignment[aRow])
                geneList <- unique(unlist(strsplit(geneList, "  ")));
                geneSymbols[aRow] <- geneList[2];

                #       transcript IDs start with "ENST" (argument of ensembl)
                #
                ENST <- grep(ensembl, geneList);
                if(length(ENST) == 0)  next;
                
                if(length(ENST)==1) {   targetID[aRow] <- geneList[ENST];
                } else {  targetID[aRow] <- paste(geneList[ENST], collapse=";");}

                 message(geneSymbols[aRow], "  ", targetID[aRow]);
                 }

                probeSetID <- as.character(annotation$probeset_id);
                transcriptID <- as.character(annotation$transcript_cluster_id);
                annot.out <- data.frame(probeSetID, transcriptID,
                                geneSymbols=geneSymbols, targetID=targetID);

        return(annot.out);
}





