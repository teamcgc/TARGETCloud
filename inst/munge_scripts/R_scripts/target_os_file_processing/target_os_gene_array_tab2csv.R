#
#   Covert tab delimited text table that contains gene expression data 
#   from Affymetrix microarray to a csv file for Google BigQuery
#
#   Sample ID from RNAseq experiment
#
#   > as.character(unique(bigData$sampleID))
#    [1] "TARGET-40-0A4HLD-01A-01R" "TARGET-40-0A4HMC-01A-01R"
#    [3] "TARGET-40-0A4HX8-01A-01R" "TARGET-40-0A4HXS-01A-01R"
#    [5] "TARGET-40-0A4HY5-01A-01R" "TARGET-40-0A4I0Q-01A-01R"
#    [7] "TARGET-40-0A4I0W-01A-01R" "TARGET-40-0A4I3S-01A-01R"
#    [9] "TARGET-40-0A4I42-01A-01R" "TARGET-40-0A4I48-01A-01R"
#   [11] "TARGET-40-0A4I4E-01A-01R" "TARGET-40-0A4I4M-01A-01R"
#   [13] "TARGET-40-0A4I4O-01A-01R" "TARGET-40-0A4I5B-01A-01R"
#   [15] "TARGET-40-0A4I65-01A-01R" "TARGET-40-0A4I6O-01A-01R"
#   [17] "TARGET-40-0A4I9K-01A-01R" "TARGET-40-PAKFVX-01A-01R"
#   [19] "TARGET-40-PAKUZU-01A-01R" "TARGET-40-PAKXLD-01A-01R"
#   [21] "TARGET-40-PAKZZK-01A-01R" "TARGET-40-PALECC-01A-01R"
#   [23] "TARGET-40-PALFYN-01A-01R" "TARGET-40-PALHRL-01A-01R"
#   [25] "TARGET-40-PALKDP-01A-01R" "TARGET-40-PALKGN-01A-01R"
#   [27] "TARGET-40-PALWWX-01A-01R" "TARGET-40-PALZGU-01A-01R"
#   [29] "TARGET-40-PAMEKS-01A-01R" "TARGET-40-PAMHLF-01A-01R"
#   [31] "TARGET-40-PAMHYN-01A-01R" "TARGET-40-PAMJXS-01A-01R"
#   [33] "TARGET-40-PAMLKS-01A-01R" "TARGET-40-PAMRHD-01A-01R"
#   [35] "TARGET-40-PAMTCM-01A-01R" "TARGET-40-PAMYYJ-01A-01R"
#   [37] "TARGET-40-PANGPE-01A-01R" "TARGET-40-PANGRW-01A-01R"
#   [39] "TARGET-40-PANMIG-01A-01R" "TARGET-40-PANPUM-01A-01R"
#   [41] "TARGET-40-PANSEN-01A-01R" "TARGET-40-PANVJJ-01A-01R"
#   [43] "TARGET-40-PANXSC-01A-01R" "TARGET-40-PANZHX-01A-01R"
#   [45] "TARGET-40-PANZZJ-01A-01R" "TARGET-40-PAPFLB-01A-01R"
#   [47] "TARGET-40-PAPIJR-01A-01R" "TARGET-40-PAPKWD-01A-01R"
#   [49] "TARGET-40-PAPNVD-01A-01R" "TARGET-40-PAPWWC-01A-01R"
#   [51] "TARGET-40-PAPXGT-01A-01R" "TARGET-40-PARBGW-01A-01R"
#   [53] "TARGET-40-PARDAX-01A-01R" "TARGET-40-PARFTG-01A-01R"
#   [55] "TARGET-40-PARGTM-01A-01R" "TARGET-40-PARJXU-01A-01R"
#   [57] "TARGET-40-PARKAF-01A-01R" "TARGET-40-PASEBY-01A-01R"
#   [59] "TARGET-40-PASEFS-01A-01R" "TARGET-40-PASFCV-01A-01R"
#   [61] "TARGET-40-PASKZZ-01A-01R" "TARGET-40-PASNZV-01A-01R"
#   [63] "TARGET-40-PASRNE-01A-01R" "TARGET-40-PASSLM-01A-01R"
#   [65] "TARGET-40-PASUUH-01A-01R" "TARGET-40-PASYUK-01A-01R"
#   [67] "TARGET-40-PATAWV-01A-01R" "TARGET-40-PATEEM-01A-01R"
#   [69] "TARGET-40-PATJVI-01A-01R" "TARGET-40-PATKSS-01A-01R"
#   [71] "TARGET-40-PATMIF-01A-01R" "TARGET-40-PATMPU-01A-01R"
#   [73] "TARGET-40-PATMXR-01A-01R" "TARGET-40-PATPBS-01A-01R"
#   [75] "TARGET-40-PATUXZ-01A-01R" "TARGET-40-PAUBIT-01A-01R"
#   [77] "TARGET-40-PAUTWB-01A-01R" "TARGET-40-PAUTYB-01A-01R"
#   [79] "TARGET-40-PAUUML-01A-01R" "TARGET-40-PAUVUL-01A-01R"
#   [81] "TARGET-40-PAUXPZ-01A-01R" "TARGET-40-PAUYTT-01A-01R"
#   [83] "TARGET-40-PAVALD-01A-01R" "TARGET-40-PAVCLP-01A-01R"
#   [85] "TARGET-40-PAVDTY-01A-01R" "TARGET-40-PAVECB-01A-01R"


setwd("/data/CCRBioinfo/Henry/BigQuery.Data")

geneExp <- read.table("gene_core_rma_summary_annot.txt", 
    header=TRUE, sep="\t", quote="")
#

dim(geneExp)                        #   [1] 22011    91
totalProbes  <- nrow(geneExp);
totalSamples <- ncol(geneExp)-2;
platform     <- "HuEx.1_0.st.v2.01.1";
probeSet     <- as.character(geneExp[,1]);

#   Reformat the sample IDs
#
sampleNames <- colnames(geneExp)[-c(1,2)]
sampleNames <- gsub(".{1,}HuEx.1_0.st.v2.01.1_.", "TARGET-40-", sampleNames)
sampleNames <- gsub("..CEL", "", sampleNames)
sampleNames <- gsub("\\.", "-", sampleNames)
sampleNames <- gsub("01A01R", "-01A-01R", sampleNames)
sampleNames <- gsub("01A02R", "-01A-02R", sampleNames)
sampleNames <- gsub("_", "", sampleNames)
sampleNames <- gsub("-$", "", sampleNames)

shortName <- which(nchar(sampleNames)==16) 
shortName                       #   [1] 61 62 63 64 68 69 82 85 86
sampleNames[shortName]

#   [1] "TARGET-40-PALECC" "TARGET-40-PANGPE" 
#   [3] "TARGET-40-PAPFLB" "TARGET-40-PASNZV"
#   [5] "TARGET-40-PATMIF" "TARGET-40-PATMPU" 
#   [7] "TARGET-40-PATMXR" "TARGET-40-PATPBS"
#   [9] "TARGET-40-PATUXZ"

sampleNames[shortName] <- paste0(sampleNames[shortName], "-01A-01R")
sampleNames[shortName]

#   [1] "TARGET-40-PALECC-01A-01R" "TARGET-40-PANGPE-01A-01R"
#   [3] "TARGET-40-PAPFLB-01A-01R" "TARGET-40-PASNZV-01A-01R"
#   [5] "TARGET-40-PATMIF-01A-01R" "TARGET-40-PATMPU-01A-01R"
#   [7] "TARGET-40-PATMXR-01A-01R" "TARGET-40-PATPBS-01A-01R"
#   [9] "TARGET-40-PATUXZ-01A-01R"


annotation  <- as.character(geneExp[,2]);
geneNamess <- annotation;
targetID    <- annotation;
for(aRow in seq_along(annotation))
{
    if(annotation[aRow] == "") next;

    geneList <- unique(unlist(strsplit(annotation[aRow], " // ")));
    ENST <- grep("ENST", geneList);

    #   annotation are all target ID
    #   =================================
    if(length(ENST) == length(geneList))
    {
        geneNames[aRow] <- "";
        if(length(geneList) > 1) {
            targetID[aRow] <- paste(geneList, collapse=";");
        } else { targetID[aRow] <- geneList; }

    #   annotation are all geneNames (length(ENST) == 0)
    #   =================================================
    } else if (length(ENST) == 0){
        targetID[aRow] <- "";
        if(length(geneList) > 1) {
            geneNames[aRow] <- paste(geneList, collapse=";");
        } else {
            geneNames[aRow] <- geneList;
        }

    #   Mixed geneNames and targetIDs
    #   ================================
    } else {
        isoforms <- geneList[ENST];
        if(length(isoforms) > 1) {
            targetID[aRow] <- paste(isoforms, collapse=";");
        } else { targetID[aRow] <- isoforms; }
        
        genes <- geneList[-ENST];
        if(length(genes) > 1) {
            geneNames[aRow] <- paste(genes, collapse=";");
        } else { geneNames[aRow] <- genes;}
    }
}


#   Converte probeBySample table to sampleThenSample table
#   Columns will be:
#
#   sampleID
#   probeSet
#   geneSymbols
#   targetID
#   log2Value
#   ======================================================

totalProbes  <- nrow(geneExp);
totalSamples <- ncol(geneExp)-2;
platform     <- "HuEx.1_0.st.v2.01.1";
probeSet     <- as.character(geneExp[,1]);

sampleID    <- rep(sampleNames, each=totalProbes);
probeSetID  <- rep(probeSet, times=totalSamples);
geneSymbols <- rep(geneNames, times=totalSamples);
targetID    <- rep(targetID, times=totalSamples);
expValues   <- as.numeric(as.matrix(geneExp[,3:ncol(geneExp)]));

target_os_affy_exp <- data.frame(sampleID, probeSetID, geneSymbols,
                                targetID, expValues)
save(target_os_affy_exp, file="target_os_affy_exp.RData")
write.table(target_os_affy_exp, file="target_os_affy_exp.csv",
        sep=",", quote=FALSE, col.names=TRUE, row.names=FALSE)



#   data check
#   =================================
length(sampleID)    #   [1] 1958979
length(probeSetID)  #   [1] 1958979
length(targetID)    #   [1] 1958979
length(geneSymbols) #   [1] 1958979
length(expValues)   #   [1] 1958979

for(x in 1:totalSamples)
{
    first <- (x-1) * totalProbes + 1;
    last  <- x*totalProbes;
    columns <- x + 2;
    
    same <- sum(expValues[first:last] == geneExp[,columns])
    message(x, ": ", sampleNames[x], ": ", same)
}

#   1: TARGET-40-PATKSS-01A-01R: 22011
#   2: TARGET-40-PAUTWB-01A-01R: 22011
#   3: TARGET-40-PAUTYB-01A-01R: 22011
#   4: TARGET-40-PAUUML-01A-01R: 22011
#   5: TARGET-40-PAUVUL-01A-01R: 22011
#   6: TARGET-40-PAUXPZ-01A-01R: 22011
#   7: TARGET-40-PAUYTT-01A-01R: 22011
#   8: TARGET-40-PAVALD-01A-01R: 22011
#   9: TARGET-40-PAVCLP-01A-01R: 22011
#   10: TARGET-40-PAVDTY-01A-01R: 22011
#   11: TARGET-40-PAVECB-01A-01R: 22011
#   12: TARGET-40-PAUBIT-01A-01R: 22011
#   13: TARGET-40-PATXFN-01A-01R: 22011
#   14: TARGET-40-PAPNVD-01A-01R: 22011
#   15: TARGET-40-PAPXGT-01A-01R: 22011
#   16: TARGET-40-PASEBY-01A-01R: 22011
#   17: TARGET-40-PASFCV-01A-01R: 22011
#   18: TARGET-40-PASKZZ-01A-01R: 22011
#   19: TARGET-40-PASSLM-01A-01R: 22011
#   20: TARGET-40-PASEFS-01A-01R: 22011
#   21: TARGET-40-PANVJJ-01A-01R: 22011
#   22: TARGET-40-PAMHLF-01A-01R: 22011
#   23: TARGET-40-PASRNE-01A-01R: 22011
#   24: TARGET-40-PATJVI-01A-01R: 22011
#   25: TARGET-40-PAMRHD-01A-01R: 22011
#   26: TARGET-40-PASYUK-01A-01R: 22011
#   27: TARGET-40-PATAWV-01A-01R: 22011
#   28: TARGET-40-PANMIG-01A-01R: 22011
#   29: TARGET-40-PAKXLD-01A-01R: 22011
#   30: TARGET-40-PAMEKS-01A-01R: 22011
#   31: TARGET-40-PAMLKS-01A-01R: 22011
#   32: TARGET-40-PARDAX-01A-01R: 22011
#   33: TARGET-40-PARFTG-01A-01R: 22011
#   34: TARGET-40-PARGTM-01A-01R: 22011
#   35: TARGET-40-PARJXU-01A-01R: 22011
#   36: TARGET-40-PARKAF-01A-01R: 22011
#   37: TARGET-40-PANGRW-01A-02R: 22011
#   38: TARGET-40-PAKZZK-01A-01R: 22011
#   39: TARGET-40-PANZZJ-01A-01R: 22011
#   40: TARGET-40-PARBGW-01A-01R: 22011
#   41: TARGET-40-PANZHX-01A-01R: 22011
#   42: TARGET-40-PAPKWD-01A-01R: 22011
#   43: TARGET-40-PAPVYW-01A-01R: 22011
#   44: TARGET-40-PAPWWC-01A-01R: 22011
#   45: TARGET-40-PAMTCM-01A-01R: 22011
#   46: TARGET-40-PATEEM-01A-01R: 22011
#   47: TARGET-40-PAMHYN-01A-01R: 22011
#   48: TARGET-40-PAMJXS-01A-01R: 22011
#   49: TARGET-40-PAMYYJ-01A-01R: 22011
#   50: TARGET-40-PANPUM-01A-01R: 22011
#   51: TARGET-40-PANSEN-01A-01R: 22011
#   52: TARGET-40-PANXSC-01A-01R: 22011
#   53: TARGET-40-PALWWX-01A-01R: 22011
#   54: TARGET-40-PALZGU-01A-01R: 22011
#   55: TARGET-40-PALHRL-01A-01R: 22011
#   56: TARGET-40-PAKFVX-01A-01R: 22011
#   57: TARGET-40-PALFYN-01A-01R: 22011
#   58: TARGET-40-PALKDP-01A-01R: 22011
#   59: TARGET-40-PALKGN-01A-01R: 22011
#   60: TARGET-40-PAKUZU-01A-01R: 22011
#   61: TARGET-40-PALECC-01A-01R: 22011
#   62: TARGET-40-PANGPE-01A-01R: 22011
#   63: TARGET-40-PAPFLB-01A-01R: 22011
#   64: TARGET-40-PASNZV-01A-01R: 22011
#   65: TARGET-40-PAPIJR-01A-01R: 22011
#   66: TARGET-40-0A4HX8-01A-01R: 22011
#   67: TARGET-40-0A4HXS-01A-01R: 22011
#   68: TARGET-40-PATMIF-01A-01R: 22011
#   69: TARGET-40-PATMPU-01A-01R: 22011
#   70: TARGET-40-0A4I0Q-01A-01R: 22011
#   71: TARGET-40-0A4I0W-01A-01R: 22011
#   72: TARGET-40-0A4I3S-01A-01R: 22011
#   73: TARGET-40-0A4I42-01A-01R: 22011
#   74: TARGET-40-0A4I4O-01A-01R: 22011
#   75: TARGET-40-0A4I5B-01A-01R: 22011
#   76: TARGET-40-0A4I48-01A-01R: 22011
#   77: TARGET-40-0A4I4E-01A-01R: 22011
#   78: TARGET-40-0A4I65-01A-01R: 22011
#   79: TARGET-40-0A4I6O-01A-01R: 22011
#   80: TARGET-40-0A4I8U-01A-01R: 22011
#   81: TARGET-40-0A4HLD-01A-01R: 22011
#   82: TARGET-40-PATMXR-01A-01R: 22011
#   83: TARGET-40-0A4HY5-01A-01R: 22011
#   84: TARGET-40-0A4I4M-01A-01R: 22011
#   85: TARGET-40-PATPBS-01A-01R: 22011
#   86: TARGET-40-PATUXZ-01A-01R: 22011
#   87: TARGET-40-0A4I9K-01A-01R: 22011
#   88: TARGET-40-0A4HMC-01A-01R: 22011
#   89: TARGET-40-PASUUH-01A-01R: 22011
