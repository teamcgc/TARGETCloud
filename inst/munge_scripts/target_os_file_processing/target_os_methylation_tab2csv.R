#
#   Transform tab-delimited files to csv format and remove redundant 
#   gene symbols in each row such as NIPA2;NIPA2;NIPA2;NIPA2 -> NIPA2
#
#

#   File list. There should be 86 files
#   ==========================================================================
#
setwd("/gpfs/gsfs2/users/CCRBioinfo/Henry/BigQuery.Data/target_os_methylation")

txt.files <- list.files(pattern="txt", full.names=TRUE,
                path="target_os_methylation_L3_txt_format")
length(txt.files)   
length(unique(txt.files))


#   each file should has same number of lines: 385291
#   =================================================================
#
for(i in seq_along(txt.files)) system(paste0("wc -l ", txt.files[i]))


for(aFile in seq_along(txt.files))
{
    message(aFile, ": ", txt.files[aFile]);
    temp <- read.table(txt.files[aFile], header=TRUE, sep="\t", quote="")
    geneSymbols <- as.character(temp$GeneSymbols)
    for(aGene in seq_len(length(geneSymbols)))
    {
        numOfGenes <- grep(";", geneSymbols[aGene]);
        if(length(numOfGenes)>0) {
            genes <- strsplit(geneSymbols[aGene], split=";")
            genes <- unique(unlist(genes));

            if(length(genes) > 1) {
                genes <- paste0(genes, ";", collapse="");
                genes <- sub(";$", "", genes);
            }
            message(aGene, ": ", geneSymbols[aGene], " -> ", genes)           
            geneSymbols[aGene] <- genes;

        } else {message(aGene, ": ", geneSymbols[aGene]);}
    }
    
    if(aFile == 1) {
        geneList <- geneSymbols;
    } else {
        if(sum(geneList == geneSymbols) != length(geneList))
            stop("Gene list different between two files.")
    }
    temp$GeneSymbols <- geneSymbols;
    
    file.name <- sub("txt$", "csv", txt.files[aFile]);
    file.name <- sub(".{1,}/", "", file.name);
    write.table(temp, file=file.name, sep=",", quote=FALSE,
            row.names=FALSE, col.names=TRUE);
}



#   Merge all files as one
#   ==================================================================


input.files <- list.files(path="target_os_methylation_L3_csv_format", 
                pattern="csv", full.names=TRUE)
length(input.files) #   [1] 86

system("wc -l target_os_methylation_L3_csv_format/*.csv")
#   33135112 total

aFile <- 1;
message(aFile, ": ", input.files[aFile]);
bigData <- read.table(input.files[aFile], header=TRUE, sep=",", quote="")
dim(bigData)    #   [1] 385291      6

#   for data check
#   ===============================================
#
totalRows <- nrow(bigData)
totalCols <- ncol(bigData)

sampleID    <- as.character(bigData$SampleID)
reporterID  <- as.character(bigData$ReporterID)
geneSymbols <- as.character(bigData$GeneSymbols)
chromosomes <- as.character(bigData$Chromosome)
positions   <- as.numeric(bigData$Position)

for(aFile in seq_len(length(input.files))[-1])
{
    message(aFile, ": ", input.files[aFile]);
    moreData <- read.table(input.files[aFile], header=TRUE, sep=",", quote="")

    if(nrow(moreData) != totalRows || ncol(moreData) != totalCols)
        stop("Bad file dimensions!\n")

    #   sampleID in each file should be unique
    #
    moreSamplesID   <- as.character(moreData$SampleID)
    if(sum(moreSamplesID == sampleID) > 0) stop("SampleID error!")

    #   ReporterID in all files should be a same set
    #
    moreReporterID  <- as.character(moreData$ReporterID)
    if(sum(moreReporterID == reporterID) != totalRows) 
        stop("ReporterID error!")

    #   Gene symbols in all files should be a same set
    #
    moreGeneSymbols <- as.character(moreData$GeneSymbols)
    if(sum(moreGeneSymbols == geneSymbols) != totalRows)
        stop("GeneSymbols error!")

    #   Chromosome names in all files should be a same set
    #
    moreChromosomes <- as.character(moreData$Chromosome)
    if(sum(moreChromosomes == chromosomes) != totalRows)
        stop("Chromosome error!")

    #   Positions in all files should be a same set
    #
    morePositions   <- as.numeric(moreData$Position)
    if(sum(morePositions == positions) != totalRows)
        stop("Position error!")

    bigData <- rbind(bigData, moreData);
}

#   data check
#   =====================================================
dim(bigData)                    #   [1] 33135026        6
totalRows*length(input.files)   #   [1] 33135026
head(bigData)

#                     SampleID ReporterID      Signal 
#   1 TARGET-40-0A4HLD-01A-01D cg00000622 0.007846173
#   2 TARGET-40-0A4HLD-01A-01D cg00000769 0.058150113
#   3 TARGET-40-0A4HLD-01A-01D cg00001245 0.004987354
#   4 TARGET-40-0A4HLD-01A-01D cg00001349 0.489639398
#   5 TARGET-40-0A4HLD-01A-01D cg00001583 0.009334729
#   6 TARGET-40-0A4HLD-01A-01D cg00001594 0.001645875
#
#     GeneSymbols Chromosome  Position
#   1       NIPA2         15  23034447
#   2       DDX55         12 124086477
#   3      MRPS25          3  15106710
#   4        MAEL          1 166958439
#   5       NR5A2          1 200011786
#   6       ROCK2          2  11484705

save(bigData, file="target_os_methylation_level_3_all_samples.RData")
write.table(bigData, file="target_os_methylation_level_3_all_samples.csv",
        sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
#


#   data check
#   ===============================================================
system("head -n 10 target_os_methylation_level_3_all_samples.csv")

#   SampleID,ReporterID,Signal,GeneSymbols,Chromosome,Position
#   TARGET-40-0A4HLD-01A-01D,cg00000622,0.00784617293426534,NIPA2,15,23034447
#   TARGET-40-0A4HLD-01A-01D,cg00000769,0.0581501125602346,DDX55,12,124086477
#   TARGET-40-0A4HLD-01A-01D,cg00001245,0.004987353794664,MRPS25,3,15106710
#   TARGET-40-0A4HLD-01A-01D,cg00001349,0.489639398165064,MAEL,1,166958439
#   TARGET-40-0A4HLD-01A-01D,cg00001583,0.00933472931363702,NR5A2,1,200011786
#   TARGET-40-0A4HLD-01A-01D,cg00001594,0.00164587548735616,ROCK2,2,11484705
#   TARGET-40-0A4HLD-01A-01D,cg00001747,0.0176092381449123,,2,105460263
#   TARGET-40-0A4HLD-01A-01D,cg00002028,0.0206798107174625,PINK1,1,20960010
#   TARGET-40-0A4HLD-01A-01D,cg00002033,0.640620725126475,LRFN1,19,39798481

system("tail -n 10 target_os_methylation_level_3_all_samples.csv")

#   TARGET-40-PAVECB-01A-01D,ch.9.93373462R,0.000216090461952016,,9,94333641
#   TARGET-40-PAVECB-01A-01D,ch.9.941347R,0.00315605303976776,,9,77851250
#   TARGET-40-PAVECB-01A-01D,ch.9.945770F,0.0110338968252401,,9,78090914
#   TARGET-40-PAVECB-01A-01D,ch.9.97139671F,0.00775076676069893,,9,98099850
#   TARGET-40-PAVECB-01A-01D,ch.9.98936572R,0.0803823412587732,,9,99896751
#   TARGET-40-PAVECB-01A-01D,ch.9.98937537R,0.00257009542162926,,9,99897716
#   TARGET-40-PAVECB-01A-01D,ch.9.98957343R,0.00903066382316203,,9,99917522
#   TARGET-40-PAVECB-01A-01D,ch.9.98959675F,0.0953693107155211,,9,99919854
#   TARGET-40-PAVECB-01A-01D,ch.9.98989607R,0.000227294175343261,,9,99949786
#   TARGET-40-PAVECB-01A-01D,ch.9.991104F,0.000468616253400212,GNA14,9,80102978