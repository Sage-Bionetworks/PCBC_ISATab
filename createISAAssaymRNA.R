## Create an mRNA ISA-Tab Sample file for all PCBC cell lines.

library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

## Generally useful columns
generalCols <- c("UID", "C4_Cell_Line_ID", "dataType", "fileType", "Cell_Line_of_Origin",
                 "Diffname_short", "Public_Data")

fileTypes <- c("bam", "fastq")

## Get all the files
queryAll <- "select * from file where benefactorId=='syn1773109' AND dataType=='mRNA'"
qr <- synQuery(queryAll, blockSize = 50)
mrnaAll <- qr$collectAll()
mrnaAll <- tbl_df(mrnaAll)
colnames(mrnaAll) <- gsub("file\\.", "", colnames(mrnaAll))

mrnaAll$fpkmType <- gsub(".*_(.*)\\.fpkm_tracking", "\\1", mrnaAll$name)


mrnaBam <- filter(mrnaAll, fileType=="bam", bamType=="mapped")
mrnaFastq <- filter(mrnaAll, fileType=="fastq")
mrnaFpkm <- filter(mrnaAll, fileType=="fpkm", fpkmType == "genes")

mrnaUse <- rbind(mrnaBam, mrnaFastq, mrnaFpkm)

intersectUIDs <- intersect(mrnaBam$UID, mrnaFastq$UID)
missingFastqUIDs <- setdiff(mrnaBam$UID, mrnaFastq$UID)
missingBamUIDs <- setdiff(mrnaFastq$UID, mrnaBam$UID)

length(intersectUIDs)
length(missingBamUIDs)
length(missingFastqUIDs)

head(mrnaFastq[, generalCols])
filter(mrnaUse, UID %in% missingBamUIDs[1:5])[, generalCols]

filter(mrnaUse, C4_Cell_Line_ID=="H9", Diffname_short=="EB")[, generalCols]

meltByFileType <- melt(mrnaUse, id.vars=c("id", "UID", "C4_Cell_Line_ID", "Diffname_short"),
                       measure.vars="fileType")

mRNABinaryStatus <- dcast(meltByFileType,  UID ~ value)[, c("UID", "fastq", "bam", "fpkm")]
write.csv(mRNABinaryStatus, file="./mRNABinaryStatus.csv")
