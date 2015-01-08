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

mrnaBam <- filter(mrnaAll, fileType=="bam", bamType=="mapped")
mrnaFastq <- filter(mrnaAll, fileType=="fastq")
mrnaFpkm <- filter(mrnaAll, fileType=="fpkm", fpkmType == "genes")
# mrnaFpkmIso <- filter(mrnaAll, fileType=="fpkm", fpkmType == "isoforms")

mrnaUse <- rbind(mrnaBam, mrnaFastq, mrnaFpkm)
castToISA <- dcast(mrnaUse, UID + C4_Cell_Line_ID + Diffname_short ~ fileType,
                   value.var="name",
                   fun.aggregate=function(x) paste(x, collapse=","))

castToISA <- filter(castToISA, bam != "", fastq != "", fpkm != "")

castToISA$fastq.source <- "Synapse"
castToISA$bam.source <- "Synapse"
castToISA$fpkm.source <- "Synapse"

castToISA$fastq.id <- mrnaUse$id[match(castToISA$fastq, mrnaUse$name)]
castToISA$bam.id <- mrnaUse$id[match(castToISA$bam, mrnaUse$name)]
castToISA$fpkm.id <- mrnaUse$id[match(castToISA$fpkm, mrnaUse$name)]

