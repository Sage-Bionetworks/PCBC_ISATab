## Create an miRNA ISA-Tab Assay file for all PCBC cell lines.

library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

sampleProcessMetadataTable <- synGet("syn3131513")
sampleProcessMetadataQuery <- synTableQuery(sprintf("SELECT * FROM %s WHERE DateDNASubmittedToCore<>'N/A'", sampleProcessMetadataTable$properties$id))
sampleProcessMetadata <- sampleProcessMetadataQuery@values

# Only include those that were actually submitted to the core and not old or test
# The old/test was breaking the matching, because no other samples have a 'not old' designation...
sampleProcessMetadata <- filter(sampleProcessMetadata,
                                DateDNASubmittedToCore != "N/A",
                                !SampleTypeNotes %in% c("old lot", "test lot"))

metadataTable <- synGet('syn3105812')
metadataQuery <- synTableQuery(sprintf("SELECT File,Sample,biologicalSampleName,Channel,TechnicalReplicate FROM %s", metadataTable$properties$id))
metadata <- metadataQuery@values

metadata <- merge(metadata, sampleProcessMetadata, by="biologicalSampleName")

## Get all the files
queryAll <- "select id,name,UID,fileType from file where benefactorId=='syn1773109' AND dataType=='methylation' AND fileType=='idat'"
qr <- synQuery(queryAll, blockSize = 250)
fileData <- qr$collectAll()
# fileData <- tbl_df(fileData)
colnames(fileData) <- gsub("file\\.", "", colnames(fileData))

fileData <- merge(fileData, metadata, by.x="name", by.y="File")
fileData$sample.name <- with(fileData, paste(C4_Cell_Line_ID, Replicate, Diffname_short))
fileData$assay.name <- with(fileData, paste(C4_Cell_Line_ID, Replicate, Diffname_short, TechnicalReplicate))

# Only EB and SC for this paper
fileData <- filter(fileData, Diffname_short %in% c("EB", "SC"))

# filter samples if they were redone based on the replicate notes
filterRedone <- function(x) {
  if (any(x$ReplicateNotes == "R")) {
    subset(x, ReplicateNotes == "R")  
  }
  else {
    x
  }
}

# select the samples that were redone, if so
fileDataReduced <- ddply(fileData, 
                         .(Sample), 
                         filterRedone)

castToISA <- dcast(fileDataReduced, 
                   sample.name + assay.name ~ Channel, 
                   value.var="name", fun.aggregate=str_c, collapse=", ")

dcast(fileDataReduced, 
                   sample.name + assay.name ~ Channel, 
                   value.var="name")

castToISA <- tbl_df(castToISA)

castToISA$Red.comment <- paste("Synapse ID", fileData$id[match(castToISA$Red, fileData$name)])
castToISA$Red.comment[castToISA$Red.comment == "Synapse ID NA"] <- ""

castToISA$Grn.comment <- paste("Synapse ID", fileData$id[match(castToISA$Grn, fileData$name)])
castToISA$Grn.comment[castToISA$Grn.comment == "Synapse ID NA"] <- ""

attrcols <- c("sample.name", "assay.name")
assaycols <- c("Red", "Red.comment",
               "Grn", "Grn.comment")

colorder <- c(attrcols, assaycols)

cols.rename <- c("Sample Name", "Assay Name",
                 "Raw Data File", "Comment",
                 "Raw Data File", "Comment")

castToISA <- castToISA[, colorder]

names(castToISA) <- cols.rename

# Add summarized file
summarizedFile <- synGet("syn2233188", downloadFile=FALSE)

castToISA$"Normalization Name[Summarize]" <- "NA"
castToISA$"Data Transformation Name[Summarize]" <- "Summarized data"
castToISA$"Derived Data File[Summarize]" <- summarizedFile$properties$name
castToISA$"Comment[Summarize]" <- paste("Synapse ID", summarizedFile$properties$id)

# Fix the duplicated column names
colnames(castToISA) <- gsub(pattern="Comment.*", "Comment", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Normalization Name.*", "Normalization Name", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Data Transformation Name.*", "Data Transformation Name", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Derived Data File.*", "Derived Data File", colnames(castToISA))

write.table(castToISA, file="a_methylation_PCBC.txt", row.names=FALSE, quote=1:ncol(castToISA), sep="\t")

synfileMe  <- File(path="createISAAssayMethylation.R", parentId="syn2814512")
synfileMe <- synStore(synfileMe)

synfile  <- File(path="a_methylation_PCBC.txt", name="PCBC Methylation Assay ISA-Tab", parentId="syn2814512")

synfile <- synStore(synfile, 
                    executed=synfileMe$properties$id,
                    used=c(sampleProcessMetadataTable$properties$id, metadataTable$properties$id))
