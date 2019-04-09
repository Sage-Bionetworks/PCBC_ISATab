## Create an miRNA ISA-Tab Assay file for all PCBC cell lines.

library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

sampleProcessMetadataTable <- synGet("syn3131513")
sampleProcessMetadataQuery <- synTableQuery(sprintf("SELECT * FROM %s", sampleProcessMetadataTable$properties$id))
sampleProcessMetadata <- sampleProcessMetadataQuery@values

metadataTable <- synGet("syn3105814")
metadataQuery <- synTableQuery(sprintf("SELECT UID,biologicalSampleName,pass_qc FROM %s", metadataTable$properties$id))
metadata <- metadataQuery@values

metadata <- subset(metadata, pass_qc == TRUE)

metadata <- merge(metadata, sampleProcessMetadata, by="biologicalSampleName")

## Get all the files
queryAll <- "select id,name,UID,fileType from file where benefactorId=='syn1773109' AND dataType=='miRNA'"
qr <- synQuery(queryAll, blockSize = 250)
fileData <- qr$collectAll()
fileData <- tbl_df(fileData)
colnames(fileData) <- gsub("file\\.", "", colnames(fileData))

fileData <- merge(fileData, metadata, by="UID")

fileData <- filter(fileData, Diffname_short %in% c("EB", "SC"))

fileData$sample.name <- with(fileData, paste(C4_Cell_Line_ID, Replicate, Diffname_short))

castToISA <- dcast(fileData, 
                   UID + sample.name ~ fileType,
                   value.var="name")

castToISA$fastq.comment <- paste("Synapse ID", fileData$id[match(castToISA$fastq, fileData$name)])
castToISA$fastq.comment[castToISA$fastq.comment == "Synapse ID NA"] <- ""

castToISA$align.trans.name <- "mirExpress Alignment"
castToISA$align.norm.name <- "NA"
castToISA$align.comment <- paste("Synapse ID", fileData$id[match(castToISA$align, fileData$name)])
castToISA$align.comment[castToISA$align.comment == "Synapse ID NA"] <- ""

castToISA$expr.trans.name <- "mirExpress Expression"
castToISA$expr.norm.name <- "NA"
castToISA$expr.comment <- paste("Synapse ID", fileData$id[match(castToISA$expr, fileData$name)])
castToISA$expr.comment[castToISA$expr.comment == "Synapse ID NA"] <- ""

castToISA <- tbl_df(castToISA)

attrcols <- c("sample.name", "UID")
assaycols <- c("fastq", "fastq.comment",
               "align.norm.name", "align.trans.name", "align", "align.comment",
               "expr.norm.name", "expr.trans.name", "expr", "expr.comment")

colorder <- c(attrcols, assaycols)

cols.rename <- c("Sample Name", "Assay Name",
                 "Raw Data File", "Comment",
                 'Normalization Name', 'Data Transformation Name', 'Derived Data File', 'Comment',
                 'Normalization Name', 'Data Transformation Name', 'Derived Data File', 'Comment')

castToISA <- castToISA[, colorder]

names(castToISA) <- cols.rename

# Add summarized file
summarizedFile <- synGet("syn2701942", downloadFile=FALSE)

castToISA$"Normalization Name[Summarize]" <- "NA"
castToISA$"Data Transformation Name[Summarize]" <- "Summarized data"
castToISA$"Derived Data File[Summarize]" <- summarizedFile$properties$name
castToISA$"Comment[Summarize]" <- paste("Synapse ID", summarizedFile$properties$id)

# Fix the duplicated column names
colnames(castToISA) <- gsub(pattern="Comment.*", "Comment", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Normalization Name.*", "Normalization Name", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Data Transformation Name.*", "Data Transformation Name", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Derived Data File.*", "Derived Data File", colnames(castToISA))

write.table(castToISA, file="a_miRNA_PCBC.txt", row.names=FALSE, quote=1:ncol(castToISA), sep="\t")

synfileMe  <- File(path="createISAAssaymiRNA.R", parentId="syn2814512")
synfileMe <- synStore(synfileMe)

synfile  <- File(path="a_miRNA_PCBC.txt", name="PCBC miRNA Assay ISA-Tab", parentId="syn2814512")

synfile <- synStore(synfile, 
                    executed=synfileMe$properties$id,
                    used=c(sampleProcessMetadataTable$properties$id, metadataTable$properties$id))
