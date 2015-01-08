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

# The raw data file
mrnaFastq <- filter(mrnaAll, fileType=="fastq")

# Derived files
mrnaBam <- filter(mrnaAll, fileType=="bam", bamType=="mapped")
mrnaFpkmGene <- filter(mrnaAll, fileType=="fpkm", fpkmType == "genes")
mrnaFpkmIso <- filter(mrnaAll, fileType=="fpkm", fpkmType == "isoform")

## Convert separate file types to a subtype attribute
mrnaBam$subType  <- "mapped"
mrnaFastq$subType <- "fastq"
mrnaFpkmGene$subType <- "genes"
mrnaFpkmIso$subType <- "isoform"

mrnaUse <- rbind(mrnaBam, mrnaFastq, mrnaFpkmGene, mrnaFpkmIso)

mrnaUse$sample.name <- with(mrnaUse, paste(C4_Cell_Line_ID, Diffname_short))

castToISA <- dcast(mrnaUse,
                   sample.name + UID ~ fileType + subType,
                   value.var="name",
                   fun.aggregate=function(x) paste(x, collapse=","))

# castToISA <- filter(castToISA, bam_mapped != "", fastq_fastq != "", fpkm_genes != "", fpkm_isoform != "")

castToISA$fastq.comment <- paste("Synapse ID", mrnaUse$id[match(castToISA$fastq_fastq, mrnaUse$name)])

castToISA$bam.trans.name <- "Alignment"
castToISA$bam.comment <- paste("Synapse ID", mrnaUse$id[match(castToISA$bam_mapped, mrnaUse$name)])

castToISA$fpkm_genes.trans.name <- "Cufflinks Genes"
castToISA$fpkm_genes.comment <- paste("Synapse ID", mrnaUse$id[match(castToISA$fpkm_genes, mrnaUse$name)])

castToISA$fpkm_isoform.trans.name <- "Cufflinks Isoforms"
castToISA$fpkm_isoform.comment <- paste("Synapse", mrnaUse$id[match(castToISA$fpkm_isoform, mrnaUse$name)])

castToISA <- tbl_df(castToISA)

attrcols <- c("sample.name", "UID")
assaycols <- c("fastq_fastq", "fastq.comment",
               "bam.trans.name", "bam_mapped", "bam.comment",
               "fpkm_genes.trans.name", "fpkm_genes", "fpkm_genes.comment",
               "fpkm_isoform.trans.name", "fpkm_isoform", "fpkm_isoform.comment")

colorder <- c(attrcols, assaycols)

cols.rename <- c("Sample Name", "Assay Name",
                 "Raw Data File", "Comment",
                 'Data Transformation Name', 'Derived Data File', 'Comment',
                 'Data Transformation Name', 'Derived Data File', 'Comment',
                 'Data Transformation Name', 'Derived Data File', 'Comment')

castToISA <- castToISA[, colorder]

names(castToISA) <- cols.rename

# Fix the duplicated column names
colnames(castToISA) <- gsub(pattern="Comment.*", "Comment", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Data Transformation Name.*", "Data Transformation Name", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Derived Data File.*", "Derived Data File", colnames(castToISA))

write.csv(castToISA, file="a_PCBC_all.csv", row.names=FALSE)

synfileMe  <- File(path="createISAAssaymRNA.R", parentId="syn2814512")
synfileMe <- synStore(synfileMe)

synfile  <- File(path="a_PCBC_all.csv", parentId="syn2814512")
synfile <- synStore(synfile, 
                    executed=synfileMe$properties$id)
