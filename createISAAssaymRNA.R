## Create an mRNA ISA-Tab Assay file for all PCBC cell lines.

library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

rnaSeqMetadataTable <- synGet("syn3104413")
rnaSeqMetadataQuery <- synTableQuery(sprintf("SELECT UID,replicate FROM %s", rnaSeqMetadataTable$properties$id))
rnaSeqMetadata <- rnaSeqMetadataQuery@values

## Generally useful columns
generalCols <- c("UID", "C4_Cell_Line_ID", "dataType", "fileType", "Diffname_short")

fileTypes <- c("bam", "fastq")

## Get all the files
queryAll <- "select id,name,UID,C4_Cell_Line_ID,Diffname_short,fileType,fileSubType,bamType,fpkmType from file where benefactorId=='syn1773109' AND dataType=='mRNA'"
qr <- synQuery(queryAll, blockSize = 250)
mrnaAll <- qr$collectAll()
mrnaAll <- tbl_df(mrnaAll)
colnames(mrnaAll) <- gsub("file\\.", "", colnames(mrnaAll))

mrnaAll <- merge(mrnaAll, rnaSeqMetadata, by="UID")

mrnaAll <- filter(mrnaAll, Diffname_short %in% c("DE", "MESO-5", "MESO-15", 
                                                 "MESO-30", "SC", "EB", "ECTO"),
                  public, pass_qc, !exclude)

# The raw data file
mrnaFastq <- filter(mrnaAll, fileType=="fastq", public)

# Derived files
mrnaBam <- filter(mrnaAll, fileType=="bam", fileSubType=="mapped")
mrnaFpkmGene <- filter(mrnaAll, fileType=="fpkm", fileSubType == "genes")
mrnaFpkmIso <- filter(mrnaAll, fileType=="fpkm", fileSubType == "isoform")

# ## Convert separate file types to a subtype attribute
# mrnaBam$subType  <- "mapped"
mrnaFastq$fileSubType <- "fastq"
# mrnaFpkmGene$subType <- "genes"
# mrnaFpkmIso$subType <- "isoform"

mrnaUse <- rbind(mrnaBam, mrnaFastq, mrnaFpkmGene, mrnaFpkmIso)

mrnaUse <- filter(mrnaUse,
                  # !is.na(replicate), replicate != "",
                  !is.na(Diffname_short), Diffname_short != "",
                  !is.na(C4_Cell_Line_ID), C4_Cell_Line_ID != ""
                  )

mrnaUse$sample.name <- with(mrnaUse, paste(C4_Cell_Line_ID, replicate, Diffname_short))

castToISA <- dcast(mrnaUse,
                   sample.name + UID ~ fileType + fileSubType,
                   value.var="name",
                   fun.aggregate=function(x) paste(x, collapse=","))

# castToISA <- filter(castToISA, bam_mapped != "", fastq_fastq != "", fpkm_genes != "", fpkm_isoform != "")

castToISA$fastq.comment <- paste("Synapse ID", mrnaUse$id[match(castToISA$fastq_fastq, mrnaUse$name)])
castToISA$fastq.comment[castToISA$fastq.comment == "Synapse ID NA"] <- ""

castToISA$bam.trans.name <- "Alignment"
castToISA$bam.norm.name <- "NA"
castToISA$bam.comment <- paste("Synapse ID", mrnaUse$id[match(castToISA$bam_mapped, mrnaUse$name)])
castToISA$bam.comment[castToISA$bam.comment == "Synapse ID NA"] <- ""

castToISA$fpkm_genes.trans.name <- "Cufflinks Genes"
castToISA$fpkm_genes.norm.name <- "NA"
castToISA$fpkm_genes.comment <- paste("Synapse ID", mrnaUse$id[match(castToISA$fpkm_genes, mrnaUse$name)])
castToISA$fpkm_genes.comment[castToISA$fpkm_genes.comment == "Synapse ID NA"] <- ""

castToISA$fpkm_isoform.trans.name <- "Cufflinks Isoforms"
castToISA$fpkm_isoform.norm.name <- "NA"
castToISA$fpkm_isoform.comment <- paste("Synapse ID", mrnaUse$id[match(castToISA$fpkm_isoform, mrnaUse$name)])
castToISA$fpkm_isoform.comment[castToISA$fpkm_isoform.comment == "Synapse ID NA"] <- ""

castToISA <- tbl_df(castToISA)

attrcols <- c("sample.name", "UID")
assaycols <- c("fastq_fastq", "fastq.comment",
               "bam.norm.name", "bam.trans.name", "bam_mapped", "bam.comment",
               "fpkm_genes.norm.name", "fpkm_genes.trans.name", "fpkm_genes", "fpkm_genes.comment",
               "fpkm_isoform.norm.name","fpkm_isoform.trans.name", "fpkm_isoform", "fpkm_isoform.comment")

colorder <- c(attrcols, assaycols)

cols.rename <- c("Sample Name", "Assay Name",
                 "Raw Data File", "Comment",
                 'Normalization Name', 'Data Transformation Name', 'Derived Data File', 'Comment',
                 'Normalization Name', 'Data Transformation Name', 'Derived Data File', 'Comment',
                 'Normalization Name', 'Data Transformation Name', 'Derived Data File', 'Comment')

castToISA <- castToISA[, colorder]

names(castToISA) <- cols.rename

# Add summarized file
summarizedFile <- synGet("syn2701943", downloadFile=FALSE)

castToISA$"Normalization Name[Summarize]" <- "NA"
castToISA$"Data Transformation Name[Summarize]" <- "Summarized data"
castToISA$"Derived Data File[Summarize]" <- summarizedFile$properties$name
castToISA$"Comment[Summarize]" <- paste("Synapse ID", summarizedFile$properties$id)

# Fix the duplicated column names
colnames(castToISA) <- gsub(pattern="Comment.*", "Comment", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Normalization Name.*", "Normalization Name", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Data Transformation Name.*", "Data Transformation Name", colnames(castToISA))
colnames(castToISA) <- gsub(pattern="Derived Data File.*", "Derived Data File", colnames(castToISA))

write.table(castToISA, file="a_mRNA_PCBC.txt", row.names=FALSE, quote=1:ncol(castToISA), sep="\t")

synfileMe  <- File(path="createISAAssaymRNA.R", parentId="syn2814512")
synfileMe <- synStore(synfileMe)

synfile  <- File(path="a_mRNA_PCBC.txt", name="PCBC mRNA Assay ISA-Tab", parentId="syn2814512")
synfile <- synStore(synfile, 
                    executed=synfileMe$properties$id,
                    used=c(rnaSeqMetadataTable$properties$id))
