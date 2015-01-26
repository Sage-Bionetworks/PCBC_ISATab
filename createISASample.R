## Create an ISA-Tab Sample file for all PCBC cell lines.

library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

## Setup the column names that form different parts of the ISA-Tab records

## descriptors
metadataTableIdVars <- c("PCBC_Cell_Line_Name", "C4_Cell_Line_ID", "High_Confidence_Donor_ID", "Diffname_short", "Replicate")


freetextCols <- c("Small_Molecules", "Other_Conditions_During_Reprogramming",
                  "Public_Data", "Originating_Lab", "Pubmed_ID",
                  "Culture_Conditions", "Reprogramming_Vector_Type")

# # Not working, moved to freetextCols
# freetextParams <- c("Culture_Conditions", # should be in regular, not free text, params, but ontology is not consistent
#                     "Reprogramming_Vector_Type"  # should be in regular, not free text, params, but ontology is not consistent
# )

# These columns should be turned into characteristics
characteristicsCols <- c("Donor_Life_Stage", "Race", "Ethnicity", "Gender", "Genotype",
                         "Cell_Type", "Host_Species", "Cell_Line_Type", "Tissue_of_Origin",
                         "Cell_Type_of_Origin", "Cell_Line_of_Origin")


## Get the metadata standard
metadataStandard <- synGet("syn2767699")
query <- paste("SELECT * FROM", metadataStandard$properties$id)
queryResult <- synTableQuery(query, loadResult=TRUE)
metadataStandard <- tbl_df(queryResult@values)

## Get the cell line metadata
cellLineMetadataTable <- synGet("syn2767694")
cellLineQuery <- paste("SELECT * FROM", cellLineMetadataTable$properties$id)
cellLineMetadata <- tbl_df(synTableQuery(cellLineQuery, loadResult=TRUE)@values)

# Get sample Metadata
sampleProcessMetadataTable <- synGet("syn3131513")
sampleProcessMetadataQuery <- synTableQuery(sprintf("SELECT C4_Cell_Line_ID,Diffname_short,Replicate FROM %s", sampleProcessMetadataTable$properties$id))
sampleProcessMetadata  <- sampleProcessMetadataQuery@values

mrnaAll <- tbl_df(left_join(sampleProcessMetadata, cellLineMetadata, by="C4_Cell_Line_ID"))

## Get the metadata table describing the cell lines
## Using annotations on mRNA raw data files
# colnamesToUse <- c(metadataTableIdVars, freetextCols, freetextParams, characteristicsCols)
colnamesToUse <- c(metadataTableIdVars, freetextCols, characteristicsCols)

# queryAll <- sprintf("select %s from file where benefactorId=='syn1773109' AND dataType=='mRNA'",
#                     paste(colnamesToUse, collapse=","))

# queryAll <- "select id,UID,C4_Cell_Line_ID,Diffname_short from file where benefactorId=='syn1773109' AND dataType=='mRNA'"
# qr <- synQuery(queryAll, blockSize=250)
# mrnaAll <- qr$collectAll()
# mrnaAll <- tbl_df(mrnaAll)

# # Clean column names
# colnames(mrnaAll) <- gsub("file\\.", "", colnames(mrnaAll))

# mrnaAll <- left_join(mrnaAll, cellLineMetadata, by="C4_Cell_Line_ID")

mrnaAll[mrnaAll == "N/A"] <- NA
mrnaAll[mrnaAll == ""] <- NA

metadataTable <- distinct(mrnaAll[, colnamesToUse])

# Only considering those with C4 ID's
metadataTable <- filter(metadataTable, C4_Cell_Line_ID != "", !is.na(C4_Cell_Line_ID))

## Only temporary until someone decides it's fate
metadataTable <- filter(metadataTable, C4_Cell_Line_ID != 'SC12-041')

## Also should be temporary - metadata is missing!
metadataTable <- filter(metadataTable, !(C4_Cell_Line_ID  %in% c('H9Hypox', 'SC14-066')))


# Columns that uniquely identify records
# Everything else is a measurement
metadataTableMelted <- tbl_df(melt(metadataTable, id.vars=metadataTableIdVars))


# These are free text, not in the ontology, and don't need Source or Accession descriptor columns
# Hence, they do not need to be merged with the metadataStandard and can be included as is.

# These should be Parameter Values
# paramCols <- c("Reprogramming_Vector_Type")

# Filter the melted metadata table for each of the sets of columns
metadataTableCharacteristics <- filter(metadataTableMelted, variable %in% characteristicsCols)
metadataTableFreetext <- filter(metadataTableMelted, variable %in% freetextCols)
# metadataTableFreetextParams <- filter(metadataTableMelted, variable %in% freetextParams)
# metadataTableParams <- filter(metadataTableMelted, variable %in% paramCols)

# Combine these tables with the metadata standard to get accessions, ontology information
mergedTableCharacteristics <- merge(metadataTableCharacteristics, metadataStandard, 
                                    by.x=c("variable", "value"),
                                    by.y=c("Category", "Metadata_Term"))

# mergedTableParams <- merge(metadataTableParams, metadataStandard, 
#                            by.x=c("variable", "value"),
#                            by.y=c("Category", "Metadata_Term"))

# Do stuff with the characteristics columns
# These will end up with prefix of "Characteristics[]" in the column name

# Get value columns
mergedTableCharacteristicsTermValues <- dcast(mergedTableCharacteristics,
                                              C4_Cell_Line_ID + Replicate + Diffname_short ~ variable,
                                              value.var="value")

rownames(mergedTableCharacteristicsTermValues) <- with(mergedTableCharacteristicsTermValues,
                                                       paste(C4_Cell_Line_ID,
                                                             Replicate,
                                                             Diffname_short,
                                                             sep=" "))
                                                       
mergedTableCharacteristicsTermValues$C4_Cell_Line_ID <- NULL
mergedTableCharacteristicsTermValues$Diffname_short <- NULL
mergedTableCharacteristicsTermValues$Replicate <- NULL

colnames(mergedTableCharacteristicsTermValues) <- paste("Characteristics[", 
                                                        colnames(mergedTableCharacteristicsTermValues), 
                                                        "]", sep="")

# Get Source columns
mergedTableCharacteristicsTermSource <- dcast(mergedTableCharacteristics,
                                              C4_Cell_Line_ID + Replicate + Diffname_short ~ variable,
                                              value.var="Source_Ontology")


rownames(mergedTableCharacteristicsTermSource) <- with(mergedTableCharacteristicsTermSource,
                                                       paste(C4_Cell_Line_ID,
                                                             Replicate,
                                                             Diffname_short,
                                                             sep=" "))

mergedTableCharacteristicsTermSource$C4_Cell_Line_ID <- NULL
mergedTableCharacteristicsTermSource$Diffname_short <- NULL
mergedTableCharacteristicsTermSource$Replicate <- NULL

colnames(mergedTableCharacteristicsTermSource) <- paste("Term Source REF[", 
                                                        colnames(mergedTableCharacteristicsTermSource),
                                                        "]", sep="")

# Get Accession columns
mergedTableCharacteristicsTermAccessions <- dcast(mergedTableCharacteristics,
                                                  C4_Cell_Line_ID + Replicate + Diffname_short ~ variable,
                                                  value.var="URL_Xref")

rownames(mergedTableCharacteristicsTermAccessions) <- with(mergedTableCharacteristicsTermAccessions,
                                                           paste(C4_Cell_Line_ID,
                                                                 Replicate,
                                                                 Diffname_short,
                                                                 sep=" "))

mergedTableCharacteristicsTermAccessions$C4_Cell_Line_ID <- NULL
mergedTableCharacteristicsTermAccessions$Diffname_short <- NULL
mergedTableCharacteristicsTermAccessions$Replicate <- NULL

colnames(mergedTableCharacteristicsTermAccessions) <- paste("Term Accession Number[", 
                                             colnames(mergedTableCharacteristicsTermAccessions), 
                                             "]", sep="")

# Combine the value, source, and accession columns together
combinedCharacteristicsCols <- cbind(mergedTableCharacteristicsTermValues,
                                     mergedTableCharacteristicsTermSource,
                                     mergedTableCharacteristicsTermAccessions)

# Reorder columns so for each characteristic, the source and accession follow it.
reorderedCols <- c(sapply(1:(ncol(combinedCharacteristicsCols)/3), 
                          function(x) seq(from=x, 
                                          to=ncol(combinedCharacteristicsCols),
                                          by=ncol(combinedCharacteristicsCols) / 3)))

combinedCharacteristicsCols <- combinedCharacteristicsCols[, reorderedCols]

# # Parameter columns
# 
# # These columns should be turned into parameters
# 
# # Get value columns
# mergedTableParamsTermValues <- dcast(mergedTableParams,
#                                      C4_Cell_Line_ID ~ variable,
#                                      value.var="value")
# 
# rownames(mergedTableParamsTermValues) <- mergedTableParamsTermValues$C4_Cell_Line_ID
# mergedTableParamsTermValues$C4_Cell_Line_ID <- NULL
# colnames(mergedTableParamsTermValues) <- paste("Parameter Value[", 
#                                          colnames(mergedTableParamsTermValues), 
#                                          "]", sep="")
# 
# # Get Source columns
# mergedTableParamsTermSource <- dcast(mergedTableParams, 
#                                      C4_Cell_Line_ID ~ variable,
#                                      value.var="Source_Ontology")
#
# rownames(mergedTableParamsTermSource) <- mergedTableParamsTermSource$C4_Cell_Line_ID
# mergedTableParamsTermSource$C4_Cell_Line_ID <- NULL
# colnames(mergedTableParamsTermSource) <- paste("Term Source REF[", 
#                                          colnames(mergedTableParamsTermSource), 
#                                          "]", sep="")
# 
# # Get Accession columns
# mergedTableParamsTermAccessions <- dcast(mergedTableParams,
#                                          C4_Cell_Line_ID ~ variable,
#                                          value.var="URL_Xref")
# 
# rownames(mergedTableParamsTermAccessions) <- mergedTableParamsTermAccessions$C4_Cell_Line_ID
# mergedTableParamsTermAccessions$C4_Cell_Line_ID <- NULL
# colnames(mergedTableParamsTermAccessions) <- paste("Term Accession Number[", 
#                                              colnames(mergedTableParamsTermAccessions), 
#                                              "]", sep="")
# 
# combinedParamsCols <- cbind(mergedTableParamsTermValues,
#                                      mergedTableParamsTermSource,
#                                      mergedTableParamsTermAccessions)
# 
# reorderedCols <- c(sapply(1:(ncol(combinedParamsCols)/3), 
#                           function(x) seq(from=x, 
#                                           to=ncol(combinedParamsCols),
#                                           by=ncol(combinedParamsCols) / 3)))
# 
# combinedParamsCols <- combinedParamsCols[, reorderedCols]

# Free text columns
# These do not need to be merged with the metadata standard
# There is only the value column, so no reordering is necessary either!

# Free text as characteristics
freetextTermValues <- dcast(metadataTableFreetext, 
                            C4_Cell_Line_ID + Replicate + Diffname_short ~ variable,
                            value.var="value")

rownames(freetextTermValues) <- with(freetextTermValues,
                                     paste(C4_Cell_Line_ID,
                                           Replicate,
                                           Diffname_short,
                                           sep=" "))

freetextTermValues$C4_Cell_Line_ID <- NULL
freetextTermValues$Diffname_short <- NULL
freetextTermValues$Replicate <- NULL

colnames(freetextTermValues) <- paste("Characteristics[", 
                                      colnames(freetextTermValues), 
                                      "]", sep="")

# # Free text as parameters
# freetextParamValues <- dcast(metadataTableFreetextParams, 
#                             C4_Cell_Line_ID + Replicate + Diffname_short ~ variable,
#                             value.var="value")
# 
# rownames(freetextParamValues) <- with(freetextParamValues,
#                                       paste(C4_Cell_Line_ID,
#                                             Replicate,
#                                             Diffname_short,
#                                             sep=" "))
# 
# freetextParamValues$C4_Cell_Line_ID <- NULL
# freetextParamValues$Diffname_short <- NULL
# freetextParamValues$Replicate <- NULL
# 
# colnames(freetextParamValues) <- paste("Parameter Value[", 
#                                       colnames(freetextParamValues), 
#                                       "]", sep="")

# When params are working again, use this!
#combinedCols <- cbind(combinedCharacteristicsCols, combinedParamsCols, freetextTermValues, freetextParamValues)

# only use these
useRows <- rownames(combinedCharacteristicsCols)
combinedCols <- cbind(combinedCharacteristicsCols[useRows, ], freetextTermValues[useRows, ])
# combinedCols <- cbind(combinedCols[useRows, ], freetextParamValues[useRows, ])

# Add back the row names, this will become the Source Name column
combinedCols <- cbind(rownames(combinedCols), combinedCols)

# For easier output, convert to a tbldf; this allows duplicated column names
# Why ISA-Tab would allow duplicated column names is beyond me...

combinedCols <- tbl_df(combinedCols)

# Fix the column names for ISA-Tab format
colnames(combinedCols)[1] <- "Source Name"
colnames(combinedCols) <- gsub(pattern="Term Source REF.*", "Term Source REF", colnames(combinedCols))
colnames(combinedCols) <- gsub(pattern="Term Accession Number.*", "Term Accession Number", colnames(combinedCols))

## Output, including saving and uploading the script and output file to Synapse.
## Records provenance to include the Synapse tables queried to create this file.

write.table(combinedCols, file="s_PCBC.txt", row.names=FALSE, quote=1:ncol(combinedCols), sep="\t")

synfileMe  <- File(path="createISASample.R", parentId="syn2814512")
synfileMe <- synStore(synfileMe)

synfile  <- File(path="s_PCBC.txt", name="PCBC Study ISA-Tab", parentId="syn2814512")
synfile <- synStore(synfile,
                    executed=synfileMe$properties$id,
                    used=c(metadataStandard$properties$id, cellLineMetadataTable$properties$id))
