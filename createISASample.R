## Create an ISA-Tab Sample file for all PCBC cell lines.

library(dplyr)
library(reshape2)
library(synapseClient)
synapseLogin()

## Get the metadata standard
metadataStandardSynId <- "syn2767699"
query <- paste("SELECT * FROM", metadataStandardSynId)
queryResult <- synTableQuery(query, loadResult=TRUE)
metadataStandard <- tbl_df(queryResult@values)

## Get the metadata table describing the cell lines
metadataTableSynId <- "syn2767694"
query <- paste("SELECT * FROM", metadataTableSynId)
queryResult <- synTableQuery(query, loadResult=TRUE)

metadataTable <- tbl_df(queryResult@values)

# Only considering those with C4 ID's
metadataTable <- filter(metadataTable, C4_Cell_Line_ID != "")

## Only temporary until someone decides it's fate
metadataTable <- filter(metadataTable, C4_Cell_Line_ID != 'SC12-041')

# Columns that uniquely identify records
# Everything else is a measurement

metadataTableIdVars <- c("PCBC_Cell_Line_Name", "C4_Cell_Line_ID", "High_Confidence_Donor_ID")
metadataTableMelted <- tbl_df(melt(metadataTable, id.vars=metadataTableIdVars))


# These are free text, not in the ontology, and don't need Source or Accession descriptor columns
# Hence, they do not need to be merged with the metadataStandard and can be included as is.

freetextCols <- c("Small_Molecules", "Other_Conditions_During_Reprogramming",
                  "Public_Data", "Originating_Lab", "Pubmed_ID")

freetextParams <- c("Culture_Conditions", # should be in regular, not free text, params, but ontology is not consistent
                    "Reprogramming_Vector_Type"  # should be in regular, not free text, params, but ontology is not consistent
                    )

# These columns should be turned into characteristics
characteristicsCols <- c("Donor_Life_Stage", "Race", "Ethnicity", "Gender", "Genotype",
                         "Cell_Type", "Host_Species", "Cell_Line_Type", "Tissue_of_Origin",
                         "Cell_Type_of_Origin", "Cell_Line_of_Origin")

# These should be Parameter Values
# paramCols <- c("Reprogramming_Vector_Type")

# Filter the melted metadata table for each of the sets of columns
metadataTableCharacteristics <- filter(metadataTableMelted, variable %in% characteristicsCols)
metadataTableFreetext <- filter(metadataTableMelted, variable %in% freetextCols)
metadataTableFreetextParams <- filter(metadataTableMelted, variable %in% freetextParams)
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
                                              C4_Cell_Line_ID ~ variable,
                                              value.var="value")

rownames(mergedTableCharacteristicsTermValues) <- mergedTableCharacteristicsTermValues$C4_Cell_Line_ID
mergedTableCharacteristicsTermValues$C4_Cell_Line_ID <- NULL
colnames(mergedTableCharacteristicsTermValues) <- paste("Characteristics[", 
                                         colnames(mergedTableCharacteristicsTermValues), 
                                         "]", sep="")

# Get Source columns
mergedTableCharacteristicsTermSource <- dcast(mergedTableCharacteristics,
                                              C4_Cell_Line_ID ~ variable,
                                              value.var="Source_Ontology")

rownames(mergedTableCharacteristicsTermSource) <- mergedTableCharacteristicsTermSource$C4_Cell_Line_ID
mergedTableCharacteristicsTermSource$C4_Cell_Line_ID <- NULL
colnames(mergedTableCharacteristicsTermSource) <- paste("Term Source REF[", 
                                             colnames(mergedTableCharacteristicsTermSource), 
                                             "]", sep="")

# Get Accession columns
mergedTableCharacteristicsTermAccessions <- dcast(mergedTableCharacteristics,
                                                  C4_Cell_Line_ID ~ variable,
                                                  value.var="URL_Xref")

rownames(mergedTableCharacteristicsTermAccessions) <- mergedTableCharacteristicsTermAccessions$C4_Cell_Line_ID
mergedTableCharacteristicsTermAccessions$C4_Cell_Line_ID <- NULL
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
                            C4_Cell_Line_ID ~ variable,
                            value.var="value")

freetextTermValues$C4_Cell_Line_ID <- NULL
colnames(freetextTermValues) <- paste("Characteristics[", 
                                      colnames(freetextTermValues), 
                                      "]", sep="")

# Free text as parameters
freetextParamValues <- dcast(metadataTableFreetextParams, 
                            C4_Cell_Line_ID ~ variable,
                            value.var="value")

freetextParamValues$C4_Cell_Line_ID <- NULL
colnames(freetextParamValues) <- paste("Parameter Value[", 
                                      colnames(freetextParamValues), 
                                      "]", sep="")

# When params are working again, use this!
#combinedCols <- cbind(combinedCharacteristicsCols, combinedParamsCols, freetextTermValues, freetextParamValues)

combinedCols <- cbind(combinedCharacteristicsCols, freetextTermValues, freetextParamValues)

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

write.csv(combinedCols, file="s_PCBC_all.csv", row.names=FALSE)

synfileMe  <- File(path="createISASample.R", parentId="syn2814512")
synfileMe <- synStore(synfileMe)

synfile  <- File(path="s_PCBC_all.csv", parentId="syn2814512")
synfile <- synStore(synfile, 
                    executed=synfileMe$properties$id,
                    used=c(metadataStandardSynId, metadataTableSynId))
