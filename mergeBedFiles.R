library(data.table)

# Function to merge two BED files
# Filenames are auto-created for a given chromosome

chr <- "22"

prevDataTableDir <- "TP73_datatable_20240303"
outputDataDir <- "."
updateDataDir <- "."



mergeTwoBedFiles <- function(chr) {
    file1.name <- file.path(prevDataTableDir,paste0("TP73_datatable_",chr,".bed.gz"))
    file2.name <- file.path(updateDataDir,paste0("TP73_MA0861.1_bidirect_",chr,".combined.bed"))
    outputFile.name <- file.path(outputDataDir,paste0("TP73_datatable_",chr,".bed.gz"))

    if (file.exists(outputFile.name)) {
        cat("I: Found existing file for chromosome", chr, "...\n")
        next
    }

    # Read both BED files
    bed1 <- fread(file1.name, sep = "\t", header = TRUE)
    bed2 <- fread(file2.name, sep = "\t", header = TRUE)
    bed2.sub <- bed2[1:nrow(bed1),]
    # Check that first 3 columns (chrom, start, end) are identical
    if (!all.equal(bed1[, 1:3], bed2.sub[, 1:3])) {
        stop(paste("The first three columns (chrom, start, end) do not match for chromosome", chr))
    }
    # Combine the first three columns with the remaining columns from both files
    combined <- cbind(bed2.sub, bed1.sub[, -c(1:18), with = FALSE])
    fwrite(combined, outputFile.name, sep = "\t", col.names = TRUE)
    cat("I: Merged BED files successfully. Output written to", outputFile.name, "\n")
    invisible(combined)
}

for(chr in c(2)
#,"X","Y" - yet missing
) {

    cat("I: Merging BED files for chromosome", chr, "...\n")
    merged_data <- mergeTwoBedFiles(chr)
    cat("I: Merged data for chromosome", chr, ":\n")
    cat("   Dimensions: "); print(dim(merged_data))
    cat("   First 19 colnames: "); print(colnames(merged_data)[1:19])
    rm(merged_data)
    gc()
}