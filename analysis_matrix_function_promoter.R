library(data.table)

# Function to read all BED files matching the pattern into a list of data tables
readPromoterBedFiles <- function() {
    # List all files matching the pattern
    bedFiles <- list.files(path="GeneLists", pattern = "*.promoter.bed", full.names = TRUE)
    
    if (length(bedFiles) == 0) {
        stop("No BED files matching the pattern were found.")
    }
    
    # Read each BED file into a data table
    bedDataTables <- lapply(bedFiles, function(file) {
        fread(file, colClasses = c("Chr"="character", "From"="integer", "To"="integer",
                                   "Gene"="character", "Score"="integer", "Strand"="integer"))
    })
    
    # Assign file names as names of the list
    names(bedDataTables) <- basename(bedFiles)
    
    return(bedDataTables)
}

# Function to check if rows of one BED file overlap with rows of another BED file
checkBedOverlaps <- function(bed1, bed2) {
    # Ensure the data tables have the required columns
    if (!all(c("Chr", "From", "To") %in% colnames(bed1)) || !all(c("Chr", "From", "To") %in% colnames(bed2))) {
        stop("Both BED files must have columns: Chr (chromosome), From (start), End (end).")
    }
    
    # Sort both BED files by chromosome and start position
    setkey(bed1, Chr, From, To)
    setkey(bed2, Chr, From, To)
    
    # Perform a non-equi join to check for overlaps
    overlaps <- foverlaps(
        x = bed1[, .(Chr, From_start = From, To_end = To)], 
        y = bed2[, .(Chr, From_start = From, To_end = To)], 
        by.x = c("Chr", "From_start", "To_end"), 
        by.y = c("Chr", "From_start", "To_end"), 
        type = "any", 
        nomatch = 0L
    )
    
    # Return a logical vector indicating if each row in bed1 overlaps with any row in bed2
    result <- rep(FALSE, nrow(bed1))
    result[overlaps[, xid]] <- TRUE
    return(result)
}

# Example usage
promoterBedTables <- readPromoterBedFiles()

# Example usage:
# bed1 <- fread("path/to/bed1.bed")
# bed2 <- fread("path/to/bed2.bed")
# overlaps <- checkBedOverlaps(bed1, bed2)

