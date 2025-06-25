library(data.table)

# Function to read all BED files matching the pattern into a list of data tables
readBedFiles <- function(path="GeneLists", pattern="*.bed") {
    # List all files matching the pattern
    bedFiles <- list.files(path=path, pattern = pattern, full.names = TRUE)

    if (length(bedFiles) == 0) {
        stop("No BED files matching the pattern were found.")
    }

    # Read each BED file into a data table
    bedDataTables <- lapply(bedFiles, function(file) {
        dt <- fread(file)  # Read without specifying colClasses
        h <- c("Chr", "From", "To", "Gene", "Score", "Strand")
        if (ncol(dt) > length(h)) {
            if (ncol(dt) < 6) {
                stop("BED file must have at least 6 columns: Chr, From, To, Gene, Score, Strand.")
            }
            if (grepl(pattern="cutandrun", x=file)) {
                combinations <- sort(as.vector(outer(outer(c("pos", "tp73"),c("saos2","skmel29_2"),paste,sep="_"), c("TA", "DN", "GFP"),paste,sep="_")))
                h <- c(h, combinations)  # Add specific combinations for cutandrun files
            }
            if (grepl(pattern="tp73bs", x=file)) {
                combinations <- sort(as.vector(outer(outer(c("pos", "tp73"),c("saos2","skmel29_2"),paste,sep="_"), c("TA", "DN", "GFP"),paste,sep="_")))
                h <- c(h, "num.tfbs")  # Add specific combinations for cutandrun files
            }
            if (ncol(dt) > length(h)) {
                h <- c(h, paste0("V", (length(h) + 1):ncol(dt)))  # Fill with Vn for extra columns
            }
            if (ncol(dt) < length(h)) {
                cat("D: h: ", paste(h,collapse=",",sep=""),"\n",sep="")
                print(dt[1:2,])  # Print first two rows for debugging
                stop("BED file has fewer columns than expected. Expected at least ", length(h), "columns, but found ", ncol(dt), ".")
            }
        }
        setnames(dt, h)  # Rename columns
        dt[, `:=`(Chr = as.character(Chr), From = as.integer(From), To = as.integer(To),
                  Gene = as.character(Gene), Score = as.integer(Score),
                  Strand=as.character(Strand)) ]  # Ensure correct types
       return(dt)
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

    # Add an identifier column to bed2 for tracking overlaps
    bed2[, id := .I]

    # Sort both BED files by chromosome and start position
    setkey(bed1, Chr, From, To)
    setkey(bed2, Chr, From, To)

    # Perform a non-equi join to check for overlaps
    overlaps <- foverlaps(
        x = bed1[, .(Chr, From_start = From, To_end = To)],
        y = bed2[, .(Chr, From_start = From, To_end = To, id)],
        by.x = c("Chr", "From_start", "To_end"),
        by.y = c("Chr", "From_start", "To_end"),
        type = "any",
        nomatch = 0L
    )

    # Return the positions in bed2 where overlaps occur
    return(overlaps$id)
}

# Example usage
promoterBedTables <- readBedFiles(pattern="*.promoter.*bed")
utrBedTables <- readBedFiles(pattern="*.utr.*bed")


# Example usage:
# bed1 <- fread("path/to/bed1.bed")
# bed2 <- fread("path/to/bed2.bed")
# overlaps <- checkBedOverlaps(bed1, bed2)

