library(ggplot2)
library(data.table)

# Function to read .bed files and return a data table
readBedFile <- function(file) {
    dt <- fread(file, col.names = c("Chr", "Start", "End", "Name", "Score", "Strand"))
    return(dt)
}

# Function to find binding sites within 100 bp of TP73
findNearbySites <- function(tf_sites, tp73_sites, buffer = 100) {
    tf_sites[, IsNearbyTP73 := FALSE]
    for (i in 1:nrow(tp73_sites)) {
        tp73_start <- tp73_sites[i, Start]
        tp73_end <- tp73_sites[i, End]
        tf_sites[Start <= tp73_end + buffer & End >= tp73_start - buffer, IsNearbyTP73 := TRUE]
    }
    return(tf_sites)
}


tfbs.directory <- "output_RelativeRisk_20250217/"
find.data.for.TF.on.chromosome <- function(name,chr) {
    tf.file.name.pattern <- paste0("^",name,"_.*bidirect_",chr,".bed.gz")
    tf.file.name.path <- paste0(tfbs.directory,chr)
    tf.file.name <- list.files(path=tf.file.name.path,
                            pattern=tf.file.name.pattern, full.names = TRUE)
    if (length(tf.file.name) == 0) {
        stop("E: No TFBS file found for ", name, " on chromosome ", chr,
            ", tried '",tf.file.name.path,"' with pattern '",tf.file.name.pattern,"'\n",sep = "")
        return(NULL)
    }
    if (length(tf.file.name) > 1) {
        stop("E: Multiple files found for ", name, " on chromosome ", chr)
        return(NULL)
    }
    cat("I: Found TFBS file for ", name, " on chromosome ", chr, ": '", tf.file.name, "'\n", sep = "")
    if (file.exists(tf.file.name)) {
        tfbs <- readBedFile(tf.file.name)
        tfbs[, Chr := chr]
        return(tfbs)
    } else {
        stop("E: No TFBS file found for ", name, " on chromosome ", chr, ", tried '",tf.file.name,"'\n", sep = "")
        return(NULL)
    }
}

find.data.for.TF.on.chromosome.test <- find.data.for.TF.on.chromosome("TP73", "22")

find.data.for.TF.on.all.chromosomes <- function(name) {
    tfbs.list <- lapply(c(1:22,"X", "Y"), function(chr) {
        tfbs <- find.data.for.TF.on.chromosome(name, chr)
        return(tfbs)
    })
    #tfbs.all <- rbindlist(tfbs.list, use.names = TRUE, fill = TRUE)
    #return(tfbs.all)
}

find.data.for.TF.on.all.chromosomes.test <- find.data.for.TF.on.all.chromosomes("TP73")

transcription.factor.selection <- c("TP73", "HEY1", "RORA_MA0071.1", "E2F1", "SP1")

# Example input: List of genes and transcription start sites
genes <- data.frame(
    Gene = c("NFKBIA", "IL10RA", "IL4R", "LIMA1", "CDKN1A", "CDKN2A", "BBC3", "BAX", "FOXJ1", "RFX2", "MDM2", "WWOX", "FHIT"),
    Chr = c("14",     "11",     "16",    "12"),
    TSS = c(35405249 , 117985941, 27313426, 50283520),
    strand = c("-", "+", "+","-")
)

tfbs <- lapply(unique(c("TP73",transcription.factor.selection)),find.data.for.TF.on.all.chromosomes)

# Find nearby sites for each transcription factor
tf_binding_sites <- lapply(tf_binding_sites, function(tf_sites) {
    findNearbySites(tf_sites, tp73_sites)
})

# Combine all transcription factor binding sites into one data table
all_tf_sites <- rbindlist(tf_binding_sites, idcol = "TF")

# Plot the results
ggplot() +
    # Plot genes as horizontal lines
    geom_segment(data = genes, aes(x = TSS - 500, xend = TSS + 500, y = Gene, yend = Gene), color = "black") +
    # Plot transcription factor binding sites
    geom_point(data = all_tf_sites, aes(x = (Start + End) / 2, y = TF, color = IsNearbyTP73), size = 2) +
    # Highlight TP73 binding regions
    geom_rect(data = tp73_sites, aes(xmin = Start - 1, xmax = End + 1, ymin = -Inf, ymax = Inf), alpha = 0.2, fill = "red") +
    labs(title = "Transcription Factor Binding Sites", x = "Genomic Position", y = "Gene/TF") +
    theme_minimal() +
    theme(legend.position = "bottom")