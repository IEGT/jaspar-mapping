# Scripts to determine the coverage of the chromosome with
# CUT&RUN data.

path_to_cutandrun <- "cutandrun_20250602_noDuplicates"
cutandrun_files <- list.files(path_to_cutandrun, pattern = "tp73_skmel29_2_.*clipped.clean.bed", full.names = TRUE)

# Load the necessary libraries

f <- function(fn) {

    if (!file.exists(fn)) {
        stop(paste("File not found:", x))
    }

    d <- read.table(fn, header = FALSE)

    # Calculate the coverage for each chromosome
    coverage_per_chromosome <- aggregate(d$V3 - d$V2, by = list(Chromosome = d$V1), FUN = sum)
    max_per_chromosome <- aggregate(d$V3, by = list(Chromosome = d$V1), FUN = max)

    sapply(c(1:22, "X", "Y"), function(chr) {
    if (chr %in% coverage_per_chromosome$Chromosome) {
        coverage <- coverage_per_chromosome[coverage_per_chromosome$Chromosome == chr, 2]
        max_value <- max_per_chromosome[max_per_chromosome$Chromosome == chr, 2]
        return(c(bp.covered=coverage, chr.length=max_value, fraction = coverage / max_value))
    } else {
        return(c(0, 0))
    }
    }) -> coverage_results

    total_fraction_covered <- sum(coverage_results["bp.covered",])/sum(coverage_results["chr.length",])
}

sapply(cutandrun_files, f)
