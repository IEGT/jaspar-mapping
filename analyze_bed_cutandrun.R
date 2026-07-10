# Calculate the fraction of chromosome bases covered by CUT&RUN intervals.
#
# Statistical question
# --------------------
# For each CUT&RUN BED file, this script estimates the probability that a base
# selected uniformly from a chromosome is covered by at least one CUT&RUN
# interval. This is the background coverage against which the fraction of
# candidate motif positions with CUT&RUN support can be compared.
#
# A base is counted at most once. Summing `end - start` over all reads would
# count bases repeatedly where reads overlap or where duplicate intervals are
# present. The script therefore takes the genomic union of the intervals
# before counting covered bases.
#
# Coordinate and denominator conventions
# --------------------------------------
# CUT&RUN input is interpreted as BED: chromosome, 0-based start, and exclusive
# end. Thus [10, 15) covers five bases. Chromosome lengths are not inferred
# from the reads, because the largest observed read coordinate says only where
# the final read happened to occur. Lengths come from columns 1 and 2 of the
# samtools FASTA index (`.fai`) belonging to the scanned reference genome.
#
# By default the denominator contains chromosomes 1-22, X, and Y. A chromosome
# contributes its complete indexed length even when it has no CUT&RUN reads;
# mitochondrial DNA and unplaced/alternative sequences are excluded unless a
# different chromosome vector is passed explicitly.
#
# Output
# ------
# `calculate_cutandrun_coverage()` returns one row per selected chromosome and
# a final `total` row. The total is a length-weighted fraction:
#
#     sum(union-covered bases) / sum(indexed chromosome lengths)
#
# It is deliberately not the arithmetic mean of chromosome fractions.
# `run_cutandrun_coverage()` adds the source filename and combines all matching
# files into one long table, making chromosome 1 and total baselines directly
# available for calibration.
#
# Running the script
# ------------------
# The default files correspond to the Makefile settings. Build the index with
# `make genome_index`, then run `Rscript analyze_bed_cutandrun.R`. Environment
# variables can override either input location, for example:
#
#     CUTANDRUN_DIR=/data/cutandrun \
#     GENOME_FAI=/data/GRCh38.fasta.fai \
#     Rscript analyze_bed_cutandrun.R

default_chromosomes <- c(as.character(1:22), "X", "Y")

# BED and FASTA files may spell the same chromosome as `chr1` and `1`.
# Internally both are represented as `1`; no other sequence-name rewriting is
# performed.
canonical_chromosome <- function(x) {
    sub("^chr", "", as.character(x), ignore.case = TRUE)
}

# Read the reference lengths for exactly the requested chromosomes.
#
# A standard `.fai` has five columns. Only sequence name and sequence length
# are loaded; byte offsets and line-width fields are ignored. Missing requested
# chromosomes are errors rather than zero-length entries, since silently
# shrinking the denominator would inflate the estimated background coverage.
# Duplicate names after `chr` normalization are also rejected because their
# intended denominator would be ambiguous.
read_fasta_index_lengths <- function(fai_file,
                                     chromosomes = default_chromosomes) {
    if (!file.exists(fai_file)) {
        stop("FASTA index not found: ", fai_file)
    }

    fai <- read.delim(
        fai_file,
        header = FALSE,
        sep = "\t",
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE,
        colClasses = c("character", "numeric", "NULL", "NULL", "NULL")
    )
    if (ncol(fai) != 2L) {
        stop("FASTA index must contain sequence names and lengths: ", fai_file)
    }

    chromosomes <- canonical_chromosome(chromosomes)
    sequence_names <- canonical_chromosome(fai[[1L]])
    selected <- sequence_names %in% chromosomes
    sequence_names <- sequence_names[selected]
    sequence_lengths <- fai[[2L]][selected]

    if (anyDuplicated(sequence_names)) {
        duplicated_names <- unique(sequence_names[duplicated(sequence_names)])
        stop("Duplicate chromosome names after removing the 'chr' prefix: ",
             paste(duplicated_names, collapse = ", "))
    }
    if (any(!is.finite(sequence_lengths) | sequence_lengths <= 0)) {
        stop("FASTA index contains invalid chromosome lengths: ", fai_file)
    }

    lengths <- setNames(sequence_lengths, sequence_names)
    missing_chromosomes <- setdiff(chromosomes, names(lengths))
    if (length(missing_chromosomes) != 0L) {
        stop("FASTA index is missing chromosome(s): ",
             paste(missing_chromosomes, collapse = ", "))
    }

    lengths[chromosomes]
}

# Read only the three BED columns needed for coverage. `fread(select = 1:3)`
# avoids retaining names, scores, strands, or other columns from potentially
# large read files. Sequence names outside the selected denominator are
# ignored, matching the chromosome subset chosen from the FASTA index.
#
# An empty BED file is valid and represents zero covered bases on every
# selected chromosome.
read_cutandrun_intervals <- function(bed_file, chromosomes) {
    if (!file.exists(bed_file)) {
        stop("CUT&RUN BED file not found: ", bed_file)
    }
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("The data.table package is required to read CUT&RUN BED files.")
    }
    if (file.info(bed_file)$size == 0) {
        return(data.frame(
            chromosome = character(),
            start = numeric(),
            end = numeric()
        ))
    }

    intervals <- data.table::fread(
        bed_file,
        header = FALSE,
        select = 1:3,
        colClasses = list(character = 1L, numeric = 2:3),
        data.table = FALSE,
        showProgress = interactive()
    )
    names(intervals) <- c("chromosome", "start", "end")
    intervals$chromosome <- canonical_chromosome(intervals$chromosome)
    intervals[intervals$chromosome %in% chromosomes, , drop = FALSE]
}

# Return the size of the union of a set of half-open intervals.
#
# 1. Sort intervals by start and then end.
# 2. `cummax(end)` records the furthest end reached by any interval so far.
# 3. A new merged block starts only when the next start lies strictly beyond
#    that furthest end. Overlapping, nested, duplicate, and directly adjacent
#    intervals therefore remain in one block. Merging adjacent intervals does
#    not alter the number of covered bases.
# 4. For each block, subtract its first start from its furthest end and sum the
#    resulting non-overlapping lengths.
#
# This vectorized reduction avoids allocating one value per genomic base and
# does not require a Bioconductor interval package.
union_covered_bases <- function(start, end) {
    if (length(start) == 0L) {
        return(0)
    }

    interval_order <- order(start, end)
    start <- start[interval_order]
    end <- end[interval_order]
    cumulative_end <- cummax(end)

    starts_new_block <- c(
        TRUE,
        start[-1L] > cumulative_end[-length(cumulative_end)]
    )
    block_starts <- which(starts_new_block)
    block_ends <- c(block_starts[-1L] - 1L, length(start))

    sum(cumulative_end[block_ends] - start[block_starts])
}

# Calculate chromosome-level and total coverage for one CUT&RUN BED file.
#
# Before reducing intervals, enforce the BED contract: coordinates must be
# finite integers, starts must be non-negative, ends must be greater than
# starts, and no end may exceed the corresponding `.fai` chromosome length.
# Rejecting malformed intervals is preferable to clipping them silently,
# because clipping could conceal a reference-build or chromosome-name mismatch.
#
# Every chromosome named in `chromosome_lengths` receives a row, including a
# zero-coverage row when the BED file contains no selected intervals for it.
calculate_cutandrun_coverage <- function(bed_file, chromosome_lengths) {
    if (is.null(names(chromosome_lengths)) || length(chromosome_lengths) == 0L) {
        stop("chromosome_lengths must be a non-empty named vector")
    }

    chromosomes <- names(chromosome_lengths)
    intervals <- read_cutandrun_intervals(bed_file, chromosomes)

    if (nrow(intervals) != 0L) {
        interval_chromosome_lengths <- unname(
            chromosome_lengths[intervals$chromosome]
        )
        invalid <- !is.finite(intervals$start) |
            !is.finite(intervals$end) |
            intervals$start < 0 |
            intervals$end <= intervals$start |
            intervals$start != floor(intervals$start) |
            intervals$end != floor(intervals$end) |
            intervals$end > interval_chromosome_lengths
        if (any(invalid)) {
            first_invalid <- which(invalid)[1L]
            stop(
                "Invalid or out-of-bounds BED interval at selected row ",
                first_invalid, ": ",
                intervals$chromosome[first_invalid], ":",
                intervals$start[first_invalid], "-",
                intervals$end[first_invalid]
            )
        }
    }

    covered_bases <- vapply(chromosomes, function(chromosome) {
        selected <- intervals$chromosome == chromosome
        union_covered_bases(intervals$start[selected], intervals$end[selected])
    }, numeric(1L))

    result <- data.frame(
        chromosome = chromosomes,
        bp.covered = unname(covered_bases),
        chr.length = unname(chromosome_lengths),
        stringsAsFactors = FALSE
    )
    result$fraction <- result$bp.covered / result$chr.length

    total <- data.frame(
        chromosome = "total",
        bp.covered = sum(result$bp.covered),
        chr.length = sum(result$chr.length),
        fraction = sum(result$bp.covered) / sum(result$chr.length),
        stringsAsFactors = FALSE
    )
    rbind(result, total)
}

# Retain the original scalar helper contract for callers interested only in
# the genome-wide baseline. It selects the already weighted `total` row;
# per-chromosome values remain available from calculate_cutandrun_coverage().
f <- function(fn, chromosome_lengths) {
    coverage <- calculate_cutandrun_coverage(fn, chromosome_lengths)
    coverage$fraction[coverage$chromosome == "total"]
}

# Discover the experiment files, load the reference denominator once, and
# calculate coverage independently for every file. The returned long table has
# columns:
#
#   file, chromosome, bp.covered, chr.length, fraction
#
# The filename pattern preserves the original analysis scope: TP73 CUT&RUN in
# SK-MEL-29, with one `.clipped.clean.bed` input per condition.
run_cutandrun_coverage <- function(
        path_to_cutandrun = "cutandrun_20250602_noDuplicates",
        fai_file = "Homo_sapiens.GRCh38.dna.primary_assembly.fasta.fai",
        chromosomes = default_chromosomes) {
    if (!dir.exists(path_to_cutandrun)) {
        stop("CUT&RUN directory not found: ", path_to_cutandrun)
    }

    cutandrun_files <- list.files(
        path_to_cutandrun,
        pattern = "^tp73_skmel29_2_.*\\.clipped\\.clean\\.bed$",
        full.names = TRUE
    )
    if (length(cutandrun_files) == 0L) {
        stop("No matching CUT&RUN BED files found in: ", path_to_cutandrun)
    }

    chromosome_lengths <- read_fasta_index_lengths(fai_file, chromosomes)
    coverage_tables <- lapply(cutandrun_files, function(cutandrun_file) {
        coverage <- calculate_cutandrun_coverage(
            cutandrun_file,
            chromosome_lengths
        )
        coverage$file <- basename(cutandrun_file)
        coverage[, c("file", "chromosome", "bp.covered", "chr.length",
                     "fraction")]
    })
    result <- do.call(rbind, coverage_tables)
    rownames(result) <- NULL
    result
}

# Execute only when invoked with Rscript. When this file is sourced by a test or
# another analysis, the functions are defined without scanning large files.
if (sys.nframe() == 0L) {
    path_to_cutandrun <- Sys.getenv(
        "CUTANDRUN_DIR",
        unset = "cutandrun_20250602_noDuplicates"
    )
    genome_fai <- Sys.getenv(
        "GENOME_FAI",
        unset = "Homo_sapiens.GRCh38.dna.primary_assembly.fasta.fai"
    )
    print(run_cutandrun_coverage(path_to_cutandrun, genome_fai))
}
