#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: export_bigwig_chrom_bedgraph.R --input FILE.bigWig",
        "       --output FILE.bedGraph [--chrom CHROM]",
        "",
        "Export one chromosome from a BigWig as sorted, non-overlapping bedGraph.",
        "Only finite signal values greater than zero are retained. Coordinates in",
        "the bedGraph output are BED 0-based half-open coordinates.",
        "",
        "Options:",
        "  --input FILE    Source BigWig",
        "  --output FILE   New bedGraph; an existing file is never replaced",
        "  --chrom CHROM   Chromosome with or without chr prefix (default: 1)",
        "  -h, --help      Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) {
    usage()
}

values <- list(chrom = "1")
index <- 1L
while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (!option %in% c("--input", "--output", "--chrom")) {
        writeLines(paste("E: unknown argument:", option), con = stderr())
        usage(2L)
    }
    index <- index + 1L
    if (index > length(arguments)) {
        writeLines(paste("E:", option, "requires a value"), con = stderr())
        usage(2L)
    }
    values[[substring(option, 3L)]] <- arguments[[index]]
    index <- index + 1L
}

if (is.null(values$input) || is.null(values$output)) {
    writeLines("E: --input and --output are required", con = stderr())
    usage(2L)
}
if (!file.exists(values$input)) {
    stop("BigWig input does not exist: ", values$input, call. = FALSE)
}
if (file.exists(values$output)) {
    stop("refusing to replace existing output: ", values$output, call. = FALSE)
}

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
})

normalize_chrom <- function(chrom) sub("^chr", "", chrom, ignore.case = TRUE)

bigwig <- BigWigFile(values$input)
available <- seqlevels(seqinfo(bigwig))
matches <- available[normalize_chrom(available) == normalize_chrom(values$chrom)]
if (length(matches) != 1L) {
    stop(
        "BigWig does not contain exactly one chromosome matching ",
        values$chrom,
        call. = FALSE
    )
}

chrom <- matches[[1L]]
chrom_length <- seqlengths(seqinfo(bigwig))[[chrom]]
if (is.na(chrom_length) || chrom_length <= 0) {
    stop("BigWig chromosome has no usable length: ", chrom, call. = FALSE)
}

ranges <- import(
    bigwig,
    which = GRanges(chrom, IRanges(start = 1L, end = chrom_length))
)
scores <- mcols(ranges)$score
ranges <- ranges[is.finite(scores) & scores > 0]

dir.create(dirname(values$output), recursive = TRUE, showWarnings = FALSE)
export(ranges, values$output, format = "bedGraph")
message(
    "I: exported ", length(ranges), " positive ", chrom,
    " intervals to ", values$output
)
