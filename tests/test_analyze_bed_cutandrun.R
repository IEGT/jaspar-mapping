#!/usr/bin/env Rscript

source("analyze_bed_cutandrun.R")

main <- function() {
    test_directory <- tempfile("jaspar-cutandrun-coverage-")
    dir.create(test_directory)
    on.exit(unlink(test_directory, recursive = TRUE), add = TRUE)

    fai_file <- file.path(test_directory, "genome.fasta.fai")
    writeLines(c(
        "1\t100\t0\t80\t81",
        "2\t50\t102\t80\t81",
        "X\t50\t153\t80\t81"
    ), fai_file)

    bed_file <- file.path(
        test_directory,
        "tp73_skmel29_2_TA_R1.clipped.clean.bed"
    )
    writeLines(c(
        "chr1\t0\t10\tread1\t0\t+",
        "chr1\t5\t15\tread2\t0\t+",
        "chr1\t20\t30\tread3\t0\t-",
        "chr2\t40\t50\tread4\t0\t+",
        "chr2\t40\t50\tduplicate\t0\t+",
        "chrUn\t0\t100\tunselected\t0\t+"
    ), bed_file)

    chromosome_lengths <- read_fasta_index_lengths(
        fai_file,
        chromosomes = c("1", "2", "X")
    )
    stopifnot(identical(unname(chromosome_lengths), c(100, 50, 50)))

    coverage <- calculate_cutandrun_coverage(bed_file, chromosome_lengths)
    stopifnot(identical(coverage$chromosome, c("1", "2", "X", "total")))
    stopifnot(identical(coverage$bp.covered, c(25, 10, 0, 35)))
    stopifnot(identical(coverage$chr.length, c(100, 50, 50, 200)))
    stopifnot(isTRUE(all.equal(
        coverage$fraction,
        c(0.25, 0.20, 0, 35 / 200)
    )))
    stopifnot(isTRUE(all.equal(f(bed_file, chromosome_lengths), 35 / 200)))

    reported <- run_cutandrun_coverage(
        test_directory,
        fai_file,
        chromosomes = c("1", "2", "X")
    )
    stopifnot(identical(reported$chromosome, c("1", "2", "X", "total")))
    stopifnot(isTRUE(all.equal(reported$fraction, coverage$fraction)))

    invalid_bed <- file.path(test_directory, "invalid.bed")
    writeLines("1\t90\t101", invalid_bed)
    stopifnot(inherits(
        try(calculate_cutandrun_coverage(
            invalid_bed,
            chromosome_lengths
        ), silent = TRUE),
        "try-error"
    ))

    cat("CUT&RUN coverage tests passed.\n")
}

main()
