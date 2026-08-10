#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: plot_h3k4me3_metaprofile.R --input PROFILE.tsv --output FIGURE.png",
        "       [--table-output PROFILE_RATIO.tsv] [--pseudocount VALUE]",
        "",
        "Plot the GFP-only aggregate H3K4me3/input profile around TP73 anchors.",
        "The input is the small profile table emitted by",
        "build_h3k4me3_anchor_signal.py. TA/DN and cofactor labels are not used,",
        "so this figure can guide window selection without optimizing an outcome.",
        "",
        "Options:",
        "  --input FILE       GFP profile TSV",
        "  --output FILE      PNG figure",
        "  --table-output FILE  Optional derived ratio TSV",
        "  --pseudocount F    Per-base plotting pseudocount (default: 1e-6)",
        "  -h, --help         Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) usage()
values <- list(pseudocount = 1e-6)
index <- 1L
while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (!option %in% c("--input", "--output", "--table-output", "--pseudocount")) {
        writeLines(paste("E: unknown argument:", option), con = stderr())
        usage(2L)
    }
    index <- index + 1L
    if (index > length(arguments)) usage(2L)
    key <- gsub("-", "_", substring(option, 3L), fixed = TRUE)
    values[[key]] <- arguments[[index]]
    index <- index + 1L
}
if (is.null(values$input) || is.null(values$output)) usage(2L)
if (!file.exists(values$input)) stop("profile input not found", call. = FALSE)
values$pseudocount <- suppressWarnings(as.numeric(values$pseudocount))
if (!is.finite(values$pseudocount) || values$pseudocount <= 0) {
    stop("--pseudocount must be positive", call. = FALSE)
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

profile <- fread(values$input, sep = "\t")
required <- c(
    "series_id", "cell_line", "channel", "offset_start_bp", "offset_end_bp",
    "anchors_profiled", "mean_signal_per_bp"
)
if (nrow(profile) == 0L || length(setdiff(required, names(profile))) > 0L ||
    !setequal(unique(profile$channel), c("h3k4me3", "input"))) {
    stop("profile table is empty or incomplete", call. = FALSE)
}
ratio <- dcast(
    profile,
    series_id + cell_line + offset_start_bp + offset_end_bp + anchors_profiled ~
        channel,
    value.var = "mean_signal_per_bp"
)
if (ratio[, any(!is.finite(h3k4me3) | !is.finite(input) |
                 h3k4me3 < 0 | input < 0)]) {
    stop("profile contains invalid signal", call. = FALSE)
}
ratio[, offset_midpoint_bp := (offset_start_bp + offset_end_bp) / 2]
ratio[, log2_h3k4me3_input := log2(
    (h3k4me3 + values$pseudocount) / (input + values$pseudocount)
)]
ratio[, far_flank_median_log2_ratio := median(
    log2_h3k4me3_input[abs(offset_midpoint_bp) >= 1500]
), by = series_id]
ratio[, far_flank_centered_log2_ratio :=
          log2_h3k4me3_input - far_flank_median_log2_ratio]
setorder(ratio, series_id, offset_start_bp)

figure <- ggplot(ratio, aes(offset_midpoint_bp, far_flank_centered_log2_ratio)) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey65") +
    geom_vline(
        xintercept = c(-1000, -500, -150, 0, 150, 500, 1000),
        linewidth = 0.25, linetype = "dotted", colour = "grey75"
    ) +
    geom_line(linewidth = 0.75, colour = "#2166ac") +
    facet_wrap(vars(cell_line, series_id), ncol = 1, scales = "free_y") +
    labs(
        x = "Offset from TP73 motif-alignment centre (bp)",
        y = expression(log[2] * " H3K4me3/input relative to far flanks"),
        title = "GFP-only H3K4me3 profile around TP73 motif anchors",
        subtitle = paste0(
            "Deterministically spaced anchors; plotting pseudocount = ",
            format(values$pseudocount, scientific = TRUE)
        )
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = "grey70")
    )

dir.create(dirname(values$output), recursive = TRUE, showWarnings = FALSE)
ggsave(values$output, figure, width = 9, height = 6, dpi = 160)
if (!is.null(values$table_output)) {
    dir.create(dirname(values$table_output), recursive = TRUE, showWarnings = FALSE)
    fwrite(ratio, values$table_output, sep = "\t")
}
message("I: wrote GFP-only H3K4me3 metaprofile: ", values$output)
