#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: compare_tp73_cofactor_summaries.R --input LABEL=FILE [--input ...]",
        "       --output-prefix PATH",
        "",
        "Compare TP73/cofactor summaries produced by",
        "summarize_tp73_patz1_cutandrun_threshold.R. The output table retains",
        "three policies per cofactor: context score >= 0, approximately 50% of",
        "the TP73 >= 0 starts, and the exploratory strongest all-sample gate.",
        "The plot uses only the approximately 50%-retention policy so raw score",
        "scales from different motif models are not compared as equivalent.",
        "",
        "Options:",
        "  --input LABEL=FILE    Label and *_selected.tsv file (repeatable)",
        "  --output-prefix PATH  Output basename for .tsv and .png files",
        "  -h, --help            Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) {
    usage()
}

input_specs <- character()
output_prefix <- NULL
index <- 1L
while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (!option %in% c("--input", "--output-prefix")) {
        writeLines(paste("E: unknown argument:", option), con = stderr())
        usage(2L)
    }
    index <- index + 1L
    if (index > length(arguments)) {
        writeLines(paste("E:", option, "requires a value"), con = stderr())
        usage(2L)
    }
    if (option == "--input") {
        input_specs <- c(input_specs, arguments[[index]])
    } else {
        output_prefix <- arguments[[index]]
    }
    index <- index + 1L
}

if (length(input_specs) < 2L) {
    stop("at least two --input LABEL=FILE arguments are required", call. = FALSE)
}
if (is.null(output_prefix) || !nzchar(output_prefix)) {
    stop("--output-prefix is required", call. = FALSE)
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

parse_input <- function(specification, input_order) {
    separator <- regexpr("=", specification, fixed = TRUE)[[1L]]
    if (separator <= 1L || separator >= nchar(specification)) {
        stop("--input must have the form LABEL=FILE: ", specification,
             call. = FALSE)
    }
    label <- substr(specification, 1L, separator - 1L)
    path <- substr(specification, separator + 1L, nchar(specification))
    if (!file.exists(path)) {
        stop("input does not exist: ", path, call. = FALSE)
    }

    source <- fread(path)
    required <- c(
        "selection", "rule", "context_threshold", "selected",
        "matched_tp73_threshold", "matched_tp73_selected",
        "minimum_incremental_support_ratio", "median_incremental_support_ratio",
        "minimum_incremental_depth_ratio", "median_incremental_depth_ratio"
    )
    missing <- setdiff(required, names(source))
    if (length(missing) > 0L) {
        stop(path, " lacks columns: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }

    baseline_count <- source[rule == "gate", max(as.numeric(selected))]
    strategies <- list(
        score_zero = source[
            rule == "gate" & abs(context_threshold) < 1e-9 &
                grepl("; best ", selection, fixed = TRUE)
        ],
        half_retention = source[
            rule == "gate" &
                grepl("gate retaining about 50%", selection, fixed = TRUE)
        ],
        exploratory = source[
            rule == "gate" &
                grepl("Exploratory strongest all-sample", selection, fixed = TRUE)
        ]
    )
    bad <- names(strategies)[vapply(strategies, nrow, integer(1L)) != 1L]
    if (length(bad) > 0L) {
        stop(path, " does not contain exactly one row for: ",
             paste(bad, collapse = ", "), call. = FALSE)
    }

    result <- rbindlist(lapply(names(strategies), function(strategy) {
        row <- copy(strategies[[strategy]])
        row[, strategy_id := strategy]
        row
    }), use.names = TRUE, fill = TRUE)
    result[, `:=`(
        cofactor = label,
        input_order = input_order,
        baseline_tp73_starts = baseline_count,
        retained_percent = 100 * as.numeric(selected) / baseline_count
    )]
    result
}

inputs <- rbindlist(Map(parse_input, input_specs, seq_along(input_specs)),
                    use.names = TRUE, fill = TRUE)
input_labels <- unique(inputs[, .(input_order, cofactor)])
if (anyDuplicated(input_labels$cofactor)) {
    stop("--input labels must be unique", call. = FALSE)
}

strategy_labels <- c(
    score_zero = "Cofactor score >= 0",
    half_retention = "Approximately 50% retained",
    exploratory = "Exploratory strongest all-sample gate"
)
inputs[, strategy := unname(strategy_labels[strategy_id])]
inputs[, strategy := factor(strategy, levels = unname(strategy_labels))]
setorder(inputs, input_order, strategy)

output_columns <- c(
    "cofactor", "strategy_id", "strategy", "context_threshold", "selected",
    "baseline_tp73_starts", "retained_percent", "matched_tp73_threshold",
    "matched_tp73_selected", "minimum_incremental_support_ratio",
    "median_incremental_support_ratio", "minimum_incremental_depth_ratio",
    "median_incremental_depth_ratio"
)
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
fwrite(inputs[, ..output_columns], paste0(output_prefix, ".tsv"), sep = "\t")

half <- inputs[strategy_id == "half_retention"]
plot_source <- rbindlist(list(
    half[, .(cofactor, input_order, metric = "Strict immersion",
             statistic = "Median",
             improvement_percent = 100 * (median_incremental_support_ratio - 1))],
    half[, .(cofactor, input_order, metric = "Strict immersion",
             statistic = "Minimum across samples",
             improvement_percent = 100 * (minimum_incremental_support_ratio - 1))],
    half[, .(cofactor, input_order, metric = "Maximum depth",
             statistic = "Median",
             improvement_percent = 100 * (median_incremental_depth_ratio - 1))],
    half[, .(cofactor, input_order, metric = "Maximum depth",
             statistic = "Minimum across samples",
             improvement_percent = 100 * (minimum_incremental_depth_ratio - 1))]
))
cofactor_levels <- rev(unique(inputs[order(input_order), cofactor]))
plot_source[, cofactor := factor(cofactor, levels = cofactor_levels)]
plot_source[, metric := factor(metric, levels = c("Strict immersion", "Maximum depth"))]
plot_source[, statistic := factor(
    statistic, levels = c("Median", "Minimum across samples")
)]

comparison_plot <- ggplot(
    plot_source,
    aes(x = improvement_percent, y = cofactor, color = statistic, shape = statistic)
) +
    geom_vline(xintercept = 0, color = "#777777", linewidth = 0.45) +
    geom_line(aes(group = interaction(cofactor, metric)), color = "#B8B8B8",
              linewidth = 0.7) +
    geom_point(size = 3) +
    facet_wrap(~metric, nrow = 1, scales = "free_x") +
    scale_color_manual(values = c("Median" = "#0072B2",
                                  "Minimum across samples" = "#D55E00"), name = NULL) +
    scale_shape_manual(values = c("Median" = 16,
                                  "Minimum across samples" = 17),
                       name = NULL) +
    labs(
        title = "Nearby motif scores add information to TP73 sequence score",
        subtitle = paste0(
            "Each gate retains approximately 50% of TP73 >= 0 starts; ",
            "effects are relative to an exactly count-matched TP73-only rule"
        ),
        x = "Change relative to count-matched TP73-only rule (%)",
        y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA)
    )
ggsave(paste0(output_prefix, ".png"), comparison_plot,
       width = 10.5, height = 5.4, dpi = 180, bg = "white")

message("I: Wrote ", normalizePath(paste0(output_prefix, ".tsv"), mustWork = TRUE))
message("I: Wrote ", normalizePath(paste0(output_prefix, ".png"), mustWork = TRUE))
