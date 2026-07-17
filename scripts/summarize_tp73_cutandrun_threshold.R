#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: summarize_tp73_cutandrun_threshold.R --run-root DIR",
        "       [--output-prefix PATH] [--minimum-effect FRACTION]",
        "       [--risk-threshold SCORE]",
        "",
        "Compare chromosome-wide anti-p73 and matched IgG calibration curves for",
        "risk_p0, risk_p1, and log_odds_p1. The script writes a compact aggregate",
        "table, a threshold recommendation table, and three PNG figures.",
        "",
        "The data-derived lower bound is the lowest 0.2-score threshold at which",
        "all TA/DN samples show the requested relative advantage over IgG for both",
        "strict-immersion probability and baseline-normalized maximum depth.",
        "",
        "Options:",
        "  --run-root DIR          Directory containing risk_p0, risk_p1,",
        "                          log_odds_p1 and matching igg_* directories",
        "  --output-prefix PATH    Output basename (default: DIR/tp73_threshold)",
        "  --minimum-effect N      Required anti-p73/IgG advantage (default: 0.05)",
        "  --risk-threshold SCORE  Rounded operational risk threshold (default: 0)",
        "  -h, --help              Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) {
    usage()
}

values <- list(minimum_effect = 0.05, risk_threshold = 0)
index <- 1L
while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (!option %in% c(
        "--run-root", "--output-prefix", "--minimum-effect", "--risk-threshold"
    )) {
        writeLines(paste("E: unknown argument:", option), con = stderr())
        usage(2L)
    }
    index <- index + 1L
    if (index > length(arguments)) {
        writeLines(paste("E:", option, "requires a value"), con = stderr())
        usage(2L)
    }
    name <- gsub("-", "_", substring(option, 3L), fixed = TRUE)
    values[[name]] <- arguments[[index]]
    index <- index + 1L
}

if (is.null(values$run_root)) {
    writeLines("E: --run-root is required", con = stderr())
    usage(2L)
}
values$minimum_effect <- as.numeric(values$minimum_effect)
values$risk_threshold <- as.numeric(values$risk_threshold)
if (!is.finite(values$minimum_effect) || values$minimum_effect < 0) {
    stop("--minimum-effect must be a finite non-negative number", call. = FALSE)
}
if (!is.finite(values$risk_threshold)) {
    stop("--risk-threshold must be finite", call. = FALSE)
}

run_root <- normalizePath(values$run_root, mustWork = TRUE)
if (is.null(values$output_prefix)) {
    values$output_prefix <- file.path(run_root, "tp73_threshold")
}
output_prefix <- values$output_prefix
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

modes <- c("risk_p0", "risk_p1", "log_odds_p1")
mode_labels <- c(
    risk_p0 = "Log2 relative risk, pseudocount 0",
    risk_p1 = "Log2 relative risk, pseudocount 1",
    log_odds_p1 = "Log odds, pseudocount 1"
)

read_curves <- function(prefix) {
    result <- rbindlist(lapply(modes, function(mode) {
        path <- file.path(run_root, paste0(prefix, mode), "threshold_curve.tsv")
        if (!file.exists(path)) {
            stop("missing threshold curve: ", path, call. = FALSE)
        }
        table <- fread(path)
        table[, mode := mode]
        table
    }), use.names = TRUE)
    result
}

anti <- read_curves("")
igg <- read_curves("igg_")

baseline_depth <- function(table, output_name) {
    baseline <- table[
        order(-selected_motifs),
        .SD[1L],
        by = .(mode, sample_id)
    ][, .(mode, sample_id, baseline_depth = mean_effective_depth)]
    setnames(baseline, "baseline_depth", output_name)
    baseline
}

anti_baseline <- baseline_depth(anti, "anti_baseline_depth")
igg_baseline <- baseline_depth(igg, "igg_baseline_depth")

columns <- c(
    "mode", "sample_id", "threshold", "selected_motifs",
    "precision_enrichment", "mean_effective_depth"
)
pairs <- merge(
    anti[, ..columns],
    igg[, ..columns],
    by = c("mode", "sample_id", "threshold", "selected_motifs"),
    suffixes = c("_anti", "_igg")
)
pairs <- merge(pairs, anti_baseline, by = c("mode", "sample_id"))
pairs <- merge(pairs, igg_baseline, by = c("mode", "sample_id"))
pairs[, condition := sub("^.*_", "", sample_id)]
pairs <- pairs[condition %in% c("TA", "DN") & is.finite(threshold)]
pairs[, support_ratio := precision_enrichment_anti / precision_enrichment_igg]
pairs[, anti_depth_enrichment :=
    mean_effective_depth_anti / anti_baseline_depth]
pairs[, igg_depth_enrichment :=
    mean_effective_depth_igg / igg_baseline_depth]
pairs[, depth_ratio := anti_depth_enrichment / igg_depth_enrichment]

aggregate <- pairs[, .(
    selected_motifs = selected_motifs[[1L]],
    minimum_support_ratio = min(support_ratio),
    median_support_ratio = median(support_ratio),
    maximum_support_ratio = max(support_ratio),
    minimum_depth_ratio = min(depth_ratio),
    median_depth_ratio = median(depth_ratio),
    maximum_depth_ratio = max(depth_ratio)
), by = .(mode, threshold)]
setorder(aggregate, mode, threshold)

required_ratio <- 1 + values$minimum_effect
lower_bounds <- aggregate[
    selected_motifs >= 1000 &
        minimum_support_ratio >= required_ratio &
        minimum_depth_ratio >= required_ratio,
    .SD[1L],
    by = mode
]
lower_bounds[, recommendation_kind := "data_derived_lower_bound"]

operational_risk_p1 <- aggregate[
    mode == "risk_p1",
    .SD[which.min(abs(threshold - values$risk_threshold))]
]
target_count <- operational_risk_p1$selected_motifs[[1L]]
operational <- rbindlist(lapply(modes, function(current_mode) {
    candidates <- aggregate[mode == current_mode]
    if (current_mode %in% c("risk_p0", "risk_p1")) {
        candidates[which.min(abs(threshold - values$risk_threshold))]
    } else {
        candidates[which.min(abs(selected_motifs - target_count))]
    }
}))
operational[, recommendation_kind := "rounded_operational_threshold"]

recommendations <- rbindlist(list(lower_bounds, operational), use.names = TRUE)
recommendations[, minimum_effect_required := values$minimum_effect]
recommendations[, mode_label := unname(mode_labels[mode])]
setcolorder(recommendations, c(
    "recommendation_kind", "mode", "mode_label", "threshold",
    "selected_motifs", "minimum_effect_required",
    "minimum_support_ratio", "median_support_ratio",
    "minimum_depth_ratio", "median_depth_ratio"
))

fwrite(aggregate, paste0(output_prefix, "_aggregate.tsv"), sep = "\t")
fwrite(recommendations, paste0(output_prefix, "_recommendations.tsv"), sep = "\t")

plot_data <- aggregate[
    selected_motifs >= 1000 & (
        (mode %in% c("risk_p0", "risk_p1") & threshold >= -20 & threshold <= 14) |
        (mode == "log_odds_p1" & threshold >= 0 & threshold <= 55)
    )
]
plot_data[, mode_label := factor(
    unname(mode_labels[mode]),
    levels = unname(mode_labels[modes])
)]
operational_plot <- copy(operational)
operational_plot[, mode_label := factor(
    unname(mode_labels[mode]),
    levels = unname(mode_labels[modes])
)]

support_plot <- ggplot(
    plot_data,
    aes(x = threshold)
) +
    geom_hline(yintercept = 1, color = "#777777", linewidth = 0.4) +
    geom_hline(
        yintercept = required_ratio, color = "#777777",
        linewidth = 0.4, linetype = "dashed"
    ) +
    geom_line(
        aes(y = minimum_support_ratio, linetype = "Least-enriched sample"),
        color = "#56B4E9", linewidth = 0.7
    ) +
    geom_line(
        aes(y = median_support_ratio, linetype = "Median"),
        color = "#0072B2", linewidth = 0.9
    ) +
    geom_point(
        data = operational_plot,
        aes(x = threshold, y = median_support_ratio),
        inherit.aes = FALSE, color = "#D55E00", size = 2
    ) +
    facet_wrap(~mode_label, scales = "free_x", nrow = 1) +
    scale_linetype_manual(
        values = c("Median" = "solid", "Least-enriched sample" = "dashed"),
        name = NULL
    ) +
    coord_cartesian(ylim = c(0.98, 1.45)) +
    labs(
        title = "TP73 score and strict CUT&RUN immersion",
        subtitle = "Median and least-enriched of six TA/DN anti-p73 vs matched IgG comparisons",
        x = "Score threshold",
        y = "Anti-p73 / IgG normalized immersion enrichment"
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
    )

depth_plot <- ggplot(
    plot_data,
    aes(x = threshold)
) +
    geom_hline(yintercept = 1, color = "#777777", linewidth = 0.4) +
    geom_hline(
        yintercept = required_ratio, color = "#777777",
        linewidth = 0.4, linetype = "dashed"
    ) +
    geom_line(
        aes(y = minimum_depth_ratio, linetype = "Least-enriched sample"),
        color = "#6CC4A1", linewidth = 0.7
    ) +
    geom_line(
        aes(y = median_depth_ratio, linetype = "Median"),
        color = "#007F5F", linewidth = 0.9
    ) +
    geom_point(
        data = operational_plot,
        aes(x = threshold, y = median_depth_ratio),
        inherit.aes = FALSE, color = "#D55E00", size = 2
    ) +
    facet_wrap(~mode_label, scales = "free_x", nrow = 1) +
    scale_linetype_manual(
        values = c("Median" = "solid", "Least-enriched sample" = "dashed"),
        name = NULL
    ) +
    coord_cartesian(ylim = c(0.98, 1.7)) +
    labs(
        title = "TP73 score and maximum CUT&RUN depth",
        subtitle = "Depth is normalized to each track's unthresholded mean before anti-p73/IgG comparison",
        x = "Score threshold",
        y = "Anti-p73 / IgG normalized depth enrichment"
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
    )

tradeoff_plot <- ggplot(
    aggregate[selected_motifs >= 1000],
    aes(
        x = selected_motifs, y = median_support_ratio,
        color = factor(mode, levels = modes)
    )
) +
    geom_hline(yintercept = 1, color = "#777777", linewidth = 0.4) +
    geom_hline(
        yintercept = required_ratio, color = "#777777",
        linewidth = 0.4, linetype = "dashed"
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(
        data = operational,
        aes(x = selected_motifs, y = median_support_ratio, color = mode),
        inherit.aes = FALSE, size = 2.5
    ) +
    scale_x_log10(labels = scales::label_number(big.mark = ",")) +
    scale_color_manual(
        values = c("#CC79A7", "#0072B2", "#D55E00"),
        labels = unname(mode_labels[modes]),
        name = NULL
    ) +
    labs(
        title = "Specificity versus chromosome-1 storage",
        subtitle = "Orange points mark the rounded risk threshold and equal-count alternatives",
        x = "Retained chromosome-1 motif-model starts (log scale)",
        y = "Median anti-p73 / IgG immersion enrichment"
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )

ggsave(
    paste0(output_prefix, "_support.png"),
    support_plot, width = 14, height = 5.2, dpi = 180
)
ggsave(
    paste0(output_prefix, "_depth.png"),
    depth_plot, width = 14, height = 5.2, dpi = 180
)
ggsave(
    paste0(output_prefix, "_storage_tradeoff.png"),
    tradeoff_plot, width = 8.5, height = 5.5, dpi = 180
)

print(recommendations)
