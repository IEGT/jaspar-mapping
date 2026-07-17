#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: summarize_tp73_patz1_cutandrun_threshold.R --run-root DIR",
        "       [--output-prefix PATH] [--context-label NAME]",
        "       [--anchor-threshold SCORE]",
        "       [--rescue-lower-threshold SCORE] [--minimum-effect FRACTION]",
        "",
        "Summarize a streaming TP73/cofactor joint-score calibration. For every TP73",
        "alignment span, the input must contain the highest orientation-collapsed",
        "cofactor score within the configured context radius. Two rules are compared:",
        "",
        "  gate:   TP73 >= anchor threshold AND best cofactor >= context threshold",
        "  rescue: TP73 >= anchor threshold OR",
        "          rescue lower <= TP73 < anchor threshold AND best cofactor >= threshold",
        "",
        "Each rule is also compared with a TP73-only score boundary at exactly the",
        "same expected storage count by interpolating within its tied 0.2-score bin.",
        "Selected rows include the strongest uniformly beneficial and uniformly",
        "detrimental gates when all six TA/DN comparisons agree on the direction.",
        "TA and DN anti-p73 tracks are evaluated against",
        "their matched IgG tracks; GFP is retained in raw inputs but not thresholding.",
        "",
        "Expected directories below --run-root:",
        "  joint_risk_p1/      anti-p73 streaming output",
        "  igg_joint_risk_p1/  matched IgG streaming output",
        "",
        "Options:",
        "  --run-root DIR                 Joint calibration root (required)",
        "  --output-prefix PATH           Output basename (default: DIR/tp73_LABEL)",
        "  --context-label NAME           Cofactor label used in tables/plots (default: PATZ1)",
        "  --anchor-threshold SCORE       TP73 operational cutoff (default: 0)",
        "  --rescue-lower-threshold SCORE Lower TP73 rescue bound (default: -2)",
        "  --minimum-effect FRACTION      Reference effect size (default: 0.05)",
        "  -h, --help                     Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) {
    usage()
}

values <- list(
    context_label = "PATZ1",
    anchor_threshold = 0,
    rescue_lower_threshold = -2,
    minimum_effect = 0.05
)
index <- 1L
while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (!option %in% c(
        "--run-root", "--output-prefix", "--context-label", "--anchor-threshold",
        "--rescue-lower-threshold", "--minimum-effect"
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
if (!nzchar(values$context_label) || grepl("[\t\r\n]", values$context_label)) {
    stop("--context-label must be non-empty and may not contain tabs or newlines",
         call. = FALSE)
}
for (name in c("anchor_threshold", "rescue_lower_threshold", "minimum_effect")) {
    values[[name]] <- as.numeric(values[[name]])
    if (!is.finite(values[[name]])) {
        stop("--", gsub("_", "-", name), " must be finite", call. = FALSE)
    }
}
if (values$rescue_lower_threshold >= values$anchor_threshold) {
    stop("--rescue-lower-threshold must be below --anchor-threshold", call. = FALSE)
}
if (values$minimum_effect < 0) {
    stop("--minimum-effect must be non-negative", call. = FALSE)
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

run_root <- normalizePath(values$run_root, mustWork = TRUE)
if (is.null(values$output_prefix)) {
    context_slug <- tolower(gsub("[^A-Za-z0-9]+", "_", values$context_label))
    context_slug <- gsub("^_+|_+$", "", context_slug)
    values$output_prefix <- file.path(run_root, paste0("tp73_", context_slug))
}
output_prefix <- values$output_prefix
context_label <- values$context_label
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

anti_directory <- file.path(run_root, "joint_risk_p1")
igg_directory <- file.path(run_root, "igg_joint_risk_p1")
for (directory in c(anti_directory, igg_directory)) {
    for (name in c("joint_score_histogram.tsv", "threshold_curve.tsv")) {
        path <- file.path(directory, name)
        if (!file.exists(path)) {
            stop("missing calibration input: ", path, call. = FALSE)
        }
    }
}

reverse_cumsum_2d <- function(input) {
    result <- matrix(0, nrow = nrow(input), ncol = ncol(input))
    for (row in rev(seq_len(nrow(input)))) {
        for (column in rev(seq_len(ncol(input)))) {
            value <- input[row, column]
            if (row < nrow(input)) value <- value + result[row + 1L, column]
            if (column < ncol(input)) value <- value + result[row, column + 1L]
            if (row < nrow(input) && column < ncol(input)) {
                value <- value - result[row + 1L, column + 1L]
            }
            result[row, column] <- value
        }
    }
    result
}

joint_curve <- function(path) {
    source <- fread(path)
    required <- c(
        "sample_id", "anchor_bin_index", "anchor_threshold",
        "context_bin_index", "context_threshold", "n_total", "n_supported",
        "effective_depth_sum", "n_same_orientation", "n_opposite_orientation",
        "n_supported_same_orientation", "n_supported_opposite_orientation",
        "mean_signed_center_distance_bp", "mean_absolute_center_distance_bp"
    )
    missing <- setdiff(required, names(source))
    if (length(missing) > 0L) {
        stop("joint histogram lacks columns: ", paste(missing, collapse = ", "),
             call. = FALSE)
    }
    if (source[, anyDuplicated(paste(sample_id, anchor_bin_index,
                                     context_bin_index, sep = ":"))] != 0L) {
        stop("joint histogram contains duplicate sample/bin rows", call. = FALSE)
    }

    anchor_bins <- sort(unique(source$anchor_bin_index))
    context_bins <- sort(unique(source$context_bin_index))
    anchor_thresholds <- source[, .(anchor_threshold = anchor_threshold[[1L]]),
                                by = anchor_bin_index]
    context_thresholds <- source[, .(context_threshold = context_threshold[[1L]]),
                                 by = context_bin_index]
    setkey(anchor_thresholds, anchor_bin_index)
    setkey(context_thresholds, context_bin_index)
    dimensions <- c(length(anchor_bins), length(context_bins))
    flat_indexes <- arrayInd(seq_len(prod(dimensions)), .dim = dimensions)

    make_matrix <- function(table, column) {
        result <- matrix(0, nrow = dimensions[[1L]], ncol = dimensions[[2L]])
        positions <- cbind(
            match(table$anchor_bin_index, anchor_bins),
            match(table$context_bin_index, context_bins)
        )
        result[positions] <- as.numeric(table[[column]])
        reverse_cumsum_2d(result)
    }

    rbindlist(lapply(unique(source$sample_id), function(sample) {
        table <- source[sample_id == sample]
        selected <- make_matrix(table, "n_total")
        supported <- make_matrix(table, "n_supported")
        depth <- make_matrix(table, "effective_depth_sum")
        same <- make_matrix(table, "n_same_orientation")
        opposite <- make_matrix(table, "n_opposite_orientation")
        supported_same <- make_matrix(table, "n_supported_same_orientation")
        supported_opposite <- make_matrix(table, "n_supported_opposite_orientation")
        table[, signed_distance_sum := mean_signed_center_distance_bp * n_total]
        table[, absolute_distance_sum := mean_absolute_center_distance_bp * n_total]
        signed_distance <- make_matrix(table, "signed_distance_sum")
        absolute_distance <- make_matrix(table, "absolute_distance_sum")

        result <- data.table(
            sample_id = sample,
            anchor_bin_index = anchor_bins[flat_indexes[, 1L]],
            context_bin_index = context_bins[flat_indexes[, 2L]],
            selected = as.vector(selected),
            supported = as.vector(supported),
            depth = as.vector(depth),
            same_orientation = as.vector(same),
            opposite_orientation = as.vector(opposite),
            supported_same_orientation = as.vector(supported_same),
            supported_opposite_orientation = as.vector(supported_opposite),
            signed_distance_sum = as.vector(signed_distance),
            absolute_distance_sum = as.vector(absolute_distance)
        )
        result <- anchor_thresholds[result, on = "anchor_bin_index"]
        result <- context_thresholds[result, on = "context_bin_index"]
        result[selected > 0]
    }))
}

read_normalization <- function(directory, prefix) {
    curve <- fread(file.path(directory, "threshold_curve.tsv"))
    result <- curve[order(-selected_motifs), .SD[1L], by = sample_id,
                    .SDcols = c("support_prevalence", "mean_effective_depth")]
    setnames(result,
             c("support_prevalence", "mean_effective_depth"),
             paste0(prefix, c("_prevalence", "_baseline_depth")))
    result
}

normalization <- merge(
    read_normalization(anti_directory, "anti"),
    read_normalization(igg_directory, "igg"),
    by = "sample_id"
)

anti_curve <- joint_curve(file.path(anti_directory, "joint_score_histogram.tsv"))
igg_curve <- joint_curve(file.path(igg_directory, "joint_score_histogram.tsv"))
joint_keys <- c(
    "sample_id", "anchor_bin_index", "anchor_threshold",
    "context_bin_index", "context_threshold", "selected"
)
joint_pairs <- merge(
    anti_curve, igg_curve, by = joint_keys, suffixes = c("_anti", "_igg")
)
joint_pairs <- merge(joint_pairs, normalization, by = "sample_id")
joint_pairs[, condition := sub("^.*_", "", sample_id)]

calculate_ratios <- function(table) {
    table[, support_ratio :=
        ((supported_anti / selected) / anti_prevalence) /
        ((supported_igg / selected) / igg_prevalence)]
    table[, depth_ratio :=
        ((depth_anti / selected) / anti_baseline_depth) /
        ((depth_igg / selected) / igg_baseline_depth)]
    table
}
joint_pairs <- calculate_ratios(joint_pairs)

one_dimensional_pairs <- function() {
    columns <- c(
        "sample_id", "threshold", "selected_motifs",
        "supported_selected_motifs", "precision_enrichment",
        "mean_effective_depth"
    )
    anti <- fread(file.path(anti_directory, "threshold_curve.tsv"))[, ..columns]
    igg <- fread(file.path(igg_directory, "threshold_curve.tsv"))[, ..columns]
    result <- merge(
        anti, igg,
        by = c("sample_id", "threshold", "selected_motifs"),
        suffixes = c("_anti", "_igg")
    )
    result <- merge(result, normalization, by = "sample_id")
    result[, condition := sub("^.*_", "", sample_id)]
    result[, support_ratio := precision_enrichment_anti / precision_enrichment_igg]
    result[, depth_ratio :=
        (mean_effective_depth_anti / anti_baseline_depth) /
        (mean_effective_depth_igg / igg_baseline_depth)]
    result[condition %in% c("TA", "DN") & is.finite(threshold)]
}
one_dimensional <- one_dimensional_pairs()
one_counts <- unique(one_dimensional[, .(threshold, selected_motifs)])
if (one_counts[, anyDuplicated(threshold)] != 0L) {
    stop("TP73-only candidate counts differ among samples", call. = FALSE)
}

matched_tp73_threshold <- function(selected) {
    index <- which.min(abs(as.numeric(one_counts$selected_motifs) -
                           as.numeric(selected)))
    list(
        threshold = one_counts$threshold[[index]],
        selected = as.numeric(one_counts$selected_motifs[[index]])
    )
}

interpolated_tp73_comparator <- function(sample, selected) {
    candidates <- one_dimensional[sample_id == sample][order(selected_motifs)]
    counts <- as.numeric(candidates$selected_motifs)
    target <- as.numeric(selected)
    lower_index <- max(which(counts <= target))
    upper_index <- min(which(counts >= target))
    lower <- candidates[lower_index]
    upper <- candidates[upper_index]
    if (counts[[upper_index]] == counts[[lower_index]]) {
        fraction <- 0
    } else {
        fraction <- (target - counts[[lower_index]]) /
            (counts[[upper_index]] - counts[[lower_index]])
    }
    interpolate <- function(lower_value, upper_value) {
        lower_value + fraction * (upper_value - lower_value)
    }
    anti_supported <- interpolate(
        lower$supported_selected_motifs_anti,
        upper$supported_selected_motifs_anti
    )
    igg_supported <- interpolate(
        lower$supported_selected_motifs_igg,
        upper$supported_selected_motifs_igg
    )
    anti_depth <- interpolate(
        lower$mean_effective_depth_anti * counts[[lower_index]],
        upper$mean_effective_depth_anti * counts[[upper_index]]
    )
    igg_depth <- interpolate(
        lower$mean_effective_depth_igg * counts[[lower_index]],
        upper$mean_effective_depth_igg * counts[[upper_index]]
    )
    data.table(
        sample_id = sample,
        selected = target,
        tp73_comparator_lower_threshold = lower$threshold,
        tp73_comparator_upper_threshold = upper$threshold,
        tp73_comparator_lower_count = counts[[lower_index]],
        tp73_comparator_upper_count = counts[[upper_index]],
        tp73_only_support_ratio =
            ((anti_supported / target) / lower$anti_prevalence) /
            ((igg_supported / target) / lower$igg_prevalence),
        tp73_only_depth_ratio =
            ((anti_depth / target) / lower$anti_baseline_depth) /
            ((igg_depth / target) / lower$igg_baseline_depth)
    )
}

attach_matched_comparison <- function(table) {
    mapping <- unique(table[, .(context_threshold, selected)])
    nearest <- lapply(mapping$selected, matched_tp73_threshold)
    mapping[, matched_tp73_threshold := vapply(nearest, `[[`, numeric(1L),
                                                "threshold")]
    mapping[, matched_tp73_selected := vapply(nearest, `[[`, numeric(1L),
                                               "selected")]
    result <- merge(table, mapping, by = c("context_threshold", "selected"))
    comparator <- rbindlist(lapply(unique(result$selected), function(selected) {
        rbindlist(lapply(unique(result$sample_id), function(sample) {
            interpolated_tp73_comparator(sample, selected)
        }))
    }))
    result <- merge(result, comparator, by = c("sample_id", "selected"))
    result[, incremental_support_ratio := support_ratio / tp73_only_support_ratio]
    result[, incremental_depth_ratio := depth_ratio / tp73_only_depth_ratio]
    result
}

tolerance <- 1e-9
minimum_context_threshold <- min(joint_pairs$context_threshold)
gate <- joint_pairs[
    abs(anchor_threshold - values$anchor_threshold) < tolerance &
        condition %in% c("TA", "DN")
]
gate[, `:=`(
    same_anti = same_orientation_anti,
    opposite_anti = opposite_orientation_anti,
    supported_same_anti = supported_same_orientation_anti,
    supported_opposite_anti = supported_opposite_orientation_anti
)]

baseline <- joint_pairs[
    abs(anchor_threshold - values$anchor_threshold) < tolerance &
        abs(context_threshold - minimum_context_threshold) < tolerance &
        condition %in% c("TA", "DN"),
    .(
        sample_id, condition,
        baseline_selected = selected,
        baseline_supported_anti = supported_anti,
        baseline_supported_igg = supported_igg,
        baseline_depth_anti = depth_anti,
        baseline_depth_igg = depth_igg,
        baseline_same_anti = same_orientation_anti,
        baseline_opposite_anti = opposite_orientation_anti,
        baseline_supported_same_anti = supported_same_orientation_anti,
        baseline_supported_opposite_anti = supported_opposite_orientation_anti,
        baseline_signed_distance_sum_anti = signed_distance_sum_anti,
        baseline_absolute_distance_sum_anti = absolute_distance_sum_anti
    )
]

lower <- joint_pairs[
    abs(anchor_threshold - values$rescue_lower_threshold) < tolerance &
        condition %in% c("TA", "DN"),
    .(
        sample_id, condition, context_threshold,
        lower_selected = selected,
        lower_supported_anti = supported_anti,
        lower_supported_igg = supported_igg,
        lower_depth_anti = depth_anti,
        lower_depth_igg = depth_igg,
        lower_same_anti = same_orientation_anti,
        lower_opposite_anti = opposite_orientation_anti,
        lower_supported_same_anti = supported_same_orientation_anti,
        lower_supported_opposite_anti = supported_opposite_orientation_anti,
        lower_signed_distance_sum_anti = signed_distance_sum_anti,
        lower_absolute_distance_sum_anti = absolute_distance_sum_anti
    )
]
upper <- gate[, .(
    sample_id, condition, context_threshold,
    upper_selected = selected,
    upper_supported_anti = supported_anti,
    upper_supported_igg = supported_igg,
    upper_depth_anti = depth_anti,
    upper_depth_igg = depth_igg,
    upper_same_anti = same_orientation_anti,
    upper_opposite_anti = opposite_orientation_anti,
    upper_supported_same_anti = supported_same_orientation_anti,
    upper_supported_opposite_anti = supported_opposite_orientation_anti,
    upper_signed_distance_sum_anti = signed_distance_sum_anti,
    upper_absolute_distance_sum_anti = absolute_distance_sum_anti
)]

rescue <- merge(lower, upper, by = c("sample_id", "condition", "context_threshold"))
rescue <- merge(rescue, baseline, by = c("sample_id", "condition"))
rescue <- merge(rescue, normalization, by = "sample_id")
rescue[, selected := baseline_selected + lower_selected - upper_selected]
for (metric in c(
    "supported_anti", "supported_igg", "depth_anti", "depth_igg",
    "same_anti", "opposite_anti", "supported_same_anti",
    "supported_opposite_anti", "signed_distance_sum_anti",
    "absolute_distance_sum_anti"
)) {
    rescue[, (metric) :=
        get(paste0("baseline_", metric)) +
        get(paste0("lower_", metric)) -
        get(paste0("upper_", metric))]
}
rescue <- calculate_ratios(rescue)

gate <- attach_matched_comparison(gate)
rescue <- attach_matched_comparison(rescue)

sample_output_columns <- c(
    "sample_id", "condition", "context_threshold", "selected",
    "matched_tp73_threshold", "matched_tp73_selected",
    "tp73_comparator_lower_threshold", "tp73_comparator_upper_threshold",
    "tp73_comparator_lower_count", "tp73_comparator_upper_count",
    "support_ratio", "depth_ratio",
    "tp73_only_support_ratio", "tp73_only_depth_ratio",
    "incremental_support_ratio", "incremental_depth_ratio"
)
fwrite(gate[, ..sample_output_columns],
       paste0(output_prefix, "_gate_sample_curve.tsv"), sep = "\t")
fwrite(rescue[, ..sample_output_columns],
       paste0(output_prefix, "_rescue_sample_curve.tsv"), sep = "\t")

summarize_rule <- function(table, rule) {
    result <- table[, .(
        selected = selected[[1L]],
        matched_tp73_threshold = matched_tp73_threshold[[1L]],
        matched_tp73_selected = matched_tp73_selected[[1L]],
        minimum_support_ratio = min(support_ratio),
        median_support_ratio = median(support_ratio),
        maximum_support_ratio = max(support_ratio),
        minimum_depth_ratio = min(depth_ratio),
        median_depth_ratio = median(depth_ratio),
        maximum_depth_ratio = max(depth_ratio),
        minimum_incremental_support_ratio = min(incremental_support_ratio),
        median_incremental_support_ratio = median(incremental_support_ratio),
        maximum_incremental_support_ratio = max(incremental_support_ratio),
        minimum_incremental_depth_ratio = min(incremental_depth_ratio),
        median_incremental_depth_ratio = median(incremental_depth_ratio),
        maximum_incremental_depth_ratio = max(incremental_depth_ratio),
        same_orientation_fraction =
            (same_anti[[1L]] / selected[[1L]]),
        median_supported_same_orientation_fraction = median(
            supported_same_anti / supported_anti, na.rm = TRUE
        ),
        mean_signed_center_distance_bp =
            signed_distance_sum_anti[[1L]] / selected[[1L]],
        mean_absolute_center_distance_bp =
            absolute_distance_sum_anti[[1L]] / selected[[1L]]
    ), by = context_threshold]
    result[, rule := rule]
    setcolorder(result, c("rule", setdiff(names(result), "rule")))
    setorder(result, context_threshold)
    result
}

gate_summary <- summarize_rule(gate, "gate")
rescue_summary <- summarize_rule(rescue, "rescue")
fwrite(gate_summary, paste0(output_prefix, "_gate_curve.tsv"), sep = "\t")
fwrite(rescue_summary, paste0(output_prefix, "_rescue_curve.tsv"), sep = "\t")

nearest_row <- function(table, target) {
    table[which.min(abs(context_threshold - target))]
}
fraction_row <- function(table, baseline_count, fraction) {
    table[which.min(abs(as.numeric(selected) - baseline_count * fraction))]
}

baseline_count <- as.numeric(gate_summary[selected == max(selected), selected[[1L]]])
selected_rows <- rbindlist(list(
    transform(nearest_row(gate_summary, minimum_context_threshold),
              selection = paste0("TP73 >= 0; no effective ", context_label, " gate")),
    transform(nearest_row(gate_summary, 0),
              selection = paste0("TP73 >= 0; best ", context_label, " >= 0")),
    transform(fraction_row(gate_summary, baseline_count, 0.50),
              selection = paste0("TP73 >= 0; ", context_label,
                                 " gate retaining about 50%")),
    transform(fraction_row(gate_summary, baseline_count, 0.25),
              selection = paste0("TP73 >= 0; ", context_label,
                                 " gate retaining about 25%")),
    transform(nearest_row(rescue_summary, 0),
              selection = paste0("TP73 >= 0 plus [-2,0) rescued by ",
                                 context_label, " >= 0"))
), use.names = TRUE, fill = TRUE)

consistent_gate <- gate_summary[
    selected >= 1000 & minimum_incremental_support_ratio >= 1 &
        minimum_incremental_depth_ratio >= 1
]
if (nrow(consistent_gate) > 0L) {
    consistent_gate[, score := pmin(
        median_incremental_support_ratio, median_incremental_depth_ratio
    )]
    selected_rows <- rbindlist(list(
        selected_rows,
        transform(consistent_gate[which.max(score)][, score := NULL],
                  selection = paste0("Exploratory strongest all-sample ",
                                     context_label, " gate"))
    ), use.names = TRUE, fill = TRUE)
}
detrimental_gate <- gate_summary[
    selected >= 1000 & maximum_incremental_support_ratio <= 1 &
        maximum_incremental_depth_ratio <= 1
]
if (nrow(detrimental_gate) > 0L) {
    detrimental_gate[, score := pmax(
        median_incremental_support_ratio, median_incremental_depth_ratio
    )]
    selected_rows <- rbindlist(list(
        selected_rows,
        transform(detrimental_gate[which.min(score)][, score := NULL],
                  selection = paste0("Exploratory strongest all-sample detrimental ",
                                     context_label, " gate"))
    ), use.names = TRUE, fill = TRUE)
}
consistent_rescue <- rescue_summary[
    selected >= 1000 & minimum_incremental_support_ratio >= 1 &
        minimum_incremental_depth_ratio >= 1
]
if (nrow(consistent_rescue) > 0L) {
    consistent_rescue[, score := pmin(
        median_incremental_support_ratio, median_incremental_depth_ratio
    )]
    selected_rows <- rbindlist(list(
        selected_rows,
        transform(consistent_rescue[which.max(score)][, score := NULL],
                  selection = paste0("Exploratory strongest all-sample ",
                                     context_label, " rescue"))
    ), use.names = TRUE, fill = TRUE)
}
detrimental_rescue <- rescue_summary[
    selected >= 1000 & maximum_incremental_support_ratio <= 1 &
        maximum_incremental_depth_ratio <= 1
]
if (nrow(detrimental_rescue) > 0L) {
    detrimental_rescue[, score := pmax(
        median_incremental_support_ratio, median_incremental_depth_ratio
    )]
    selected_rows <- rbindlist(list(
        selected_rows,
        transform(detrimental_rescue[which.min(score)][, score := NULL],
                  selection = paste0("Exploratory strongest all-sample detrimental ",
                                     context_label, " rescue"))
    ), use.names = TRUE, fill = TRUE)
}
selected_rows[, minimum_effect_reference := values$minimum_effect]
setcolorder(selected_rows, c("selection", setdiff(names(selected_rows), "selection")))
fwrite(selected_rows, paste0(output_prefix, "_selected.tsv"), sep = "\t")

plot_source <- rbindlist(list(gate_summary, rescue_summary), use.names = TRUE)
plot_source <- plot_source[selected >= 1000 & context_threshold >= -10]
rule_labels <- c(
    gate = "Gate among TP73 >= 0",
    rescue = "Rescue TP73 [-2,0)"
)
plot_source[, rule_label := factor(
    unname(rule_labels[rule]), levels = unname(rule_labels)
)]

absolute_plot <- rbindlist(list(
    plot_source[, .(rule_label, context_threshold, metric = "Strict immersion",
                    statistic = "Median", value = median_support_ratio)],
    plot_source[, .(rule_label, context_threshold, metric = "Strict immersion",
                    statistic = "Least-enriched sample", value = minimum_support_ratio)],
    plot_source[, .(rule_label, context_threshold, metric = "Maximum depth",
                    statistic = "Median", value = median_depth_ratio)],
    plot_source[, .(rule_label, context_threshold, metric = "Maximum depth",
                    statistic = "Least-enriched sample", value = minimum_depth_ratio)]
))
absolute_plot[, metric := factor(
    metric, levels = c("Strict immersion", "Maximum depth")
)]
absolute <- ggplot(
    absolute_plot,
    aes(x = context_threshold, y = value, color = rule_label, linetype = statistic)
) +
    geom_hline(yintercept = 1, color = "#777777", linewidth = 0.4) +
    geom_hline(yintercept = 1 + values$minimum_effect, color = "#777777",
               linewidth = 0.4, linetype = "dotted") +
    geom_line(linewidth = 0.8) +
    facet_wrap(~metric, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c("#0072B2", "#009E73"), name = NULL) +
    scale_linetype_manual(values = c("Median" = "solid",
                                     "Least-enriched sample" = "dashed"),
                          name = NULL) +
    labs(
        title = paste0("TP73 CUT&RUN evidence after adding the strongest nearby ",
                       context_label, " score"),
        subtitle = paste0("Six TA/DN anti-p73 tracks normalized to matched IgG; ",
                          context_label, " within +/-150 bp"),
        x = paste0("Minimum score of the strongest nearby ", context_label,
                   " alignment"),
        y = "Anti-p73 / IgG normalized enrichment"
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA)
    )
ggsave(paste0(output_prefix, "_absolute_effect.png"), absolute,
       width = 10, height = 5.3, dpi = 170)

incremental_plot <- rbindlist(list(
    plot_source[, .(rule_label, context_threshold, metric = "Strict immersion",
                    statistic = "Median", value = median_incremental_support_ratio)],
    plot_source[, .(rule_label, context_threshold, metric = "Strict immersion",
                    statistic = "Minimum across samples",
                    value = minimum_incremental_support_ratio)],
    plot_source[, .(rule_label, context_threshold, metric = "Strict immersion",
                    statistic = "Maximum across samples",
                    value = maximum_incremental_support_ratio)],
    plot_source[, .(rule_label, context_threshold, metric = "Maximum depth",
                    statistic = "Median", value = median_incremental_depth_ratio)],
    plot_source[, .(rule_label, context_threshold, metric = "Maximum depth",
                    statistic = "Minimum across samples",
                    value = minimum_incremental_depth_ratio)],
    plot_source[, .(rule_label, context_threshold, metric = "Maximum depth",
                    statistic = "Maximum across samples",
                    value = maximum_incremental_depth_ratio)]
))
incremental_plot[, metric := factor(
    metric, levels = c("Strict immersion", "Maximum depth")
)]
incremental <- ggplot(
    incremental_plot,
    aes(x = context_threshold, y = value, color = rule_label, linetype = statistic)
) +
    geom_hline(yintercept = 1, color = "#777777", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~metric, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c("#0072B2", "#009E73"), name = NULL) +
    scale_linetype_manual(values = c(
        "Median" = "solid", "Minimum across samples" = "dashed",
        "Maximum across samples" = "dotdash"
    ), name = NULL) +
    labs(
        title = paste0("How does ", context_label,
                       " change TP73 prediction at the same storage size?"),
        subtitle = "TP73-only comparator is count-matched by interpolation within one tied 0.2-score bin",
        x = paste0("Minimum score of the strongest nearby ", context_label,
                   " alignment"),
        y = "Joint rule / matched-size TP73-only rule"
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA)
    )
ggsave(paste0(output_prefix, "_matched_storage.png"), incremental,
       width = 10, height = 5.3, dpi = 170, bg = "white")

writeLines(c(
    paste0("I: TP73 >= ", values$anchor_threshold,
           " baseline starts: ", format(baseline_count, scientific = FALSE)),
    paste0("I: Wrote ", paste0(output_prefix, "_gate_curve.tsv")),
    paste0("I: Wrote ", paste0(output_prefix, "_rescue_curve.tsv")),
    paste0("I: Wrote ", paste0(output_prefix, "_gate_sample_curve.tsv")),
    paste0("I: Wrote ", paste0(output_prefix, "_rescue_sample_curve.tsv")),
    paste0("I: Wrote ", paste0(output_prefix, "_selected.tsv")),
    paste0("I: Wrote ", paste0(output_prefix, "_absolute_effect.png")),
    paste0("I: Wrote ", paste0(output_prefix, "_matched_storage.png"))
), con = stderr())
