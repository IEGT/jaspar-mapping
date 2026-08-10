#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: analyze_tp73_cofactor_enrichment.R --anchor-evidence FILE",
        "       --cofactor-maxima FILE --thresholds FILE --output-prefix PATH",
        "       [options]",
        "",
        "Describe cofactor enrichment/depletion across TP73-score and normalized",
        "CUT&RUN-depth strata, and fit one matched anti/control occupancy test per",
        "motif. Positive cofactor calls use score >= T. Negative references use",
        "strict score < N (including an absent sparse hit); intermediate scores",
        "are excluded from the corresponding contrast.",
        "",
        "Required:",
        "  --anchor-evidence FILE  TP73 anchor/CUT&RUN feature Parquet",
        "  --cofactor-maxima FILE  Rectangular schema-v4 cofactor maxima Parquet",
        "  --thresholds FILE       TSV: motif_id and positive_threshold",
        "  --output-prefix PATH    Basename for generated TSV files",
        "",
        "Options:",
        "  --negative-reference-thresholds LIST",
        "                          Cumulative strict negative references",
        "                          (default: -1,0)",
        "  --primary-negative-reference N",
        "                          Reference used for the one primary occupancy",
        "                          test per motif (default: -1)",
        "  --tp73-score-breaks LIST",
        "                          Half-open stratum boundaries (default:",
        "                          -5,-1,0,1,2,5,10,15,Inf)",
        "  --depth-quantiles LIST  Positive-depth quantiles in addition to strict",
        "                          immersion (default: 0.5,0.75,0.9,0.95,0.99)",
        "  --block-size BP         Genomic cluster size for primary standard",
        "                          errors (default: 5000000)",
        "  --spline-df N           TP73 score spline degrees of freedom (default: 4)",
        "  --minimum-class-fraction N",
        "                          Flag/skip a primary class below this fraction",
        "                          (default: 0.01)",
        "  --duckdb COMMAND        DuckDB CLI executable (default: duckdb)",
        "  --source-commit HEX     Source commit supplied by the orchestrator",
        "                          (default: unknown)",
        "  --source-dirty BOOL     Whether tracked/staged source differed",
        "                          (default: false)",
        "  -h, --help              Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) {
    usage()
}

values <- list(
    negative_reference_thresholds = "-1,0",
    primary_negative_reference = -1,
    tp73_score_breaks = "-5,-1,0,1,2,5,10,15,Inf",
    depth_quantiles = "0.5,0.75,0.9,0.95,0.99",
    block_size = 5000000,
    spline_df = 4L,
    minimum_class_fraction = 0.01,
    duckdb = "duckdb",
    source_commit = "unknown",
    source_dirty = "false"
)
known <- c(
    "--anchor-evidence", "--cofactor-maxima", "--thresholds",
    "--output-prefix", "--negative-reference-thresholds",
    "--primary-negative-reference", "--tp73-score-breaks",
    "--depth-quantiles", "--block-size", "--spline-df",
    "--minimum-class-fraction", "--duckdb", "--source-commit",
    "--source-dirty"
)
index <- 1L
while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (!option %in% known) {
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

required <- c(
    "anchor_evidence", "cofactor_maxima", "thresholds", "output_prefix"
)
missing <- required[vapply(required, function(name) is.null(values[[name]]),
                           logical(1))]
if (length(missing) > 0L) {
    writeLines(paste(
        "E: missing required options:",
        paste(paste0("--", gsub("_", "-", missing)), collapse = ", ")
    ), con = stderr())
    usage(2L)
}

parse_numbers <- function(specification, allow_infinite = FALSE) {
    parsed <- suppressWarnings(as.numeric(strsplit(
        specification, ",", fixed = TRUE
    )[[1L]]))
    valid <- if (allow_infinite) !is.na(parsed) else is.finite(parsed)
    if (length(parsed) == 0L || any(!valid)) {
        stop("invalid numeric list: ", specification, call. = FALSE)
    }
    parsed
}

negative_references <- sort(unique(parse_numbers(
    values$negative_reference_thresholds
)))
values$primary_negative_reference <- suppressWarnings(as.numeric(
    values$primary_negative_reference
))
if (!is.finite(values$primary_negative_reference) ||
    !values$primary_negative_reference %in% negative_references) {
    stop(
        "--primary-negative-reference must occur in ",
        "--negative-reference-thresholds", call. = FALSE
    )
}
tp73_breaks <- parse_numbers(values$tp73_score_breaks, allow_infinite = TRUE)
if (length(tp73_breaks) < 2L || any(diff(tp73_breaks) <= 0) ||
    any(is.infinite(tp73_breaks[-length(tp73_breaks)])) ||
    !is.infinite(tail(tp73_breaks, 1L)) || tail(tp73_breaks, 1L) < 0) {
    stop(
        "--tp73-score-breaks must be strictly increasing, with only the final ",
        "boundary equal to Inf", call. = FALSE
    )
}
depth_quantiles <- sort(unique(parse_numbers(values$depth_quantiles)))
if (any(depth_quantiles <= 0 | depth_quantiles >= 1)) {
    stop("--depth-quantiles values must be in (0,1)", call. = FALSE)
}
values$block_size <- suppressWarnings(as.numeric(values$block_size))
values$spline_df <- suppressWarnings(as.integer(values$spline_df))
values$minimum_class_fraction <- suppressWarnings(as.numeric(
    values$minimum_class_fraction
))
if (!is.finite(values$block_size) || values$block_size <= 300) {
    stop("--block-size must exceed the 150 bp context radius", call. = FALSE)
}
if (is.na(values$spline_df) || values$spline_df < 2L) {
    stop("--spline-df must be at least 2", call. = FALSE)
}
if (!is.finite(values$minimum_class_fraction) ||
    values$minimum_class_fraction < 0 ||
    values$minimum_class_fraction >= 0.5) {
    stop("--minimum-class-fraction must be in [0,0.5)", call. = FALSE)
}
if (!grepl("^([0-9a-f]{40}|unknown)$", values$source_commit)) {
    stop("--source-commit must be 40 lowercase hex digits or unknown", call. = FALSE)
}
if (!values$source_dirty %in% c("true", "false")) {
    stop("--source-dirty must be true or false", call. = FALSE)
}
values$source_dirty <- values$source_dirty == "true"

suppressPackageStartupMessages(library(data.table))

values$anchor_evidence <- normalizePath(values$anchor_evidence, mustWork = TRUE)
values$cofactor_maxima <- normalizePath(values$cofactor_maxima, mustWork = TRUE)
values$thresholds <- normalizePath(values$thresholds, mustWork = TRUE)
duckdb_path <- Sys.which(values$duckdb)
if (!nzchar(duckdb_path)) {
    stop("DuckDB CLI not found: ", values$duckdb, call. = FALSE)
}
output_prefix <- values$output_prefix
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

sql_string <- function(value) {
    paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
}

sql_identifier <- function(value) {
    paste0('"', gsub('"', '""', value, fixed = TRUE), '"')
}

duckdb_fread <- function(query) {
    sql_path <- tempfile("tp73-enrichment-", fileext = ".sql")
    on.exit(unlink(sql_path), add = TRUE)
    writeLines(query, sql_path, useBytes = TRUE)
    command <- paste(
        shQuote(duckdb_path),
        "-light-mode -csv -nullvalue NA :memory: <", shQuote(sql_path)
    )
    fread(cmd = command, na.strings = "NA", showProgress = FALSE)
}

thresholds <- fread(values$thresholds, sep = "\t", na.strings = "NA")
expected_threshold_columns <- c("motif_id", "positive_threshold")
missing_threshold_columns <- setdiff(expected_threshold_columns, names(thresholds))
if (length(missing_threshold_columns) > 0L) {
    stop(
        "threshold table lacks columns: ",
        paste(missing_threshold_columns, collapse = ", "), call. = FALSE
    )
}
for (optional_column in c(
    "factor_name", "positive_threshold_source", "selection_semantics"
)) {
    if (!optional_column %in% names(thresholds)) {
        thresholds[, (optional_column) := NA_character_]
    }
}
thresholds <- thresholds[, .(
    motif_id = as.character(motif_id),
    factor_name = as.character(factor_name),
    positive_threshold = as.numeric(positive_threshold),
    positive_threshold_source = as.character(positive_threshold_source),
    selection_semantics = as.character(selection_semantics)
)]
if (nrow(thresholds) == 0L || any(!nzchar(thresholds$motif_id)) ||
    any(!is.finite(thresholds$positive_threshold)) ||
    anyDuplicated(thresholds$motif_id)) {
    stop("threshold table is empty, duplicated, or non-finite", call. = FALSE)
}
if (any(thresholds$positive_threshold < max(negative_references))) {
    stop(
        "every positive threshold must be at least the largest negative ",
        "reference threshold", call. = FALSE
    )
}
setorder(thresholds, motif_id)

description <- duckdb_fread(paste0(
    "DESCRIBE SELECT * FROM read_parquet(",
    sql_string(values$anchor_evidence), ");"
))
anchor_columns <- as.character(description$column_name)
base_anchor_columns <- c("chrom", "anchor_start", "anchor_end", "anchor_score")
missing_anchor_columns <- setdiff(base_anchor_columns, anchor_columns)
if (length(missing_anchor_columns) > 0L) {
    stop(
        "anchor evidence lacks columns: ",
        paste(missing_anchor_columns, collapse = ", "), call. = FALSE
    )
}
anti_support_columns <- sort(grep(
    "^supported_anti_", anchor_columns, value = TRUE
))
if (length(anti_support_columns) > 0L) {
    sample_ids <- sub("^supported_anti_", "", anti_support_columns)
    control_support_columns <- paste0("supported_control_", sample_ids)
    anti_depth_columns <- paste0("depth_anti_", sample_ids)
    control_depth_columns <- paste0("depth_control_", sample_ids)
    evidence_column_scheme <- "supported_anti_and_control"
} else {
    anti_support_columns <- sort(grep(
        "^supported_tp73_", anchor_columns, value = TRUE
    ))
    sample_ids <- sub("^supported_tp73_", "", anti_support_columns)
    control_support_columns <- paste0(
        "supported_negative_control_", sample_ids
    )
    anti_depth_columns <- paste0("depth_tp73_", sample_ids)
    control_depth_columns <- paste0("depth_negative_control_", sample_ids)
    evidence_column_scheme <- "supported_tp73_and_negative_control"
}
if (length(sample_ids) == 0L || anyDuplicated(sample_ids)) {
    stop("anchor evidence contains no unique anti-p73 samples", call. = FALSE)
}
samples <- data.table(
    sample_id = sample_ids,
    anti_support = anti_support_columns,
    control_support = control_support_columns,
    anti_depth = anti_depth_columns,
    control_depth = control_depth_columns
)
required_evidence_columns <- unique(c(
    base_anchor_columns, samples$anti_support, samples$control_support,
    samples$anti_depth, samples$control_depth
))
missing_evidence_columns <- setdiff(required_evidence_columns, anchor_columns)
if (length(missing_evidence_columns) > 0L) {
    stop(
        "anchor evidence lacks matched support/depth columns: ",
        paste(missing_evidence_columns, collapse = ", "), call. = FALSE
    )
}

anchor_select <- c(
    "CAST(chrom AS VARCHAR) AS chrom", "anchor_start", "anchor_end",
    "anchor_score", vapply(
        setdiff(required_evidence_columns, base_anchor_columns),
        sql_identifier, character(1)
    )
)
message("I: reading TP73 anchors and matched CUT&RUN support/depth")
anchors <- duckdb_fread(paste0(
    "SELECT ", paste(anchor_select, collapse = ", "), "\n",
    "FROM read_parquet(", sql_string(values$anchor_evidence), ")\n",
    "ORDER BY chrom, anchor_start, anchor_end;"
))
if (nrow(anchors) == 0L || anchors[, anyDuplicated(
    paste(chrom, anchor_start, anchor_end)
)] != 0L) {
    stop("anchor evidence is empty or has duplicate spans", call. = FALSE)
}
if (anchors[, any(!is.finite(anchor_score))]) {
    stop("anchor evidence contains non-finite TP73 scores", call. = FALSE)
}
if (min(anchors$anchor_score) < tp73_breaks[[1L]] ||
    max(anchors$anchor_score) >= tail(tp73_breaks, 1L)) {
    stop("TP73 score breaks do not cover all anchors", call. = FALSE)
}
for (sample_index in seq_len(nrow(samples))) {
    for (antibody in c("anti", "control")) {
        support_column <- samples[[paste0(antibody, "_support")]][[sample_index]]
        depth_column <- samples[[paste0(antibody, "_depth")]][[sample_index]]
        support <- as.logical(anchors[[support_column]])
        depth <- as.numeric(anchors[[depth_column]])
        if (anyNA(support) || any(!is.finite(depth)) || any(depth < 0)) {
            stop("invalid support/depth values in sample ",
                 samples$sample_id[[sample_index]], call. = FALSE)
        }
        if (any(support != (depth > 0))) {
            stop(
                "support and effective depth disagree for ", antibody, " ",
                samples$sample_id[[sample_index]], call. = FALSE
            )
        }
        set(anchors, j = support_column, value = support)
        set(anchors, j = depth_column, value = depth)
    }
}

motif_list <- paste(vapply(
    thresholds$motif_id, sql_string, character(1)
), collapse = ",")
message("I: reading interval-defined cofactor maxima")
maxima <- duckdb_fread(paste0(
    "SELECT CAST(chrom AS VARCHAR) AS chrom, anchor_start, anchor_end, ",
    "motif_id, context_score, source_score_floor, context_flank_bp, ",
    "context_distance_metric\n",
    "FROM read_parquet(", sql_string(values$cofactor_maxima), ")\n",
    "WHERE motif_id IN (", motif_list, ")\n",
    "ORDER BY motif_id, chrom, anchor_start, anchor_end;"
))
expected_maxima_columns <- c(
    "chrom", "anchor_start", "anchor_end", "motif_id", "context_score",
    "source_score_floor", "context_flank_bp", "context_distance_metric"
)
if (nrow(maxima) == 0L ||
    length(setdiff(expected_maxima_columns, names(maxima))) > 0L) {
    stop("cofactor maxima are empty or incomplete", call. = FALSE)
}
if (!setequal(unique(maxima$motif_id), thresholds$motif_id)) {
    stop("threshold motifs and cofactor-maxima motifs differ", call. = FALSE)
}
if (nrow(maxima) != nrow(anchors) * nrow(thresholds) ||
    maxima[, anyDuplicated(paste(chrom, anchor_start, anchor_end, motif_id))] != 0L) {
    stop("cofactor maxima are not rectangular and unique", call. = FALSE)
}
if (!all(maxima$context_distance_metric == "signed_interval_edge_distance")) {
    stop("cofactor maxima do not use schema-v4 interval geometry", call. = FALSE)
}
if (any(maxima$context_flank_bp != 150)) {
    stop("this analysis requires the prespecified 150 bp context", call. = FALSE)
}
source_floors <- maxima[, .(
    n_source_floors = uniqueN(source_score_floor),
    source_score_floor = suppressWarnings(min(source_score_floor))
), by = motif_id]
if (source_floors[, any(n_source_floors != 1L)] ||
    source_floors[, any(!is.finite(source_score_floor))] ||
    source_floors[, any(source_score_floor > min(negative_references))]) {
    stop(
        "each motif needs one source floor no greater than the smallest ",
        "negative reference", call. = FALSE
    )
}

format_boundary <- function(value) {
    if (is.infinite(value)) "+Inf" else format(value, trim = TRUE, scientific = FALSE)
}

score_comparison_label <- function(positive_threshold, negative_reference) {
    paste0(
        "score>=", format_boundary(positive_threshold),
        " versus score<", format_boundary(negative_reference)
    )
}
tp73_labels <- vapply(seq_len(length(tp73_breaks) - 1L), function(i) {
    paste0(
        "[", format_boundary(tp73_breaks[[i]]), ",",
        format_boundary(tp73_breaks[[i + 1L]]), ")"
    )
}, character(1))
anchors[, tp73_score_stratum := as.character(cut(
    anchor_score, breaks = tp73_breaks, labels = tp73_labels,
    right = FALSE, include.lowest = TRUE
))]
if (anchors[, anyNA(tp73_score_stratum)]) {
    stop("failed to assign a TP73 score stratum", call. = FALSE)
}

depth_manifest_rows <- list()
manifest_index <- 1L
for (sample_index in seq_len(nrow(samples))) {
    for (antibody in c("anti", "control")) {
        depth_column <- samples[[paste0(antibody, "_depth")]][[sample_index]]
        depth <- anchors[[depth_column]]
        positive_depth <- depth[depth > 0]
        if (length(positive_depth) == 0L) {
            stop("sample has no positive effective depths: ",
                 samples$sample_id[[sample_index]], call. = FALSE)
        }
        rows <- list(data.table(
            sample_id = samples$sample_id[[sample_index]],
            antibody = antibody,
            depth_tier = "strict_immersion",
            depth_tier_order = 0L,
            positive_depth_quantile = NA_real_,
            effective_depth_cutoff = 0,
            anchors_meeting_tier = sum(depth > 0),
            anchor_fraction_meeting_tier = mean(depth > 0),
            supported_fraction_meeting_tier = 1
        ))
        previous_event <- depth > 0
        for (quantile_index in seq_along(depth_quantiles)) {
            probability <- depth_quantiles[[quantile_index]]
            cutoff <- as.numeric(quantile(
                positive_depth, probability, names = FALSE, type = 7
            ))
            event <- depth > 0 & depth >= cutoff
            rows[[length(rows) + 1L]] <- data.table(
                sample_id = samples$sample_id[[sample_index]],
                antibody = antibody,
                depth_tier = paste0(
                    "positive_depth_q", sprintf("%02d", round(100 * probability))
                ),
                depth_tier_order = as.integer(quantile_index),
                positive_depth_quantile = probability,
                effective_depth_cutoff = cutoff,
                anchors_meeting_tier = sum(event),
                anchor_fraction_meeting_tier = mean(event),
                supported_fraction_meeting_tier = mean(event[depth > 0]),
                duplicates_previous_tier = all(event == previous_event)
            )
            previous_event <- event
        }
        rows[[1L]][, duplicates_previous_tier := FALSE]
        depth_manifest_rows[[manifest_index]] <- rbindlist(rows, fill = TRUE)
        manifest_index <- manifest_index + 1L
    }
}
depth_manifest <- rbindlist(depth_manifest_rows, use.names = TRUE, fill = TRUE)
setorder(depth_manifest, sample_id, antibody, depth_tier_order)
fwrite(
    depth_manifest, paste0(output_prefix, "_depth_tier_manifest.tsv"), sep = "\t"
)

jeffreys_rate <- function(events, total) {
    if (total <= 0L) NA_real_ else (events + 0.5) / (total + 1)
}

raw_rate <- function(events, total) {
    if (total <= 0L) NA_real_ else events / total
}

finite_mean <- function(value) {
    selected <- value[is.finite(value)]
    if (length(selected) == 0L) NA_real_ else mean(selected)
}

finite_median <- function(value) {
    selected <- value[is.finite(value)]
    if (length(selected) == 0L) NA_real_ else median(selected)
}

descriptive_rows <- list()
class_rows <- list()
descriptive_index <- 1L
class_index <- 1L
primary_inputs <- list()

setkey(anchors, chrom, anchor_start, anchor_end)
setkey(maxima, chrom, anchor_start, anchor_end)
for (motif_index in seq_len(nrow(thresholds))) {
    current_motif <- thresholds$motif_id[[motif_index]]
    factor_name <- thresholds$factor_name[[motif_index]]
    positive_threshold <- thresholds$positive_threshold[[motif_index]]
    positive_threshold_source <-
        thresholds$positive_threshold_source[[motif_index]]
    selection_semantics <- thresholds$selection_semantics[[motif_index]]
    motif_maxima <- maxima[motif_id == current_motif]
    motif_features <- motif_maxima[anchors, on = .(chrom, anchor_start, anchor_end)]
    if (nrow(motif_features) != nrow(anchors)) {
        stop("cofactor join changed anchor cardinality for ", current_motif,
             call. = FALSE)
    }
    score <- motif_features$context_score
    for (negative_reference in negative_references) {
        positive <- !is.na(score) & score >= positive_threshold
        negative <- is.na(score) | score < negative_reference
        intermediate <- !(positive | negative)
        if (any(positive & negative)) {
            stop("positive and negative classes overlap for ", current_motif,
                 call. = FALSE)
        }
        for (stratum in c("all", tp73_labels)) {
            stratum_selected <- if (stratum == "all") {
                rep(TRUE, nrow(anchors))
            } else {
                anchors$tp73_score_stratum == stratum
            }
            class_rows[[class_index]] <- data.table(
                motif_id = current_motif,
                factor_name = factor_name,
                positive_threshold = positive_threshold,
                positive_threshold_source = positive_threshold_source,
                selection_semantics = selection_semantics,
                negative_reference_threshold = negative_reference,
                comparison_label = score_comparison_label(
                    positive_threshold, negative_reference
                ),
                negative_reference_rule = "context_score < N or absent",
                positive_rule = "context_score >= T",
                intermediate_rule = "N <= context_score < T",
                tp73_score_stratum = stratum,
                anchors_total = sum(stratum_selected),
                anchors_positive = sum(stratum_selected & positive),
                anchors_negative_reference = sum(stratum_selected & negative),
                anchors_intermediate = sum(stratum_selected & intermediate)
            )
            class_index <- class_index + 1L
        }

        if (negative_reference == values$primary_negative_reference) {
            primary_inputs[[current_motif]] <- list(
                positive = positive, negative = negative,
                positive_threshold = positive_threshold
            )
        }

        for (sample_index in seq_len(nrow(samples))) {
            current_sample <- samples$sample_id[[sample_index]]
            anti_depth <- anchors[[samples$anti_depth[[sample_index]]]]
            control_depth <- anchors[[samples$control_depth[[sample_index]]]]
            sample_manifest <- depth_manifest[sample_id == current_sample]
            for (tier_index in sort(unique(sample_manifest$depth_tier_order))) {
                anti_manifest <- sample_manifest[
                    antibody == "anti" & depth_tier_order == tier_index
                ]
                control_manifest <- sample_manifest[
                    antibody == "control" & depth_tier_order == tier_index
                ]
                if (nrow(anti_manifest) != 1L || nrow(control_manifest) != 1L) {
                    stop("depth tier manifest is incomplete", call. = FALSE)
                }
                if (tier_index == 0L) {
                    anti_event <- anti_depth > 0
                    control_event <- control_depth > 0
                } else {
                    anti_event <- anti_depth > 0 &
                        anti_depth >= anti_manifest$effective_depth_cutoff
                    control_event <- control_depth > 0 &
                        control_depth >= control_manifest$effective_depth_cutoff
                }
                for (stratum in c("all", tp73_labels)) {
                    stratum_selected <- if (stratum == "all") {
                        rep(TRUE, nrow(anchors))
                    } else {
                        anchors$tp73_score_stratum == stratum
                    }
                    positive_selected <- stratum_selected & positive
                    negative_selected <- stratum_selected & negative
                    intermediate_selected <- stratum_selected & intermediate
                    n_positive <- sum(positive_selected)
                    n_negative <- sum(negative_selected)
                    n_intermediate <- sum(intermediate_selected)
                    anti_positive_events <- sum(anti_event & positive_selected)
                    anti_negative_events <- sum(anti_event & negative_selected)
                    control_positive_events <- sum(
                        control_event & positive_selected
                    )
                    control_negative_events <- sum(
                        control_event & negative_selected
                    )
                    anti_positive_rate <- raw_rate(
                        anti_positive_events, n_positive
                    )
                    anti_negative_rate <- raw_rate(
                        anti_negative_events, n_negative
                    )
                    control_positive_rate <- raw_rate(
                        control_positive_events, n_positive
                    )
                    control_negative_rate <- raw_rate(
                        control_negative_events, n_negative
                    )
                    anti_positive_smoothed <- jeffreys_rate(
                        anti_positive_events, n_positive
                    )
                    anti_negative_smoothed <- jeffreys_rate(
                        anti_negative_events, n_negative
                    )
                    control_positive_smoothed <- jeffreys_rate(
                        control_positive_events, n_positive
                    )
                    control_negative_smoothed <- jeffreys_rate(
                        control_negative_events, n_negative
                    )
                    anti_rr <- anti_positive_smoothed / anti_negative_smoothed
                    control_rr <- control_positive_smoothed /
                        control_negative_smoothed
                    eligible <- n_positive + n_negative
                    class_fraction <- if (eligible == 0L) 0 else
                        min(n_positive, n_negative) / eligible
                    status <- if (n_positive == 0L || n_negative == 0L) {
                        "empty_anchor_class"
                    } else if (class_fraction < values$minimum_class_fraction) {
                        "underpowered_anchor_class"
                    } else {
                        "descriptive"
                    }
                    log2_specificity <- log2(anti_rr / control_rr)
                    descriptive_rows[[descriptive_index]] <- data.table(
                        motif_id = current_motif,
                        factor_name = factor_name,
                        positive_threshold = positive_threshold,
                        positive_threshold_source = positive_threshold_source,
                        selection_semantics = selection_semantics,
                        negative_reference_threshold = negative_reference,
                        comparison_label = score_comparison_label(
                            positive_threshold, negative_reference
                        ),
                        negative_reference_rule = "context_score < N or absent",
                        tp73_score_stratum = stratum,
                        sample_id = current_sample,
                        depth_tier = anti_manifest$depth_tier,
                        depth_tier_order = tier_index,
                        anti_effective_depth_cutoff =
                            anti_manifest$effective_depth_cutoff,
                        control_effective_depth_cutoff =
                            control_manifest$effective_depth_cutoff,
                        anchors_positive = n_positive,
                        anchors_negative_reference = n_negative,
                        anchors_intermediate = n_intermediate,
                        anti_positive_events = anti_positive_events,
                        anti_negative_events = anti_negative_events,
                        control_positive_events = control_positive_events,
                        control_negative_events = control_negative_events,
                        anti_positive_probability = anti_positive_rate,
                        anti_negative_probability = anti_negative_rate,
                        control_positive_probability = control_positive_rate,
                        control_negative_probability = control_negative_rate,
                        anti_relative_risk_jeffreys = anti_rr,
                        control_relative_risk_jeffreys = control_rr,
                        anti_control_specificity_ratio_jeffreys = anti_rr /
                            control_rr,
                        log2_anti_control_specificity_ratio_jeffreys =
                            log2_specificity,
                        anti_control_risk_difference_difference =
                            (anti_positive_rate - anti_negative_rate) -
                            (control_positive_rate - control_negative_rate),
                        direction = if (!is.finite(log2_specificity)) {
                            "not_estimable"
                        } else if (log2_specificity > 0) {
                            "anti_p73_enriched"
                        } else if (log2_specificity < 0) {
                            "anti_p73_depleted"
                        } else {
                            "neutral"
                        },
                        status = status
                    )
                    descriptive_index <- descriptive_index + 1L
                }
            }
        }
    }
    message("I: summarized motif ", current_motif, " (", motif_index, "/",
            nrow(thresholds), ")")
}

class_counts <- rbindlist(class_rows, use.names = TRUE)
descriptive <- rbindlist(descriptive_rows, use.names = TRUE)
setorder(
    class_counts, motif_id, negative_reference_threshold, tp73_score_stratum
)
setorder(
    descriptive, motif_id, negative_reference_threshold, tp73_score_stratum,
    sample_id, depth_tier_order
)
fwrite(class_counts, paste0(output_prefix, "_class_counts.tsv"), sep = "\t")
fwrite(descriptive, paste0(output_prefix, "_descriptive.tsv"), sep = "\t")

macro <- descriptive[, .(
    samples_total = .N,
    samples_estimable = sum(is.finite(
        log2_anti_control_specificity_ratio_jeffreys
    )),
    mean_log2_anti_control_specificity_ratio_jeffreys = finite_mean(
        log2_anti_control_specificity_ratio_jeffreys
    ),
    median_log2_anti_control_specificity_ratio_jeffreys = finite_median(
        log2_anti_control_specificity_ratio_jeffreys
    ),
    minimum_log2_anti_control_specificity_ratio_jeffreys = suppressWarnings(
        min(log2_anti_control_specificity_ratio_jeffreys, na.rm = TRUE)
    ),
    maximum_log2_anti_control_specificity_ratio_jeffreys = suppressWarnings(
        max(log2_anti_control_specificity_ratio_jeffreys, na.rm = TRUE)
    ),
    mean_anti_control_risk_difference_difference = finite_mean(
        anti_control_risk_difference_difference
    ),
    samples_anti_p73_enriched = sum(direction == "anti_p73_enriched"),
    samples_anti_p73_depleted = sum(direction == "anti_p73_depleted"),
    support_status = if (all(status == "descriptive")) {
        "descriptive"
    } else {
        paste(sort(unique(status)), collapse = ",")
    }
), by = .(
    motif_id, factor_name, positive_threshold, positive_threshold_source,
    selection_semantics, negative_reference_threshold, comparison_label,
    negative_reference_rule, tp73_score_stratum, depth_tier, depth_tier_order
)]
macro[!is.finite(
    minimum_log2_anti_control_specificity_ratio_jeffreys
), minimum_log2_anti_control_specificity_ratio_jeffreys := NA_real_]
macro[!is.finite(
    maximum_log2_anti_control_specificity_ratio_jeffreys
), maximum_log2_anti_control_specificity_ratio_jeffreys := NA_real_]
setorder(
    macro, motif_id, negative_reference_threshold, tp73_score_stratum,
    depth_tier_order
)
fwrite(macro, paste0(output_prefix, "_macro_summary.tsv"), sep = "\t")

clustered_retained_test <- function(input) {
    positive <- input$positive
    negative <- input$negative
    eligible <- positive | negative
    eligible_total <- sum(eligible)
    positive_fraction <- if (eligible_total == 0L) 0 else
        sum(positive) / eligible_total
    negative_fraction <- if (eligible_total == 0L) 0 else
        sum(negative) / eligible_total
    status_row <- function(status, note) data.table(
        anchors_total = nrow(anchors),
        anchors_positive = sum(positive),
        anchors_negative_reference = sum(negative),
        anchors_intermediate = sum(!eligible),
        discordant_observations = NA_integer_,
        genomic_blocks = NA_integer_,
        adjusted_log_odds = NA_real_,
        adjusted_odds_ratio = NA_real_,
        block_clustered_standard_error = NA_real_,
        confidence_interval_95_lower = NA_real_,
        confidence_interval_95_upper = NA_real_,
        p_value = NA_real_,
        evaluation_status = status,
        evaluation_note = note
    )
    if (eligible_total == 0L ||
        min(positive_fraction, negative_fraction) <
            values$minimum_class_fraction) {
        return(status_row(
            "underpowered_anchor_class",
            "primary positive or negative class is below the declared fraction"
        ))
    }
    evidence_rows <- lapply(seq_len(nrow(samples)), function(sample_index) {
        anti <- anchors[[samples$anti_support[[sample_index]]]]
        control <- anchors[[samples$control_support[[sample_index]]]]
        selected <- eligible & anti != control
        data.table(
            chrom = anchors$chrom[selected],
            anchor_start = anchors$anchor_start[selected],
            anchor_score = anchors$anchor_score[selected],
            sample_id = samples$sample_id[[sample_index]],
            retained = as.integer(positive[selected]),
            outcome = as.integer(anti[selected])
        )
    })
    evidence <- rbindlist(evidence_rows, use.names = TRUE)
    if (nrow(evidence) == 0L || evidence[, uniqueN(outcome)] != 2L ||
        evidence[, uniqueN(retained)] != 2L) {
        return(status_row(
            "not_estimable", "discordant matched evidence lacks two classes"
        ))
    }
    evidence[, sample_id := factor(sample_id, levels = samples$sample_id)]
    evidence[, genomic_block := paste(
        chrom, floor(anchor_start / values$block_size), sep = ":"
    )]
    formula <- as.formula(paste0(
        "outcome ~ sample_id + splines::ns(anchor_score, df = ",
        values$spline_df, ") + retained"
    ))
    fit <- tryCatch(
        suppressWarnings(glm(
            formula, data = evidence, family = binomial(),
            control = glm.control(maxit = 50L, epsilon = 1e-8),
            x = TRUE, y = TRUE
        )),
        error = function(condition) NULL
    )
    if (is.null(fit) || !isTRUE(fit$converged) ||
        !"retained" %in% names(coef(fit)) ||
        !is.finite(coef(fit)[["retained"]])) {
        return(status_row("not_estimable", "primary GLM did not converge"))
    }
    model_matrix <- fit$x
    bread <- tryCatch(vcov(fit), error = function(condition) NULL)
    if (is.null(bread) || any(!is.finite(bread))) {
        return(status_row("not_estimable", "primary covariance is singular"))
    }
    score_rows <- model_matrix * as.numeric(fit$y - fitted(fit))
    cluster_scores <- rowsum(
        score_rows, group = evidence$genomic_block, reorder = FALSE
    )
    clusters <- nrow(cluster_scores)
    observations <- nrow(evidence)
    parameters <- ncol(model_matrix)
    if (clusters < 2L || observations <= parameters) {
        return(status_row(
            "not_estimable", "too few genomic blocks for clustered uncertainty"
        ))
    }
    correction <- (clusters / (clusters - 1)) *
        ((observations - 1) / (observations - parameters))
    clustered_covariance <- bread %*% crossprod(cluster_scores) %*% bread *
        correction
    coefficient_index <- match("retained", colnames(model_matrix))
    standard_error <- sqrt(clustered_covariance[
        coefficient_index, coefficient_index
    ])
    estimate <- unname(coef(fit)[["retained"]])
    if (!is.finite(standard_error) || standard_error <= 0) {
        return(status_row(
            "not_estimable", "block-clustered standard error is invalid"
        ))
    }
    z_value <- estimate / standard_error
    data.table(
        anchors_total = nrow(anchors),
        anchors_positive = sum(positive),
        anchors_negative_reference = sum(negative),
        anchors_intermediate = sum(!eligible),
        discordant_observations = observations,
        genomic_blocks = clusters,
        adjusted_log_odds = estimate,
        adjusted_odds_ratio = exp(estimate),
        block_clustered_standard_error = standard_error,
        confidence_interval_95_lower = exp(estimate - 1.96 * standard_error),
        confidence_interval_95_upper = exp(estimate + 1.96 * standard_error),
        p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE),
        evaluation_status = "ok",
        evaluation_note = paste(
            "discordant-pair conditioning estimates the cofactor-by-antibody",
            "interaction; SE clusters genomic blocks of",
            format(values$block_size, scientific = FALSE), "bp"
        )
    )
}

message("I: fitting one primary matched occupancy test per motif")
primary_rows <- lapply(seq_len(nrow(thresholds)), function(motif_index) {
    motif_id <- thresholds$motif_id[[motif_index]]
    result <- clustered_retained_test(primary_inputs[[motif_id]])
    result[, `:=`(
        motif_id = motif_id,
        factor_name = thresholds$factor_name[[motif_index]],
        positive_threshold = thresholds$positive_threshold[[motif_index]],
        positive_threshold_source =
            thresholds$positive_threshold_source[[motif_index]],
        selection_semantics = thresholds$selection_semantics[[motif_index]],
        negative_reference_threshold = values$primary_negative_reference,
        comparison_label = score_comparison_label(
            thresholds$positive_threshold[[motif_index]],
            values$primary_negative_reference
        ),
        negative_reference_rule = "context_score < N or absent",
        outcome = "strict_immersion",
        matched_design = "anti_p73_vs_its_matched_control",
        adjustment = paste0(
            "sample fixed effect + TP73 natural spline(df=",
            values$spline_df, ")"
        ),
        unavailable_covariates = paste(
            "GC,mappability,repeat_class,accessibility,GFP",
            "not present in anchor evidence", sep = ":"
        )
    )]
    result
})
primary <- rbindlist(primary_rows, use.names = TRUE, fill = TRUE)
primary[, q_value_bh_panel := NA_real_]
estimable <- which(primary$evaluation_status == "ok" &
                   is.finite(primary$p_value))
if (length(estimable) > 0L) {
    primary$q_value_bh_panel[estimable] <- p.adjust(
        primary$p_value[estimable], method = "BH"
    )
}
primary[, `:=`(
    requested_motifs_in_multiple_testing_scope = nrow(thresholds),
    estimable_motifs_in_multiple_testing_scope = length(estimable),
    multiple_testing_scope = "motifs in this threshold input; exploratory"
)]
setcolorder(primary, c(
    "motif_id", "factor_name", "positive_threshold",
    "positive_threshold_source", "selection_semantics",
    "negative_reference_threshold", "comparison_label",
    "negative_reference_rule", "outcome", "matched_design", "adjustment",
    "anchors_total", "anchors_positive", "anchors_negative_reference",
    "anchors_intermediate", "discordant_observations", "genomic_blocks",
    "adjusted_log_odds", "adjusted_odds_ratio",
    "block_clustered_standard_error", "confidence_interval_95_lower",
    "confidence_interval_95_upper", "p_value", "q_value_bh_panel",
    "requested_motifs_in_multiple_testing_scope",
    "estimable_motifs_in_multiple_testing_scope", "multiple_testing_scope",
    "unavailable_covariates", "evaluation_status", "evaluation_note"
))
setorder(primary, motif_id)
fwrite(primary, paste0(output_prefix, "_primary_occupancy.tsv"), sep = "\t")

run_config <- data.table(
    anchor_evidence = values$anchor_evidence,
    cofactor_maxima = values$cofactor_maxima,
    thresholds = values$thresholds,
    motifs = paste(thresholds$motif_id, collapse = ","),
    evidence_column_scheme = evidence_column_scheme,
    sample_ids = paste(samples$sample_id, collapse = ","),
    sample_count = nrow(samples),
    negative_reference_thresholds = paste(negative_references, collapse = ","),
    primary_negative_reference = values$primary_negative_reference,
    negative_reference_semantics = "strict context_score < N or absent",
    positive_semantics = "inclusive context_score >= T",
    intermediate_semantics = "excluded N <= context_score < T",
    context_flank_bp = 150,
    context_geometry = "signed_interval_edge_distance",
    tp73_score_breaks = paste(tp73_breaks, collapse = ","),
    depth_tiers = paste(
        c("strict_immersion", paste0("positive_depth_q", sprintf(
            "%02d", round(100 * depth_quantiles)
        ))), collapse = ","
    ),
    depth_normalization = "positive effective-depth quantiles within each track",
    descriptive_smoothing = "Jeffreys beta(0.5,0.5) event probabilities",
    primary_estimand = paste(
        "cofactor-by-antibody interaction via discordant matched-pair",
        "conditioning"
    ),
    primary_formula = paste0(
        "outcome ~ sample_id + ns(anchor_score,df=", values$spline_df,
        ") + retained"
    ),
    block_size_bp = values$block_size,
    spline_df = values$spline_df,
    minimum_class_fraction = values$minimum_class_fraction,
    source_commit = values$source_commit,
    source_dirty = values$source_dirty,
    source_dirty_scope = "tracked_and_staged_files_only",
    inference_status = paste(
        "exploratory chr1; selected thresholds; missing GC/mappability/repeat/",
        "accessibility/GFP and independent validation", sep = ""
    )
)
fwrite(run_config, paste0(output_prefix, "_run_config.tsv"), sep = "\t")

message("I: wrote cofactor enrichment outputs with prefix ", output_prefix)
