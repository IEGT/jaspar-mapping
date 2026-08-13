#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: analyze_tp73_cofactor_score_surface.R --anchor-evidence FILE",
        "       --cofactor-maxima FILE --thresholds FILE --output-prefix PATH",
        "       [options]",
        "",
        "Estimate how the matched anti-p73/control occupancy association changes",
        "jointly with TP73 score and ordered cofactor-score strength. Cofactor",
        "scores are divided into empirical bands above a fixed strict negative",
        "reference; no-hit anchors are negative only when the sparse source floor",
        "makes that classification observable.",
        "",
        "Required:",
        "  --anchor-evidence FILE  TP73 anchor/CUT&RUN evidence Parquet",
        "  --cofactor-maxima FILE  Rectangular schema-v4 context-maxima Parquet",
        "  --thresholds FILE       TSV with motif_id and positive_threshold",
        "  --output-prefix PATH    Basename for generated TSV files",
        "",
        "Options:",
        "  --negative-reference N Strict negative score boundary (default: -1)",
        "  --cofactor-quantiles LIST",
        "                         Positive-score band boundaries (default:",
        "                         0.5,0.75,0.9,0.95,0.99)",
        "  --tp73-grid LIST       Fixed TP73 scores for prediction (default:",
        "                         -1,0,1,2,5,10,15)",
        "  --tp73-quantiles LIST  Additional TP73 prediction quantiles (default:",
        "                         0.05,0.25,0.5,0.75,0.95)",
        "  --block-size BP        Cluster width for uncertainty (default: 5000000)",
        "  --spline-df N          TP73 natural-spline degrees of freedom (default: 4)",
        "  --minimum-band-count N Minimum anchors in every modeled band",
        "                         (default: 100)",
        "  --minimum-class-fraction N",
        "                         Minimum negative-reference fraction",
        "                         (default: 0.01)",
        "  --exclude-samples LIST Comma-separated evidence sample IDs to omit",
        "                         (default: none)",
        "  --duckdb COMMAND       DuckDB CLI executable (default: duckdb)",
        "  --source-commit HEX    Source commit supplied by an orchestrator",
        "                         (default: unknown)",
        "  --source-dirty BOOL    Whether tracked/staged source differed",
        "                         (default: false)",
        "  --inference-status TEXT",
        "                         Scientific caveat recorded verbatim",
        "  -h, --help             Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) {
    usage()
}

values <- list(
    negative_reference = -1,
    cofactor_quantiles = "0.5,0.75,0.9,0.95,0.99",
    tp73_grid = "-1,0,1,2,5,10,15",
    tp73_quantiles = "0.05,0.25,0.5,0.75,0.95",
    block_size = 5000000,
    spline_df = 4L,
    minimum_band_count = 100L,
    minimum_class_fraction = 0.01,
    exclude_samples = "",
    duckdb = "duckdb",
    source_commit = "unknown",
    source_dirty = "false",
    inference_status = paste(
        "exploratory score-response surface; motif-score bands are empirical;",
        "matched CUT&RUN strict immersion is a technical occupancy outcome"
    )
)
known <- c(
    "--anchor-evidence", "--cofactor-maxima", "--thresholds",
    "--output-prefix", "--negative-reference", "--cofactor-quantiles",
    "--tp73-grid", "--tp73-quantiles", "--block-size", "--spline-df",
    "--minimum-band-count", "--minimum-class-fraction", "--exclude-samples",
    "--duckdb", "--source-commit", "--source-dirty",
    "--inference-status"
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

required <- c("anchor_evidence", "cofactor_maxima", "thresholds", "output_prefix")
missing <- required[vapply(required, function(name) is.null(values[[name]]), logical(1))]
if (length(missing) > 0L) {
    writeLines(paste(
        "E: missing required options:",
        paste(paste0("--", gsub("_", "-", missing)), collapse = ", ")
    ), con = stderr())
    usage(2L)
}

parse_numbers <- function(specification) {
    parsed <- suppressWarnings(as.numeric(strsplit(
        specification, ",", fixed = TRUE
    )[[1L]]))
    if (length(parsed) == 0L || any(!is.finite(parsed))) {
        stop("invalid numeric list: ", specification, call. = FALSE)
    }
    parsed
}

values$negative_reference <- suppressWarnings(as.numeric(values$negative_reference))
cofactor_quantiles <- sort(unique(parse_numbers(values$cofactor_quantiles)))
fixed_tp73_grid <- sort(unique(parse_numbers(values$tp73_grid)))
tp73_quantiles <- sort(unique(parse_numbers(values$tp73_quantiles)))
values$block_size <- suppressWarnings(as.numeric(values$block_size))
values$spline_df <- suppressWarnings(as.integer(values$spline_df))
values$minimum_band_count <- suppressWarnings(as.integer(values$minimum_band_count))
values$minimum_class_fraction <- suppressWarnings(as.numeric(
    values$minimum_class_fraction
))
if (!is.finite(values$negative_reference)) {
    stop("--negative-reference must be finite", call. = FALSE)
}
if (any(cofactor_quantiles <= 0 | cofactor_quantiles >= 1) ||
    any(tp73_quantiles <= 0 | tp73_quantiles >= 1)) {
    stop("score quantiles must be in (0,1)", call. = FALSE)
}
if (length(cofactor_quantiles) == 0L || length(tp73_quantiles) < 2L ||
    min(tp73_quantiles) >= max(tp73_quantiles)) {
    stop(
        "provide at least one cofactor quantile and two distinct TP73 quantiles",
        call. = FALSE
    )
}
if (!is.finite(values$block_size) || values$block_size <= 300) {
    stop("--block-size must exceed the 150 bp context radius", call. = FALSE)
}
if (is.na(values$spline_df) || values$spline_df < 2L) {
    stop("--spline-df must be at least 2", call. = FALSE)
}
if (is.na(values$minimum_band_count) || values$minimum_band_count < 1L) {
    stop("--minimum-band-count must be positive", call. = FALSE)
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
if (!nzchar(values$inference_status) || grepl("[\r\n]", values$inference_status)) {
    stop("--inference-status must be non-empty and single-line", call. = FALSE)
}

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
    sql_path <- tempfile("tp73-score-surface-", fileext = ".sql")
    on.exit(unlink(sql_path), add = TRUE)
    writeLines(query, sql_path, useBytes = TRUE)
    command <- paste(
        shQuote(duckdb_path),
        "-light-mode -csv -nullvalue NA :memory: <", shQuote(sql_path)
    )
    fread(cmd = command, na.strings = "NA", showProgress = FALSE)
}

thresholds <- fread(values$thresholds, sep = "\t", na.strings = "NA")
if (!all(c("motif_id", "positive_threshold") %in% names(thresholds))) {
    stop("threshold table requires motif_id and positive_threshold", call. = FALSE)
}
if (!"factor_name" %in% names(thresholds)) {
    thresholds[, factor_name := NA_character_]
}
thresholds <- thresholds[, .(
    motif_id = as.character(motif_id),
    factor_name = as.character(factor_name),
    operating_threshold = as.numeric(positive_threshold)
)]
if (nrow(thresholds) == 0L || any(!nzchar(thresholds$motif_id)) ||
    any(!is.finite(thresholds$operating_threshold)) ||
    anyDuplicated(thresholds$motif_id)) {
    stop("threshold table is empty, duplicated, or non-finite", call. = FALSE)
}
setorder(thresholds, motif_id)

description <- duckdb_fread(paste0(
    "DESCRIBE SELECT * FROM read_parquet(",
    sql_string(values$anchor_evidence), ");"
))
anchor_columns <- as.character(description$column_name)
base_anchor_columns <- c("chrom", "anchor_start", "anchor_end", "anchor_score")
if (length(setdiff(base_anchor_columns, anchor_columns)) > 0L) {
    stop("anchor evidence lacks required coordinate/score columns", call. = FALSE)
}
anti_support_columns <- sort(grep("^supported_anti_", anchor_columns, value = TRUE))
if (length(anti_support_columns) > 0L) {
    sample_ids <- sub("^supported_anti_", "", anti_support_columns)
    control_support_columns <- paste0("supported_control_", sample_ids)
    evidence_column_scheme <- "supported_anti_and_control"
} else {
    anti_support_columns <- sort(grep("^supported_tp73_", anchor_columns, value = TRUE))
    sample_ids <- sub("^supported_tp73_", "", anti_support_columns)
    control_support_columns <- paste0("supported_negative_control_", sample_ids)
    evidence_column_scheme <- "supported_tp73_and_negative_control"
}
if (length(sample_ids) == 0L || anyDuplicated(sample_ids) ||
    length(setdiff(control_support_columns, anchor_columns)) > 0L) {
    stop("anchor evidence lacks matched anti-p73/control support columns", call. = FALSE)
}
samples <- data.table(
    sample_id = sample_ids,
    anti_support = anti_support_columns,
    control_support = control_support_columns
)
excluded_samples <- if (nzchar(values$exclude_samples)) {
    sort(unique(strsplit(values$exclude_samples, ",", fixed = TRUE)[[1L]]))
} else character()
if (any(!nzchar(excluded_samples)) ||
    length(setdiff(excluded_samples, samples$sample_id)) > 0L) {
    stop("--exclude-samples contains an unknown or empty sample ID", call. = FALSE)
}
samples <- samples[!sample_id %in% excluded_samples]
if (nrow(samples) == 0L) {
    stop("--exclude-samples removed every matched evidence sample", call. = FALSE)
}
anchor_select <- c(
    "CAST(chrom AS VARCHAR) AS chrom", "anchor_start", "anchor_end",
    "anchor_score", vapply(
        c(samples$anti_support, samples$control_support),
        sql_identifier, character(1)
    )
)
message("I: reading TP73 anchors and matched CUT&RUN support")
anchors <- duckdb_fread(paste0(
    "SELECT ", paste(anchor_select, collapse = ", "), "\n",
    "FROM read_parquet(", sql_string(values$anchor_evidence), ")\n",
    "ORDER BY chrom, anchor_start, anchor_end;"
))
if (nrow(anchors) == 0L || anchors[, anyDuplicated(
    paste(chrom, anchor_start, anchor_end)
)] != 0L || anchors[, any(!is.finite(anchor_score))]) {
    stop("anchor evidence is empty, duplicated, or has invalid scores", call. = FALSE)
}
for (column in c(samples$anti_support, samples$control_support)) {
    value <- as.logical(anchors[[column]])
    if (anyNA(value)) {
        stop("support column contains missing values: ", column, call. = FALSE)
    }
    set(anchors, j = column, value = value)
}
multi_chromosome <- uniqueN(anchors$chrom) > 1L

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
if (nrow(maxima) == 0L || !setequal(unique(maxima$motif_id), thresholds$motif_id) ||
    nrow(maxima) != nrow(anchors) * nrow(thresholds) ||
    maxima[, anyDuplicated(paste(chrom, anchor_start, anchor_end, motif_id))] != 0L) {
    stop("cofactor maxima are not rectangular, unique, and motif-complete", call. = FALSE)
}
if (any(maxima$context_flank_bp != 150) ||
    !all(maxima$context_distance_metric == "signed_interval_edge_distance")) {
    stop("cofactor maxima do not use the required schema-v4 150 bp geometry",
         call. = FALSE)
}
source_floors <- maxima[, .(
    n_source_floors = uniqueN(source_score_floor),
    source_score_floor = suppressWarnings(min(source_score_floor))
), by = motif_id]
if (source_floors[, any(n_source_floors != 1L)] ||
    source_floors[, any(!is.finite(source_score_floor))]) {
    stop("each motif needs exactly one finite source score floor", call. = FALSE)
}
setkey(source_floors, motif_id)
setkey(anchors, chrom, anchor_start, anchor_end)
setkey(maxima, chrom, anchor_start, anchor_end)

build_tp73_grid <- function(scores) {
    quantile_values <- as.numeric(quantile(
        scores, probs = tp73_quantiles, type = 1, names = FALSE
    ))
    candidates <- rbindlist(list(
        data.table(score = fixed_tp73_grid, source = "fixed"),
        data.table(
            score = quantile_values,
            source = paste0("q", sprintf("%02d", round(100 * tp73_quantiles)))
        )
    ))
    candidates <- candidates[score >= min(scores) & score <= max(scores)]
    candidates[, rounded_key := sprintf("%.9g", score)]
    result <- candidates[, .(
        tp73_score = score[[1L]],
        tp73_grid_source = paste(sort(unique(source)), collapse = "+")
    ), by = rounded_key]
    result[, `:=`(
        tp73_score_percentile = vapply(
            tp73_score, function(value) mean(scores <= value), numeric(1)
        ),
        tp73_grid_order = frank(tp73_score, ties.method = "dense")
    )]
    setorder(result, tp73_score)
    result[, rounded_key := NULL]
    result
}

tp73_grid <- build_tp73_grid(anchors$anchor_score)
tp73_contrast_quantiles <- range(tp73_quantiles)
tp73_contrast_scores <- as.numeric(quantile(
    anchors$anchor_score, probs = tp73_contrast_quantiles,
    type = 1, names = FALSE
))

status_rows <- list()
band_rows <- list()
surface_rows <- list()
tp73_contrast_rows <- list()
status_index <- band_index <- surface_index <- tp73_contrast_index <- 1L

empty_table <- function(columns) {
    as.data.table(setNames(
        rep(list(logical()), length(columns)), columns
    ))
}

band_columns <- c(
    "motif_id", "factor_name", "operating_threshold", "source_score_floor",
    "negative_reference_threshold", "score_band", "score_band_order",
    "score_lower_inclusive", "score_upper_exclusive", "anchors",
    "anchor_fraction", "finite_scores", "minimum_finite_score",
    "median_finite_score", "maximum_finite_score",
    "contains_operating_threshold", "score_band_definition"
)
surface_columns <- c(
    "motif_id", "factor_name", "operating_threshold", "source_score_floor",
    "negative_reference_threshold", "score_band", "score_band_order",
    "score_lower_inclusive", "score_upper_exclusive", "score_band_anchors",
    "contains_operating_threshold", "previous_score_band", "tp73_score",
    "tp73_score_percentile", "tp73_grid_source",
    "adjusted_log_odds_vs_negative", "block_clustered_se_vs_negative",
    "adjusted_odds_ratio_vs_negative",
    "confidence_interval_95_lower_vs_negative",
    "confidence_interval_95_upper_vs_negative", "p_value_vs_negative",
    "adjusted_log_odds_vs_previous_band",
    "block_clustered_se_vs_previous_band",
    "adjusted_odds_ratio_vs_previous_band",
    "confidence_interval_95_lower_vs_previous_band",
    "confidence_interval_95_upper_vs_previous_band",
    "p_value_vs_previous_band", "adjusted_log_odds_vs_weakest_positive",
    "block_clustered_se_vs_weakest_positive",
    "adjusted_odds_ratio_vs_weakest_positive",
    "confidence_interval_95_lower_vs_weakest_positive",
    "confidence_interval_95_upper_vs_weakest_positive",
    "p_value_vs_weakest_positive", "tp73_by_cofactor_interaction_df",
    "tp73_by_cofactor_interaction_chisq",
    "tp73_by_cofactor_interaction_p_value", "discordant_observations",
    "genomic_blocks", "evaluation_status"
)
tp73_contrast_columns <- c(
    "motif_id", "factor_name", "operating_threshold", "source_score_floor",
    "negative_reference_threshold", "score_band", "score_band_order",
    "score_lower_inclusive", "score_upper_exclusive", "score_band_anchors",
    "contains_operating_threshold", "low_tp73_quantile", "low_tp73_score",
    "high_tp73_quantile", "high_tp73_score",
    "adjusted_log_ratio_of_odds_ratios", "block_clustered_standard_error",
    "adjusted_ratio_of_odds_ratios", "confidence_interval_95_lower",
    "confidence_interval_95_upper", "p_value", "interpretation",
    "evaluation_status"
)

model_status <- function(status, note, bands = 0L, evidence = 0L, blocks = 0L,
                         interaction_df = 0L, interaction_statistic = NA_real_,
                         interaction_p_value = NA_real_) {
    data.table(
        evaluation_status = status,
        evaluation_note = note,
        positive_score_bands = bands,
        discordant_observations = evidence,
        genomic_blocks = blocks,
        tp73_by_cofactor_interaction_df = interaction_df,
        tp73_by_cofactor_interaction_chisq = interaction_statistic,
        tp73_by_cofactor_interaction_p_value = interaction_p_value
    )
}

clustered_covariance <- function(fit, clusters) {
    bread <- tryCatch(vcov(fit), error = function(condition) NULL)
    if (is.null(bread) || any(!is.finite(bread))) return(NULL)
    matrix <- fit$x
    score_rows <- matrix * as.numeric(fit$y - fitted(fit))
    scores <- rowsum(score_rows, group = clusters, reorder = FALSE)
    cluster_count <- nrow(scores)
    observations <- nrow(matrix)
    parameters <- ncol(matrix)
    if (cluster_count < 2L || observations <= parameters) return(NULL)
    correction <- (cluster_count / (cluster_count - 1)) *
        ((observations - 1) / (observations - parameters))
    covariance <- bread %*% crossprod(scores) %*% bread * correction
    if (any(!is.finite(covariance))) return(NULL)
    list(covariance = covariance, clusters = cluster_count)
}

wald_test <- function(fit, coefficients, covariance) {
    if (length(coefficients) == 0L) {
        return(list(df = 0L, statistic = NA_real_, p_value = NA_real_))
    }
    selected <- covariance[coefficients, coefficients, drop = FALSE]
    selected <- (selected + t(selected)) / 2
    decomposition <- eigen(selected, symmetric = TRUE)
    tolerance <- max(abs(decomposition$values)) * 1e-8
    retained <- decomposition$values > tolerance
    rank <- sum(retained)
    if (rank == 0L) {
        return(list(df = 0L, statistic = NA_real_, p_value = NA_real_))
    }
    beta <- coef(fit)[coefficients]
    projected <- crossprod(
        decomposition$vectors[, retained, drop = FALSE], beta
    )
    statistic <- sum(
        as.numeric(projected)^2 / decomposition$values[retained]
    )
    list(
        df = rank,
        statistic = statistic,
        p_value = pchisq(statistic, df = rank, lower.tail = FALSE)
    )
}

contrast_result <- function(design_difference, coefficients, covariance) {
    if (all(abs(design_difference) < 1e-12)) {
        return(c(
            log_odds = 0, standard_error = 0, odds_ratio = 1,
            lower = 1, upper = 1, p_value = 1
        ))
    }
    estimate <- as.numeric(design_difference %*% coefficients)
    variance <- as.numeric(design_difference %*% covariance %*% design_difference)
    if (!is.finite(estimate) || !is.finite(variance) || variance <= 0) {
        return(c(
            log_odds = NA_real_, standard_error = NA_real_, odds_ratio = NA_real_,
            lower = NA_real_, upper = NA_real_, p_value = NA_real_
        ))
    }
    standard_error <- sqrt(variance)
    z_value <- estimate / standard_error
    c(
        log_odds = estimate,
        standard_error = standard_error,
        odds_ratio = exp(estimate),
        lower = exp(estimate - 1.96 * standard_error),
        upper = exp(estimate + 1.96 * standard_error),
        p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE)
    )
}

for (motif_index in seq_len(nrow(thresholds))) {
    motif <- thresholds[motif_index]
    current_motif <- motif$motif_id
    motif_maxima <- maxima[motif_id == current_motif]
    source_floor <- source_floors[current_motif, source_score_floor]
    features <- motif_maxima[anchors, on = .(chrom, anchor_start, anchor_end)]
    if (nrow(features) != nrow(anchors)) {
        stop("cofactor join changed anchor cardinality for ", current_motif,
             call. = FALSE)
    }
    score <- as.numeric(features$context_score)
    negative_observable <- source_floor <= values$negative_reference
    positive_scores <- score[!is.na(score) & score >= values$negative_reference]
    negative <- if (negative_observable) {
        is.na(score) | score < values$negative_reference
    } else rep(FALSE, length(score))

    if (!negative_observable) {
        status <- model_status(
            "negative_reference_below_source_floor",
            paste0(
                "negative reference ", values$negative_reference,
                " is below retained scan floor ", source_floor,
                "; absent rows cannot define the requested negative class"
            )
        )
        status_rows[[status_index]] <- cbind(data.table(
            motif_id = current_motif, factor_name = motif$factor_name,
            operating_threshold = motif$operating_threshold,
            source_score_floor = source_floor,
            negative_reference_threshold = values$negative_reference,
            negative_reference_observable = FALSE,
            anchors_total = nrow(anchors), anchors_negative_reference = 0L,
            anchors_with_score_at_or_above_reference = length(positive_scores)
        ), status)
        status_index <- status_index + 1L
        message("I: score surface censored for ", current_motif, " (",
                motif_index, "/", nrow(thresholds), ")")
        next
    }
    if (length(positive_scores) < values$minimum_band_count ||
        sum(negative) < values$minimum_band_count ||
        mean(negative) < values$minimum_class_fraction) {
        status <- model_status(
            "underpowered_score_reference",
            paste0(
                "negative or score-at/above-reference class lacks support; ",
                "minimum negative fraction is ",
                format(values$minimum_class_fraction, trim = TRUE)
            )
        )
        status_rows[[status_index]] <- cbind(data.table(
            motif_id = current_motif, factor_name = motif$factor_name,
            operating_threshold = motif$operating_threshold,
            source_score_floor = source_floor,
            negative_reference_threshold = values$negative_reference,
            negative_reference_observable = TRUE,
            anchors_total = nrow(anchors),
            anchors_negative_reference = sum(negative),
            anchors_with_score_at_or_above_reference = length(positive_scores)
        ), status)
        status_index <- status_index + 1L
        next
    }

    quantile_cuts <- as.numeric(quantile(
        positive_scores, probs = cofactor_quantiles, type = 1, names = FALSE
    ))
    lower_bounds <- sort(unique(c(values$negative_reference, quantile_cuts)))
    lower_bounds <- lower_bounds[lower_bounds >= values$negative_reference]
    upper_bounds <- c(lower_bounds[-1L], Inf)
    positive_band_count <- length(lower_bounds)
    positive_band_ids <- sprintf("positive_%02d", seq_len(positive_band_count))
    band_levels <- c("negative", positive_band_ids)
    assigned_band <- rep(NA_character_, length(score))
    assigned_band[negative] <- "negative"
    for (current_band in seq_len(positive_band_count)) {
        selected <- !is.na(score) & score >= lower_bounds[[current_band]] &
            (is.infinite(upper_bounds[[current_band]]) |
             score < upper_bounds[[current_band]])
        assigned_band[selected] <- positive_band_ids[[current_band]]
    }
    if (anyNA(assigned_band)) {
        stop("score bands failed to partition anchors for ", current_motif,
             call. = FALSE)
    }
    assigned_band <- factor(assigned_band, levels = band_levels)
    band_counts <- tabulate(as.integer(assigned_band), nbins = length(band_levels))

    for (current_band in seq_along(band_levels)) {
        band_id <- band_levels[[current_band]]
        selected <- assigned_band == band_id
        lower <- if (band_id == "negative") -Inf else
            lower_bounds[[current_band - 1L]]
        upper <- if (band_id == "negative") values$negative_reference else
            upper_bounds[[current_band - 1L]]
        finite_selected <- score[selected & !is.na(score)]
        band_rows[[band_index]] <- data.table(
            motif_id = current_motif,
            factor_name = motif$factor_name,
            operating_threshold = motif$operating_threshold,
            source_score_floor = source_floor,
            negative_reference_threshold = values$negative_reference,
            score_band = band_id,
            score_band_order = current_band - 1L,
            score_lower_inclusive = lower,
            score_upper_exclusive = upper,
            anchors = band_counts[[current_band]],
            anchor_fraction = band_counts[[current_band]] / nrow(anchors),
            finite_scores = length(finite_selected),
            minimum_finite_score = if (length(finite_selected)) {
                min(finite_selected)
            } else NA_real_,
            median_finite_score = if (length(finite_selected)) {
                median(finite_selected)
            } else NA_real_,
            maximum_finite_score = if (length(finite_selected)) {
                max(finite_selected)
            } else NA_real_,
            contains_operating_threshold =
                motif$operating_threshold >= lower & motif$operating_threshold < upper,
            score_band_definition = if (band_id == "negative") {
                "context_score < negative reference or observable absent hit"
            } else {
                "empirical quantile band among context_score >= negative reference"
            }
        )
        band_index <- band_index + 1L
    }

    if (any(band_counts < values$minimum_band_count)) {
        status <- model_status(
            "underpowered_score_band",
            paste0("at least one modeled score band has fewer than ",
                   values$minimum_band_count, " anchors"),
            bands = positive_band_count
        )
        status_rows[[status_index]] <- cbind(data.table(
            motif_id = current_motif, factor_name = motif$factor_name,
            operating_threshold = motif$operating_threshold,
            source_score_floor = source_floor,
            negative_reference_threshold = values$negative_reference,
            negative_reference_observable = TRUE,
            anchors_total = nrow(anchors),
            anchors_negative_reference = sum(negative),
            anchors_with_score_at_or_above_reference = length(positive_scores)
        ), status)
        status_index <- status_index + 1L
        next
    }

    evidence_rows <- lapply(seq_len(nrow(samples)), function(sample_index) {
        anti <- anchors[[samples$anti_support[[sample_index]]]]
        control <- anchors[[samples$control_support[[sample_index]]]]
        selected <- anti != control
        data.table(
            chrom = anchors$chrom[selected],
            anchor_start = anchors$anchor_start[selected],
            anchor_score = anchors$anchor_score[selected],
            sample_id = samples$sample_id[[sample_index]],
            score_band = assigned_band[selected],
            outcome = as.integer(anti[selected])
        )
    })
    evidence <- rbindlist(evidence_rows, use.names = TRUE)
    evidence[, `:=`(
        sample_id = factor(sample_id, levels = samples$sample_id),
        score_band = factor(score_band, levels = band_levels),
        genomic_block = paste(
            chrom, floor(anchor_start / values$block_size), sep = ":"
        )
    )]
    if (multi_chromosome) {
        evidence[, chromosome := factor(
            chrom, levels = sort(unique(anchors$chrom))
        )]
    }
    formula <- as.formula(paste0(
        "outcome ~ sample_id",
        if (multi_chromosome) " + chromosome" else "",
        " + splines::ns(anchor_score, df = ", values$spline_df,
        ") * score_band"
    ))
    current_fit <- tryCatch(
        suppressWarnings(glm(
            formula, data = evidence, family = binomial(),
            control = glm.control(maxit = 50L, epsilon = 1e-8),
            x = TRUE, y = TRUE
        )),
        error = function(condition) NULL
    )
    if (is.null(current_fit) || !isTRUE(current_fit$converged) ||
        any(!is.finite(coef(current_fit)))) {
        status <- model_status(
            "not_estimable", "score-surface GLM did not converge",
            bands = positive_band_count, evidence = nrow(evidence),
            blocks = uniqueN(evidence$genomic_block)
        )
        status_rows[[status_index]] <- cbind(data.table(
            motif_id = current_motif, factor_name = motif$factor_name,
            operating_threshold = motif$operating_threshold,
            source_score_floor = source_floor,
            negative_reference_threshold = values$negative_reference,
            negative_reference_observable = TRUE,
            anchors_total = nrow(anchors),
            anchors_negative_reference = sum(negative),
            anchors_with_score_at_or_above_reference = length(positive_scores)
        ), status)
        status_index <- status_index + 1L
        next
    }
    robust <- clustered_covariance(current_fit, evidence$genomic_block)
    if (is.null(robust)) {
        status <- model_status(
            "not_estimable", "block-clustered covariance is unavailable",
            bands = positive_band_count, evidence = nrow(evidence),
            blocks = uniqueN(evidence$genomic_block)
        )
        status_rows[[status_index]] <- cbind(data.table(
            motif_id = current_motif, factor_name = motif$factor_name,
            operating_threshold = motif$operating_threshold,
            source_score_floor = source_floor,
            negative_reference_threshold = values$negative_reference,
            negative_reference_observable = TRUE,
            anchors_total = nrow(anchors),
            anchors_negative_reference = sum(negative),
            anchors_with_score_at_or_above_reference = length(positive_scores)
        ), status)
        status_index <- status_index + 1L
        next
    }
    interaction_coefficients <- grep(
        "splines::ns\\(anchor_score.*:score_band|score_band.*:splines::ns\\(anchor_score",
        names(coef(current_fit)), value = TRUE
    )
    interaction_test <- wald_test(
        current_fit, interaction_coefficients, robust$covariance
    )
    status <- model_status(
        "ok",
        paste(
            "ordered cofactor score bands versus strict negative reference;",
            "discordant anti-p73/control observations; SE clusters",
            format(values$block_size, scientific = FALSE), "bp genomic blocks"
        ),
        bands = positive_band_count, evidence = nrow(evidence),
        blocks = robust$clusters, interaction_df = interaction_test$df,
        interaction_statistic = interaction_test$statistic,
        interaction_p_value = interaction_test$p_value
    )
    status_rows[[status_index]] <- cbind(data.table(
        motif_id = current_motif, factor_name = motif$factor_name,
        operating_threshold = motif$operating_threshold,
        source_score_floor = source_floor,
        negative_reference_threshold = values$negative_reference,
        negative_reference_observable = TRUE,
        anchors_total = nrow(anchors),
        anchors_negative_reference = sum(negative),
        anchors_with_score_at_or_above_reference = length(positive_scores)
    ), status)
    status_index <- status_index + 1L

    reference_sample <- samples$sample_id[[1L]]
    reference_chromosome <- sort(unique(anchors$chrom))[[1L]]
    terms_without_response <- delete.response(terms(current_fit))
    coefficient_values <- coef(current_fit)
    for (current_band in seq_len(positive_band_count)) {
        band_id <- positive_band_ids[[current_band]]
        previous_band <- band_levels[[current_band]]
        metadata <- band_rows[[
            band_index - positive_band_count - 1L + current_band
        ]]
        for (grid_index in seq_len(nrow(tp73_grid))) {
            grid <- tp73_grid[grid_index]
            new_data <- data.frame(
                sample_id = factor(
                    rep(reference_sample, 4L), levels = samples$sample_id
                ),
                anchor_score = rep(grid$tp73_score, 4L),
                score_band = factor(
                    c("negative", previous_band, band_id, positive_band_ids[[1L]]),
                    levels = band_levels
                )
            )
            if (multi_chromosome) {
                new_data$chromosome <- factor(
                    rep(reference_chromosome, 4L),
                    levels = sort(unique(anchors$chrom))
                )
            }
            design <- model.matrix(
                terms_without_response, new_data,
                contrasts.arg = current_fit$contrasts,
                xlev = current_fit$xlevels
            )
            versus_negative <- contrast_result(
                design[3L, ] - design[1L, ],
                coefficient_values, robust$covariance
            )
            versus_previous <- contrast_result(
                design[3L, ] - design[2L, ],
                coefficient_values, robust$covariance
            )
            versus_weakest <- contrast_result(
                design[3L, ] - design[4L, ],
                coefficient_values, robust$covariance
            )
            surface_rows[[surface_index]] <- data.table(
                motif_id = current_motif,
                factor_name = motif$factor_name,
                operating_threshold = motif$operating_threshold,
                source_score_floor = source_floor,
                negative_reference_threshold = values$negative_reference,
                score_band = band_id,
                score_band_order = current_band,
                score_lower_inclusive = metadata$score_lower_inclusive,
                score_upper_exclusive = metadata$score_upper_exclusive,
                score_band_anchors = metadata$anchors,
                contains_operating_threshold =
                    metadata$contains_operating_threshold,
                previous_score_band = previous_band,
                tp73_score = grid$tp73_score,
                tp73_score_percentile = grid$tp73_score_percentile,
                tp73_grid_source = grid$tp73_grid_source,
                adjusted_log_odds_vs_negative = versus_negative[["log_odds"]],
                block_clustered_se_vs_negative =
                    versus_negative[["standard_error"]],
                adjusted_odds_ratio_vs_negative =
                    versus_negative[["odds_ratio"]],
                confidence_interval_95_lower_vs_negative =
                    versus_negative[["lower"]],
                confidence_interval_95_upper_vs_negative =
                    versus_negative[["upper"]],
                p_value_vs_negative = versus_negative[["p_value"]],
                adjusted_log_odds_vs_previous_band =
                    versus_previous[["log_odds"]],
                block_clustered_se_vs_previous_band =
                    versus_previous[["standard_error"]],
                adjusted_odds_ratio_vs_previous_band =
                    versus_previous[["odds_ratio"]],
                confidence_interval_95_lower_vs_previous_band =
                    versus_previous[["lower"]],
                confidence_interval_95_upper_vs_previous_band =
                    versus_previous[["upper"]],
                p_value_vs_previous_band = versus_previous[["p_value"]],
                adjusted_log_odds_vs_weakest_positive =
                    versus_weakest[["log_odds"]],
                block_clustered_se_vs_weakest_positive =
                    versus_weakest[["standard_error"]],
                adjusted_odds_ratio_vs_weakest_positive =
                    versus_weakest[["odds_ratio"]],
                confidence_interval_95_lower_vs_weakest_positive =
                    versus_weakest[["lower"]],
                confidence_interval_95_upper_vs_weakest_positive =
                    versus_weakest[["upper"]],
                p_value_vs_weakest_positive = versus_weakest[["p_value"]],
                tp73_by_cofactor_interaction_df = interaction_test$df,
                tp73_by_cofactor_interaction_chisq = interaction_test$statistic,
                tp73_by_cofactor_interaction_p_value = interaction_test$p_value,
                discordant_observations = nrow(evidence),
                genomic_blocks = robust$clusters,
                evaluation_status = "ok"
            )
            surface_index <- surface_index + 1L
        }

        contrast_data <- data.frame(
            sample_id = factor(
                rep(reference_sample, 4L), levels = samples$sample_id
            ),
            anchor_score = c(
                tp73_contrast_scores[[1L]], tp73_contrast_scores[[1L]],
                tp73_contrast_scores[[2L]], tp73_contrast_scores[[2L]]
            ),
            score_band = factor(
                c("negative", band_id, "negative", band_id),
                levels = band_levels
            )
        )
        if (multi_chromosome) {
            contrast_data$chromosome <- factor(
                rep(reference_chromosome, 4L),
                levels = sort(unique(anchors$chrom))
            )
        }
        contrast_design <- model.matrix(
            terms_without_response, contrast_data,
            contrasts.arg = current_fit$contrasts,
            xlev = current_fit$xlevels
        )
        tp73_increase <- contrast_result(
            (contrast_design[4L, ] - contrast_design[3L, ]) -
                (contrast_design[2L, ] - contrast_design[1L, ]),
            coefficient_values, robust$covariance
        )
        tp73_contrast_rows[[tp73_contrast_index]] <- data.table(
            motif_id = current_motif,
            factor_name = motif$factor_name,
            operating_threshold = motif$operating_threshold,
            source_score_floor = source_floor,
            negative_reference_threshold = values$negative_reference,
            score_band = band_id,
            score_band_order = current_band,
            score_lower_inclusive = metadata$score_lower_inclusive,
            score_upper_exclusive = metadata$score_upper_exclusive,
            score_band_anchors = metadata$anchors,
            contains_operating_threshold = metadata$contains_operating_threshold,
            low_tp73_quantile = tp73_contrast_quantiles[[1L]],
            low_tp73_score = tp73_contrast_scores[[1L]],
            high_tp73_quantile = tp73_contrast_quantiles[[2L]],
            high_tp73_score = tp73_contrast_scores[[2L]],
            adjusted_log_ratio_of_odds_ratios = tp73_increase[["log_odds"]],
            block_clustered_standard_error =
                tp73_increase[["standard_error"]],
            adjusted_ratio_of_odds_ratios = tp73_increase[["odds_ratio"]],
            confidence_interval_95_lower = tp73_increase[["lower"]],
            confidence_interval_95_upper = tp73_increase[["upper"]],
            p_value = tp73_increase[["p_value"]],
            interpretation = paste(
                "cofactor-vs-negative occupancy OR at high TP73 score divided",
                "by the same OR at low TP73 score"
            ),
            evaluation_status = "ok"
        )
        tp73_contrast_index <- tp73_contrast_index + 1L
    }
    message("I: fitted score surface for ", current_motif, " (",
            motif_index, "/", nrow(thresholds), ")")
}

status_table <- rbindlist(status_rows, use.names = TRUE, fill = TRUE)
band_table <- if (length(band_rows)) {
    rbindlist(band_rows, use.names = TRUE, fill = TRUE)
} else empty_table(band_columns)
surface_table <- if (length(surface_rows)) {
    rbindlist(surface_rows, use.names = TRUE, fill = TRUE)
} else empty_table(surface_columns)
tp73_contrast_table <- if (length(tp73_contrast_rows)) {
    rbindlist(tp73_contrast_rows, use.names = TRUE, fill = TRUE)
} else empty_table(tp73_contrast_columns)
setcolorder(band_table, band_columns)
setcolorder(surface_table, surface_columns)
setcolorder(tp73_contrast_table, tp73_contrast_columns)
setorder(status_table, motif_id)
if (nrow(band_table)) setorder(band_table, motif_id, score_band_order)
if (nrow(surface_table)) {
    setorder(surface_table, motif_id, score_band_order, tp73_score)
}
if (nrow(tp73_contrast_table)) {
    setorder(tp73_contrast_table, motif_id, score_band_order)
}
fwrite(status_table, paste0(output_prefix, "_motif_status.tsv"), sep = "\t")
fwrite(band_table, paste0(output_prefix, "_score_band.tsv"), sep = "\t")
fwrite(surface_table, paste0(output_prefix, "_score_surface.tsv"), sep = "\t")
fwrite(
    tp73_contrast_table,
    paste0(output_prefix, "_tp73_score_contrast.tsv"), sep = "\t"
)

run_config <- data.table(
    anchor_evidence = values$anchor_evidence,
    cofactor_maxima = values$cofactor_maxima,
    thresholds = values$thresholds,
    motifs = paste(thresholds$motif_id, collapse = ","),
    evidence_column_scheme = evidence_column_scheme,
    sample_ids = paste(samples$sample_id, collapse = ","),
    excluded_sample_ids = paste(excluded_samples, collapse = ","),
    negative_reference_threshold = values$negative_reference,
    negative_reference_semantics = "strict context_score < N or absent",
    score_band_semantics = paste(
        "negative reference plus empirical non-overlapping quantile bands among",
        "finite context_score >= N; positive bands are not cumulative thresholds"
    ),
    cofactor_quantiles = paste(cofactor_quantiles, collapse = ","),
    tp73_fixed_grid = paste(fixed_tp73_grid, collapse = ","),
    tp73_grid_quantiles = paste(tp73_quantiles, collapse = ","),
    model_formula = paste0(
        "outcome ~ sample_id",
        if (multi_chromosome) " + chromosome" else "",
        " + ns(anchor_score,df=", values$spline_df, ") * score_band"
    ),
    outcome = "strict_immersion anti-p73 versus matched negative control",
    discordance_conditioning = "retain anti-p73/control-discordant anchor-samples",
    context_flank_bp = 150,
    context_geometry = "signed_interval_edge_distance",
    block_size_bp = values$block_size,
    spline_df = values$spline_df,
    minimum_band_count = values$minimum_band_count,
    minimum_class_fraction = values$minimum_class_fraction,
    tp73_contrast_quantiles = paste(tp73_contrast_quantiles, collapse = ","),
    tp73_contrast_scores = paste(tp73_contrast_scores, collapse = ","),
    source_commit = values$source_commit,
    source_dirty = values$source_dirty,
    chromosomes = paste(sort(unique(anchors$chrom)), collapse = ","),
    inference_status = values$inference_status
)
fwrite(run_config, paste0(output_prefix, "_run_config.tsv"), sep = "\t")
message("I: wrote TP73/cofactor score-response outputs with prefix ", output_prefix)
