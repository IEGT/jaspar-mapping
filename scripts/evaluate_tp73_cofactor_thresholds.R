#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: evaluate_tp73_cofactor_thresholds.R --anchor-evidence FILE",
        "       --cofactor-maxima FILE --output-prefix PATH [options]",
        "",
        "Assess score thresholds for interval-defined cofactor maxima around TP73",
        "anchors. The anchor table supplies TP73 scores and matched anti-p73/control",
        "CUT&RUN labels. Five contiguous chromosome folds keep all sample copies of",
        "an anchor together. By default, each threshold is compared with one fixed",
        "negative reference: no maximum reaches -1. Scores at or above -1 but",
        "below the positive threshold are excluded rather than called negative.",
        "",
        "Options:",
        "  --anchor-evidence FILE  TP73 anchor and CUT&RUN feature Parquet (required)",
        "  --cofactor-maxima FILE  Long v4 interval-maximum Parquet (required)",
        "  --output-prefix PATH    Basename for tables and plot (required)",
        "  --thresholds LIST       Comma-separated score floors; 'auto' for every",
        "                          observed integer floor from 0 upward; or",
        "                          'auto-source-floor' for every integer from the",
        "                          recorded source floor through 0",
        "                          (default: 0,2,4,5,6,8,10,12)",
        "  --folds N               Contiguous genomic folds (default: 5)",
        "  --chrom-size BP         Chromosome span; default derives from anchors",
        "  --spline-df N           TP73 natural-spline degrees of freedom (default: 4)",
        "  --minimum-class-fraction N",
        "                          Skip fitting when either anchor class is smaller",
        "                          than N (default: 0)",
        "  --comparison-mode MODE  'fixed-negative-reference' (default) or the",
        "                          historical 'threshold-complement'",
        "  --negative-reference-threshold N",
        "                          Negative when the maximum is below this score",
        "                          or absent; the boundary is intermediate",
        "                          (default: -1)",
        "  --duckdb COMMAND        DuckDB CLI executable (default: duckdb)",
        "  --compact-output        Omit detailed per-fold tables and the plot",
        "  -h, --help              Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) {
    usage()
}

values <- list(
    thresholds = "0,2,4,5,6,8,10,12",
    folds = 5L,
    chrom_size = "auto",
    spline_df = 4L,
    minimum_class_fraction = 0,
    comparison_mode = "fixed-negative-reference",
    negative_reference_threshold = -1,
    compact_output = FALSE,
    duckdb = "duckdb"
)
known <- c(
    "--anchor-evidence", "--cofactor-maxima", "--output-prefix",
    "--thresholds", "--folds", "--chrom-size", "--spline-df",
    "--minimum-class-fraction", "--comparison-mode",
    "--negative-reference-threshold", "--duckdb"
)
index <- 1L
while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (option == "--compact-output") {
        values$compact_output <- TRUE
        index <- index + 1L
        next
    }
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

required <- c("anchor_evidence", "cofactor_maxima", "output_prefix")
missing <- required[vapply(required, function(name) is.null(values[[name]]),
                           logical(1))]
if (length(missing) > 0L) {
    writeLines(paste(
        "E: missing required options:",
        paste(paste0("--", gsub("_", "-", missing)), collapse = ", ")
    ), con = stderr())
    usage(2L)
}

for (name in c("folds", "spline_df")) {
    parsed <- suppressWarnings(as.integer(values[[name]]))
    if (is.na(parsed) || parsed <= 0L) {
        stop("--", gsub("_", "-", name), " must be a positive integer",
             call. = FALSE)
    }
    values[[name]] <- parsed
}
if (values$folds < 2L) {
    stop("--folds must be at least 2", call. = FALSE)
}
if (values$spline_df < 2L) {
    stop("--spline-df must be at least 2", call. = FALSE)
}
values$minimum_class_fraction <- suppressWarnings(as.numeric(
    values$minimum_class_fraction
))
if (!is.finite(values$minimum_class_fraction) ||
    values$minimum_class_fraction < 0 ||
    values$minimum_class_fraction >= 0.5) {
    stop("--minimum-class-fraction must be in [0, 0.5)", call. = FALSE)
}
if (!values$comparison_mode %in% c(
    "fixed-negative-reference", "threshold-complement"
)) {
    stop(
        "--comparison-mode must be 'fixed-negative-reference' or ",
        "'threshold-complement'", call. = FALSE
    )
}
values$negative_reference_threshold <- suppressWarnings(as.numeric(
    values$negative_reference_threshold
))
if (!is.finite(values$negative_reference_threshold)) {
    stop("--negative-reference-threshold must be finite", call. = FALSE)
}
threshold_specification <- tolower(values$thresholds)
threshold_mode <- switch(
    threshold_specification,
    "auto" = "observed_integer_grid",
    "auto-source-floor" = "source_floor_integer_grid",
    "explicit"
)
if (threshold_mode == "explicit") {
    thresholds <- suppressWarnings(as.numeric(strsplit(
        values$thresholds, ",", fixed = TRUE
    )[[1L]]))
    if (length(thresholds) == 0L || any(!is.finite(thresholds))) {
        stop(
            "--thresholds must be 'auto', 'auto-source-floor', or a non-empty ",
            "list of finite numbers",
            call. = FALSE
        )
    }
    thresholds <- sort(unique(thresholds))
} else {
    thresholds <- numeric()
}

suppressPackageStartupMessages({
    library(data.table)
    if (!values$compact_output) {
        library(ggplot2)
    }
})

values$anchor_evidence <- normalizePath(values$anchor_evidence, mustWork = TRUE)
values$cofactor_maxima <- normalizePath(values$cofactor_maxima, mustWork = TRUE)
duckdb_path <- Sys.which(values$duckdb)
if (!nzchar(duckdb_path)) {
    stop("DuckDB CLI not found: ", values$duckdb, call. = FALSE)
}
output_prefix <- values$output_prefix
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

sql_string <- function(value) {
    paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
}

duckdb_fread <- function(query) {
    sql_path <- tempfile("tp73-cofactor-thresholds-", fileext = ".sql")
    on.exit(unlink(sql_path), add = TRUE)
    writeLines(query, sql_path, useBytes = TRUE)
    command <- paste(
        shQuote(duckdb_path),
        "-light-mode -csv -nullvalue NA :memory: <", shQuote(sql_path)
    )
    fread(cmd = command, na.strings = "NA", showProgress = FALSE)
}

anchor_schema <- duckdb_fread(paste0(
    "DESCRIBE SELECT * FROM read_parquet(",
    sql_string(values$anchor_evidence), ");"
))
anchor_source_columns <- as.character(anchor_schema$column_name)
anti_columns <- sort(grep(
    "^supported_anti_", anchor_source_columns, value = TRUE
))
if (length(anti_columns) > 0L) {
    sample_ids <- sub("^supported_anti_", "", anti_columns)
    control_columns <- paste0("supported_control_", sample_ids)
    evidence_column_scheme <- "supported_anti_and_control"
} else {
    anti_columns <- sort(grep(
        "^supported_tp73_", anchor_source_columns, value = TRUE
    ))
    sample_ids <- sub("^supported_tp73_", "", anti_columns)
    control_columns <- paste0("supported_negative_control_", sample_ids)
    evidence_column_scheme <- "supported_tp73_and_negative_control"
}
if (length(sample_ids) == 0L ||
    any(!control_columns %in% anchor_source_columns)) {
    stop(
        "anchor evidence needs complete matched anti/control or ",
        "TP73/negative-control support columns", call. = FALSE
    )
}
samples <- data.table(
    sample_id = sample_ids,
    anti = anti_columns,
    control = control_columns
)
support_columns <- c(samples$anti, samples$control)
anchor_query <- paste0(
    "SELECT CAST(chrom AS VARCHAR) AS chrom, anchor_start, anchor_end, ",
    "anchor_score, ", paste(support_columns, collapse = ", "), "\n",
    "FROM read_parquet(", sql_string(values$anchor_evidence), ")\n",
    "ORDER BY chrom, anchor_start;"
)
message("I: reading TP73 anchors and CUT&RUN labels")
anchors <- duckdb_fread(anchor_query)
expected_anchor_columns <- c(
    "chrom", "anchor_start", "anchor_end", "anchor_score", support_columns
)
missing_anchor_columns <- setdiff(expected_anchor_columns, names(anchors))
if (length(missing_anchor_columns) > 0L) {
    stop(
        "anchor evidence lacks columns: ",
        paste(missing_anchor_columns, collapse = ", "), call. = FALSE
    )
}
if (nrow(anchors) == 0L || anchors[, uniqueN(chrom)] != 1L) {
    stop("anchor evidence must contain one non-empty chromosome", call. = FALSE)
}
if (anchors[, anyDuplicated(paste(chrom, anchor_start, anchor_end))] != 0L) {
    stop("anchor evidence contains duplicate genomic spans", call. = FALSE)
}
if (anchors[, any(!is.finite(anchor_score))]) {
    stop("anchor evidence contains non-finite TP73 scores", call. = FALSE)
}

maxima_query <- paste0(
    "SELECT CAST(chrom AS VARCHAR) AS chrom, anchor_start, anchor_end, ",
    "motif_id, context_score, source_score_floor, context_flank_bp, ",
    "capture_prefilter_center_bp, context_distance_metric\n",
    "FROM read_parquet(", sql_string(values$cofactor_maxima), ")\n",
    "ORDER BY chrom, anchor_start, motif_id;"
)
message("I: reading interval-defined cofactor maxima")
maxima <- duckdb_fread(maxima_query)
expected_maxima_columns <- c(
    "chrom", "anchor_start", "anchor_end", "motif_id", "context_score",
    "source_score_floor", "context_flank_bp",
    "capture_prefilter_center_bp", "context_distance_metric"
)
missing_maxima_columns <- setdiff(expected_maxima_columns, names(maxima))
if (length(missing_maxima_columns) > 0L) {
    stop(
        "cofactor maxima lack columns: ",
        paste(missing_maxima_columns, collapse = ", "), call. = FALSE
    )
}
if (nrow(maxima) == 0L || maxima[, anyDuplicated(
    paste(chrom, anchor_start, anchor_end, motif_id)
)] != 0L) {
    stop("cofactor maxima are empty or have duplicate keys", call. = FALSE)
}
motifs <- sort(unique(maxima$motif_id))
if (length(motifs) == 0L) {
    stop("cofactor maxima contain no motif IDs", call. = FALSE)
}
if (nrow(maxima) != nrow(anchors) * length(motifs)) {
    stop("cofactor maxima are not rectangular by anchor and motif", call. = FALSE)
}
if (!all(unique(maxima$context_distance_metric) ==
         "signed_interval_edge_distance")) {
    stop("cofactor maxima do not use schema-v4 interval geometry", call. = FALSE)
}
source_floors <- unique(maxima$source_score_floor)
if (length(source_floors) != 1L || !is.finite(source_floors[[1L]])) {
    stop("cofactor maxima must have one finite sparse source floor", call. = FALSE)
}
if (threshold_mode == "explicit" && any(thresholds < source_floors[[1L]])) {
    stop("thresholds cannot be below the sparse source score floor", call. = FALSE)
}
if (values$comparison_mode == "fixed-negative-reference") {
    if (values$negative_reference_threshold < source_floors[[1L]]) {
        stop(
            "the fixed negative-reference threshold cannot be below the sparse ",
            "source score floor", call. = FALSE
        )
    }
    if (threshold_mode == "source_floor_integer_grid") {
        stop(
            "auto-source-floor is only compatible with ",
            "--comparison-mode threshold-complement", call. = FALSE
        )
    }
    if (threshold_mode == "explicit" &&
        any(thresholds < values$negative_reference_threshold)) {
        stop(
            "fixed-reference positive thresholds must be at least ",
            "--negative-reference-threshold", call. = FALSE
        )
    }
}

motif_map <- data.table(
    motif_id = motifs,
    score_column = paste0("cofactor_score_", seq_along(motifs))
)
maxima <- merge(maxima, motif_map, by = "motif_id", sort = FALSE)
maxima_wide <- dcast(
    maxima,
    chrom + anchor_start + anchor_end ~ score_column,
    value.var = "context_score"
)
anchor_features <- merge(
    anchors, maxima_wide,
    by = c("chrom", "anchor_start", "anchor_end"),
    all.x = TRUE, sort = FALSE
)
if (nrow(anchor_features) != nrow(anchors)) {
    stop("cofactor maxima changed the anchor cardinality", call. = FALSE)
}

evidence_rows <- lapply(seq_len(nrow(samples)), function(sample_index) {
    anti <- samples$anti[[sample_index]]
    control <- samples$control[[sample_index]]
    selected <- anchor_features[[anti]] != anchor_features[[control]]
    if (anyNA(selected)) {
        stop("CUT&RUN support columns contain missing values", call. = FALSE)
    }
    result <- copy(anchor_features[selected])
    result[, `:=`(
        sample_id = samples$sample_id[[sample_index]],
        outcome = as.integer(get(anti))
    )]
    result[, c(support_columns) := NULL]
    result
})
evidence <- rbindlist(evidence_rows, use.names = TRUE)
if (nrow(evidence) == 0L || evidence[, any(!outcome %in% 0:1)]) {
    stop("discordant CUT&RUN evidence is empty or invalid", call. = FALSE)
}
evidence[, sample_id := factor(sample_id, levels = samples$sample_id)]

if (tolower(values$chrom_size) == "auto") {
    chrom_size <- max(anchor_features$anchor_end)
    chrom_size_source <- "maximum_anchor_end"
} else {
    chrom_size <- suppressWarnings(as.numeric(values$chrom_size))
    if (!is.finite(chrom_size) || chrom_size <= 0) {
        stop("--chrom-size must be positive or 'auto'", call. = FALSE)
    }
    chrom_size_source <- "command_line"
}
if (max(anchor_features$anchor_end) > chrom_size) {
    stop("--chrom-size is shorter than an anchor interval", call. = FALSE)
}
fold_width <- chrom_size / values$folds
evidence[, fold := pmin(
    floor(anchor_start / fold_width) + 1L,
    values$folds
)]
evidence[, fold := as.integer(fold)]
anchor_features[, fold := pmin(
    floor(anchor_start / fold_width) + 1L,
    values$folds
)]
anchor_features[, fold := as.integer(fold)]

fold_manifest <- data.table(
    fold = seq_len(values$folds),
    start = floor((seq_len(values$folds) - 1L) * fold_width),
    end = ceiling(seq_len(values$folds) * fold_width)
)
fold_manifest[fold == values$folds, end := as.numeric(chrom_size)]
fold_support <- evidence[, .(
    observations = .N,
    anchors = uniqueN(paste(chrom, anchor_start, anchor_end)),
    anti_only = sum(outcome == 1L),
    control_only = sum(outcome == 0L),
    samples = uniqueN(sample_id)
), by = fold]
fold_manifest <- merge(fold_manifest, fold_support, by = "fold", all.x = TRUE)
if (fold_manifest[, any(
    is.na(observations) | anti_only == 0L | control_only == 0L |
        samples != nrow(samples)
)]) {
    stop("one or more contiguous folds lack outcomes or samples", call. = FALSE)
}
fwrite(
    fold_manifest, paste0(output_prefix, "_fold_manifest.tsv"), sep = "\t"
)

roc_auc <- function(outcome, prediction) {
    positives <- as.double(sum(outcome == 1L))
    negatives <- as.double(sum(outcome == 0L))
    if (positives == 0L || negatives == 0L) return(NA_real_)
    ranks <- rank(prediction, ties.method = "average")
    (sum(ranks[outcome == 1L]) - positives * (positives + 1) / 2) /
        (positives * negatives)
}

average_precision <- function(outcome, prediction) {
    positives <- sum(outcome == 1L)
    if (positives == 0L) return(NA_real_)
    grouped <- data.table(outcome = outcome, prediction = prediction)[
        , .(positive = sum(outcome == 1L), total = .N), by = prediction
    ][order(-prediction)]
    grouped[, `:=`(
        cumulative_positive = cumsum(positive),
        cumulative_total = cumsum(total)
    )]
    sum((grouped$positive / positives) *
        (grouped$cumulative_positive / grouped$cumulative_total))
}

metric_row <- function(outcome, prediction) {
    epsilon <- 1e-15
    bounded <- pmin(pmax(prediction, epsilon), 1 - epsilon)
    data.table(
        n = length(outcome),
        anti_only_fraction = mean(outcome == 1L),
        roc_auc = roc_auc(outcome, prediction),
        average_precision = average_precision(outcome, prediction),
        log_loss = -mean(
            outcome * log(bounded) + (1 - outcome) * log(1 - bounded)
        ),
        brier_score = mean((outcome - prediction)^2)
    )
}

cell_metrics <- function(data, prediction) {
    rbindlist(lapply(seq_len(values$folds), function(held_out) {
        rbindlist(lapply(samples$sample_id, function(sample) {
            selected <- data$fold == held_out & data$sample_id == sample
            metric_row(data$outcome[selected], prediction[selected])[
                , `:=`(fold = held_out, sample_id = sample)
            ]
        }))
    }))
}

tp73_term <- sprintf(
    "splines::ns(anchor_score, df = %d)", values$spline_df
)
baseline_formula <- as.formula(paste(
    "outcome ~ sample_id +", tp73_term
))
fit_contiguous_folds <- function(data, formula, coefficient_name = NULL) {
    prediction <- rep(NA_real_, nrow(data))
    coefficient_rows <- if (is.null(coefficient_name)) NULL else
        vector("list", values$folds)
    for (held_out in seq_len(values$folds)) {
        training <- data[fold != held_out]
        testing_indexes <- which(data$fold == held_out)
        if (nrow(training) == 0L || length(testing_indexes) == 0L) {
            return(list(error = paste0(
                "fold ", held_out, " has no training or held-out observations"
            )))
        }
        if (!is.null(coefficient_name) &&
            training[, uniqueN(get(coefficient_name))] != 2L) {
            return(list(error = paste0(
                "training fold ", held_out, " has a constant indicator"
            )))
        }
        fit <- tryCatch(
            suppressWarnings(glm(
                formula, data = training, family = binomial(),
                control = glm.control(maxit = 50L, epsilon = 1e-8)
            )),
            error = function(condition) NULL
        )
        coefficient <- if (is.null(fit) || is.null(coefficient_name)) {
            NA_real_
        } else {
            unname(coef(fit)[coefficient_name])
        }
        if (is.null(fit) || !isTRUE(fit$converged) ||
            (!is.null(coefficient_name) && !is.finite(coefficient))) {
            return(list(error = paste0(
                "model did not converge in fold ", held_out
            )))
        }
        prediction[testing_indexes] <- tryCatch(
            predict(fit, newdata = data[testing_indexes], type = "response"),
            error = function(condition) rep(NA_real_, length(testing_indexes))
        )
        if (!is.null(coefficient_name)) {
            coefficient_rows[[held_out]] <- data.table(
                fold = held_out,
                retained_log_odds = coefficient,
                adjusted_odds_ratio = exp(coefficient)
            )
        }
        rm(fit, training)
        gc(verbose = FALSE)
    }
    if (any(!is.finite(prediction))) {
        return(list(error = "prediction vector is incomplete"))
    }
    list(
        prediction = prediction,
        coefficients = if (is.null(coefficient_name)) NULL else
            rbindlist(coefficient_rows),
        error = NULL
    )
}

message("I: fitting TP73-only contiguous-fold baseline")
baseline_fit <- fit_contiguous_folds(evidence, baseline_formula)
if (!is.null(baseline_fit$error)) {
    stop("TP73-only ", baseline_fit$error, call. = FALSE)
}
baseline_prediction <- baseline_fit$prediction
baseline_cells <- cell_metrics(evidence, baseline_prediction)
baseline_macro <- baseline_cells[, lapply(.SD, mean), .SDcols = c(
    "roc_auc", "average_precision", "log_loss", "brier_score"
)]
baseline_metrics <- rbindlist(list(
    metric_row(evidence$outcome, baseline_prediction)[, scope := "pooled"],
    baseline_macro[, `:=`(
        n = nrow(evidence),
        anti_only_fraction = mean(evidence$outcome == 1L),
        scope = "sample_fold_macro_mean"
    )]
), use.names = TRUE, fill = TRUE)
setcolorder(
    baseline_metrics,
    c("scope", "n", "anti_only_fraction", "roc_auc", "average_precision",
      "log_loss", "brier_score")
)
fwrite(
    baseline_metrics, paste0(output_prefix, "_baseline_metrics.tsv"), sep = "\t"
)

threshold_rows <- list()
sample_fold_rows <- list()
coefficient_rows <- list()
raw_rows <- list()
unevaluable_threshold_row <- function(motif_id, threshold, counts,
                                      discordant_observations, status, note,
                                      candidate_baseline = NULL) {
    baseline <- if (is.null(candidate_baseline)) {
        list(
            roc_auc = NA_real_, average_precision = NA_real_,
            log_loss = NA_real_, brier_score = NA_real_
        )
    } else {
        candidate_baseline
    }
    data.table(
        motif_id = motif_id,
        threshold = threshold,
        anchors_total = nrow(anchor_features),
        anchors_retained = counts$positive,
        retained_fraction = counts$positive / nrow(anchor_features),
        anchors_negative_reference = counts$negative,
        negative_reference_fraction = counts$negative / nrow(anchor_features),
        anchors_intermediate = counts$intermediate,
        intermediate_fraction = counts$intermediate / nrow(anchor_features),
        anchors_compared = counts$positive + counts$negative,
        comparison_fraction =
            (counts$positive + counts$negative) / nrow(anchor_features),
        discordant_observations_total = nrow(evidence),
        discordant_observations = discordant_observations,
        comparison_mode = values$comparison_mode,
        negative_reference_threshold = if (
            values$comparison_mode == "fixed-negative-reference"
        ) values$negative_reference_threshold else NA_real_,
        negative_reference_inclusive = if (
            values$comparison_mode == "fixed-negative-reference"
        ) FALSE else NA,
        intermediate_excluded =
            values$comparison_mode == "fixed-negative-reference",
        baseline_macro_roc_auc = baseline$roc_auc,
        augmented_macro_roc_auc = NA_real_,
        delta_macro_roc_auc = NA_real_,
        baseline_macro_average_precision = baseline$average_precision,
        augmented_macro_average_precision = NA_real_,
        delta_macro_average_precision = NA_real_,
        baseline_macro_log_loss = baseline$log_loss,
        augmented_macro_log_loss = NA_real_,
        delta_macro_log_loss = NA_real_,
        baseline_macro_brier_score = baseline$brier_score,
        augmented_macro_brier_score = NA_real_,
        delta_macro_brier_score = NA_real_,
        median_adjusted_odds_ratio = NA_real_,
        minimum_adjusted_odds_ratio = NA_real_,
        maximum_adjusted_odds_ratio = NA_real_,
        median_raw_sample_odds_ratio = NA_real_,
        samples_with_raw_odds_ratio_below_one = NA_integer_,
        samples_total = nrow(samples),
        sample_fold_cells = values$folds * nrow(samples),
        sample_fold_cells_with_roc_auc_gain = NA_integer_,
        sample_fold_cells_with_log_loss_gain = NA_integer_,
        evaluation_status = status,
        evaluation_note = note
    )
}
for (motif_index in seq_len(nrow(motif_map))) {
    motif_id <- motif_map$motif_id[[motif_index]]
    score_column <- motif_map$score_column[[motif_index]]
    anchor_scores <- anchor_features[[score_column]]
    evidence_scores <- evidence[[score_column]]
    motif_thresholds <- thresholds
    if (threshold_mode != "explicit") {
        observed_maximum <- suppressWarnings(max(anchor_scores, na.rm = TRUE))
        grid_minimum <- if (threshold_mode == "source_floor_integer_grid") {
            ceiling(source_floors[[1L]])
        } else if (values$comparison_mode == "fixed-negative-reference") {
            max(0, ceiling(values$negative_reference_threshold))
        } else {
            0
        }
        if (threshold_mode == "source_floor_integer_grid") {
            if (grid_minimum > 0) {
                stop(
                    "auto-source-floor requires a source score floor at or below 0",
                    call. = FALSE
                )
            }
            motif_thresholds <- seq.int(grid_minimum, 0)
        } else if (!is.finite(observed_maximum) || observed_maximum < grid_minimum) {
            message("I: no observed threshold at or above ", grid_minimum,
                    " for ", motif_id)
            motif_thresholds <- grid_minimum
        } else {
            motif_thresholds <- seq.int(grid_minimum, floor(observed_maximum))
        }
    }
    for (threshold in motif_thresholds) {
        message("I: fitting ", motif_id, " at score >= ", threshold)
        anchor_positive <- !is.na(anchor_scores) & anchor_scores >= threshold
        evidence_positive <- !is.na(evidence_scores) &
            evidence_scores >= threshold
        if (values$comparison_mode == "fixed-negative-reference") {
            anchor_negative <- is.na(anchor_scores) |
                anchor_scores < values$negative_reference_threshold
            evidence_negative <- is.na(evidence_scores) |
                evidence_scores < values$negative_reference_threshold
        } else {
            anchor_negative <- !anchor_positive
            evidence_negative <- !evidence_positive
        }
        anchor_intermediate <- !anchor_positive & !anchor_negative
        counts <- list(
            positive = sum(anchor_positive),
            negative = sum(anchor_negative),
            intermediate = sum(anchor_intermediate)
        )
        evidence[, retained := fifelse(
            evidence_positive, 1L,
            fifelse(evidence_negative, 0L, NA_integer_)
        )]
        candidate_evidence <- evidence[!is.na(retained)]
        if (counts$positive == 0L || counts$negative == 0L) {
            if (threshold_mode != "explicit") {
                message(
                    "I: skipping ", motif_id, " at ", threshold,
                    " because the positive or comparison class is empty"
                )
                threshold_rows[[length(threshold_rows) + 1L]] <-
                    unevaluable_threshold_row(
                        motif_id, threshold, counts, nrow(candidate_evidence),
                        "empty_comparison_class",
                        "positive or comparison anchor class is empty"
                    )
                next
            }
            stop(
                motif_id, " threshold ", threshold,
                " has an empty positive or comparison anchor class",
                call. = FALSE
            )
        }
        positive_fraction <- counts$positive / nrow(anchor_features)
        negative_fraction <- counts$negative / nrow(anchor_features)
        if (positive_fraction < values$minimum_class_fraction ||
            negative_fraction < values$minimum_class_fraction) {
            message(
                "I: skipping ", motif_id, " at ", threshold,
                " because positive/comparison support is below ",
                values$minimum_class_fraction
            )
            threshold_rows[[length(threshold_rows) + 1L]] <-
                unevaluable_threshold_row(
                    motif_id, threshold, counts, nrow(candidate_evidence),
                    "insufficient_class_support",
                    paste0(
                        "positive or comparison anchor fraction is below ",
                        "the fit floor"
                    )
                )
            next
        }
        cell_support <- candidate_evidence[, .(
            outcomes = uniqueN(outcome)
        ), by = .(fold, sample_id)]
        if (nrow(cell_support) != values$folds * nrow(samples) ||
            any(cell_support$outcomes != 2L)) {
            note <- "one or more held-out sample/fold cells lack both outcomes"
            if (threshold_mode == "explicit") {
                stop(motif_id, " threshold ", threshold, " ", note,
                     call. = FALSE)
            }
            threshold_rows[[length(threshold_rows) + 1L]] <-
                unevaluable_threshold_row(
                    motif_id, threshold, counts, nrow(candidate_evidence),
                    "model_not_evaluable", note
                )
            next
        }

        if (values$comparison_mode == "fixed-negative-reference") {
            candidate_baseline_fit <- fit_contiguous_folds(
                candidate_evidence, baseline_formula
            )
            if (!is.null(candidate_baseline_fit$error)) {
                threshold_rows[[length(threshold_rows) + 1L]] <-
                    unevaluable_threshold_row(
                        motif_id, threshold, counts, nrow(candidate_evidence),
                        "model_not_evaluable",
                        paste("TP73 baseline", candidate_baseline_fit$error)
                    )
                next
            }
            candidate_baseline_prediction <- candidate_baseline_fit$prediction
            candidate_baseline_cells <- cell_metrics(
                candidate_evidence, candidate_baseline_prediction
            )
            candidate_baseline_macro <- candidate_baseline_cells[
                , lapply(.SD, mean), .SDcols = c(
                    "roc_auc", "average_precision", "log_loss", "brier_score"
                )
            ]
        } else {
            candidate_baseline_prediction <- baseline_prediction
            candidate_baseline_cells <- baseline_cells
            candidate_baseline_macro <- baseline_macro
        }
        augmented_formula <- as.formula(paste(
            "outcome ~ sample_id +", tp73_term, "+ retained"
        ))
        augmented_fit <- fit_contiguous_folds(
            candidate_evidence, augmented_formula, "retained"
        )
        if (!is.null(augmented_fit$error)) {
            message(
                "I: skipping ", motif_id, " at ", threshold,
                ": ", augmented_fit$error
            )
            if (threshold_mode == "explicit") {
                stop(
                    motif_id, " threshold ", threshold, " ",
                    augmented_fit$error, call. = FALSE
                )
            }
            threshold_rows[[length(threshold_rows) + 1L]] <-
                unevaluable_threshold_row(
                    motif_id, threshold, counts, nrow(candidate_evidence),
                    "model_not_evaluable", augmented_fit$error,
                    candidate_baseline_macro
                )
            next
        }
        prediction <- augmented_fit$prediction
        coefficients <- augmented_fit$coefficients
        coefficients[, `:=`(
            motif_id = motif_id,
            threshold = threshold,
            comparison_mode = values$comparison_mode
        )]
        if (!values$compact_output) {
            coefficient_rows[[length(coefficient_rows) + 1L]] <- coefficients
        }

        cells <- cell_metrics(candidate_evidence, prediction)
        cells <- merge(
            cells, candidate_baseline_cells,
            by = c("fold", "sample_id"),
            suffixes = c("", "_baseline"), sort = FALSE
        )
        cells[, `:=`(
            motif_id = motif_id,
            threshold = threshold,
            comparison_mode = values$comparison_mode,
            delta_roc_auc = roc_auc - roc_auc_baseline,
            delta_average_precision =
                average_precision - average_precision_baseline,
            delta_log_loss = log_loss - log_loss_baseline,
            delta_brier_score = brier_score - brier_score_baseline
        )]
        if (!values$compact_output) {
            sample_fold_rows[[length(sample_fold_rows) + 1L]] <- cells
        }

        raw <- candidate_evidence[, .(
            retained_anti = sum(retained == 1L & outcome == 1L),
            retained_control = sum(retained == 1L & outcome == 0L),
            absent_anti = sum(retained == 0L & outcome == 1L),
            absent_control = sum(retained == 0L & outcome == 0L)
        ), by = sample_id]
        raw[, `:=`(
            motif_id = motif_id,
            threshold = threshold,
            comparison_mode = values$comparison_mode,
            raw_odds_ratio =
                ((retained_anti + 0.5) * (absent_control + 0.5)) /
                ((retained_control + 0.5) * (absent_anti + 0.5))
        )]
        if (!values$compact_output) {
            raw_rows[[length(raw_rows) + 1L]] <- raw
        }

        threshold_rows[[length(threshold_rows) + 1L]] <- data.table(
            motif_id = motif_id,
            threshold = threshold,
            anchors_total = nrow(anchor_features),
            anchors_retained = counts$positive,
            retained_fraction = counts$positive / nrow(anchor_features),
            anchors_negative_reference = counts$negative,
            negative_reference_fraction =
                counts$negative / nrow(anchor_features),
            anchors_intermediate = counts$intermediate,
            intermediate_fraction =
                counts$intermediate / nrow(anchor_features),
            anchors_compared = counts$positive + counts$negative,
            comparison_fraction =
                (counts$positive + counts$negative) / nrow(anchor_features),
            discordant_observations_total = nrow(evidence),
            discordant_observations = nrow(candidate_evidence),
            comparison_mode = values$comparison_mode,
            negative_reference_threshold = if (
                values$comparison_mode == "fixed-negative-reference"
            ) values$negative_reference_threshold else NA_real_,
            negative_reference_inclusive = if (
                values$comparison_mode == "fixed-negative-reference"
            ) FALSE else NA,
            intermediate_excluded =
                values$comparison_mode == "fixed-negative-reference",
            baseline_macro_roc_auc = candidate_baseline_macro$roc_auc,
            augmented_macro_roc_auc = mean(cells$roc_auc),
            delta_macro_roc_auc = mean(cells$delta_roc_auc),
            baseline_macro_average_precision =
                candidate_baseline_macro$average_precision,
            augmented_macro_average_precision = mean(cells$average_precision),
            delta_macro_average_precision = mean(cells$delta_average_precision),
            baseline_macro_log_loss = candidate_baseline_macro$log_loss,
            augmented_macro_log_loss = mean(cells$log_loss),
            delta_macro_log_loss = mean(cells$delta_log_loss),
            baseline_macro_brier_score = candidate_baseline_macro$brier_score,
            augmented_macro_brier_score = mean(cells$brier_score),
            delta_macro_brier_score = mean(cells$delta_brier_score),
            median_adjusted_odds_ratio = median(coefficients$adjusted_odds_ratio),
            minimum_adjusted_odds_ratio = min(coefficients$adjusted_odds_ratio),
            maximum_adjusted_odds_ratio = max(coefficients$adjusted_odds_ratio),
            median_raw_sample_odds_ratio = median(raw$raw_odds_ratio),
            samples_with_raw_odds_ratio_below_one =
                sum(raw$raw_odds_ratio < 1),
            samples_total = nrow(samples),
            sample_fold_cells = nrow(cells),
            sample_fold_cells_with_roc_auc_gain =
                sum(cells$delta_roc_auc > 0),
            sample_fold_cells_with_log_loss_gain =
                sum(cells$delta_log_loss < 0),
            evaluation_status = "evaluated",
            evaluation_note = NA_character_
        )
        rm(
            prediction, cells, coefficients, raw, candidate_evidence,
            candidate_baseline_prediction, candidate_baseline_cells,
            candidate_baseline_macro
        )
        gc(verbose = FALSE)
    }
}
evidence[, retained := NULL]

threshold_metrics <- rbindlist(threshold_rows, use.names = TRUE, fill = TRUE)
if (nrow(threshold_metrics) == 0L) {
    stop("threshold evaluator produced no candidate rows", call. = FALSE)
}
setorder(threshold_metrics, motif_id, threshold)
fwrite(
    threshold_metrics, paste0(output_prefix, "_threshold_metrics.tsv"), sep = "\t"
)
if (!values$compact_output) {
    sample_fold_metrics <- rbindlist(
        sample_fold_rows, use.names = TRUE, fill = TRUE
    )
    fold_coefficients <- rbindlist(
        coefficient_rows, use.names = TRUE, fill = TRUE
    )
    raw_association <- rbindlist(raw_rows, use.names = TRUE, fill = TRUE)
    if (nrow(sample_fold_metrics) > 0L) {
        setorder(sample_fold_metrics, motif_id, threshold, fold, sample_id)
    }
    if (nrow(fold_coefficients) > 0L) {
        setorder(fold_coefficients, motif_id, threshold, fold)
    }
    if (nrow(raw_association) > 0L) {
        setorder(raw_association, motif_id, threshold, sample_id)
    }
    fwrite(
        sample_fold_metrics,
        paste0(output_prefix, "_sample_fold_metrics.tsv"), sep = "\t"
    )
    fwrite(
        fold_coefficients,
        paste0(output_prefix, "_fold_coefficients.tsv"), sep = "\t"
    )
    fwrite(
        raw_association, paste0(output_prefix, "_raw_association.tsv"), sep = "\t"
    )
}

run_config <- data.table(
    anchor_evidence = values$anchor_evidence,
    cofactor_maxima = values$cofactor_maxima,
    motifs = paste(motifs, collapse = ","),
    evidence_column_scheme = evidence_column_scheme,
    sample_ids = paste(samples$sample_id, collapse = ","),
    sample_count = nrow(samples),
    threshold_specification = values$thresholds,
    threshold_mode = threshold_mode,
    evaluated_thresholds = paste(
        sort(unique(threshold_metrics$threshold)), collapse = ","
    ),
    folds = values$folds,
    fold_definition = "equal_width_contiguous_chromosome_spans",
    chrom_size = chrom_size,
    chrom_size_source = chrom_size_source,
    spline_df = values$spline_df,
    minimum_class_fraction = values$minimum_class_fraction,
    comparison_mode = values$comparison_mode,
    negative_reference_threshold = if (
        values$comparison_mode == "fixed-negative-reference"
    ) values$negative_reference_threshold else NA_real_,
    negative_reference_inclusive = if (
        values$comparison_mode == "fixed-negative-reference"
    ) FALSE else NA,
    intermediate_score_rule = if (
        values$comparison_mode == "fixed-negative-reference"
    ) paste0(
        values$negative_reference_threshold,
        " <= context_score < positive_threshold"
    ) else "none",
    intermediate_handling = if (
        values$comparison_mode == "fixed-negative-reference"
    ) "excluded_from_positive_vs_negative_fit" else
        "assigned_to_threshold_complement",
    compact_output = values$compact_output,
    evaluated_candidate_count = sum(
        threshold_metrics$evaluation_status == "evaluated"
    ),
    unevaluable_candidate_count = sum(
        threshold_metrics$evaluation_status != "evaluated"
    ),
    source_score_floor = source_floors[[1L]],
    context_flank_bp = paste(unique(maxima$context_flank_bp), collapse = ","),
    capture_prefilter_center_bp =
        paste(unique(maxima$capture_prefilter_center_bp), collapse = ","),
    context_distance_metric = paste(
        unique(maxima$context_distance_metric), collapse = ","
    ),
    target = "discordant_anti_p73_only_vs_matched_control_only_immersion",
    threshold_feature = if (
        values$comparison_mode == "fixed-negative-reference"
    ) "positive_threshold_vs_fixed_negative_reference" else
        "retained_context_hit_indicator",
    evidence_note =
        "sample-fold cells are descriptive correlated summaries, not replicates"
)
fwrite(run_config, paste0(output_prefix, "_run_config.tsv"), sep = "\t")

if (values$compact_output ||
    !any(threshold_metrics$evaluation_status == "evaluated")) {
    message("I: wrote TP73 cofactor threshold sweep: ", output_prefix)
    quit(status = 0L)
}
plot_data <- rbindlist(list(
    threshold_metrics[, .(
        motif_id, threshold,
        metric = "Retained anchor fraction",
        value = retained_fraction
    )],
    threshold_metrics[, .(
        motif_id, threshold,
        metric = "Delta ROC AUC vs TP73",
        value = delta_macro_roc_auc
    )],
    threshold_metrics[, .(
        motif_id, threshold,
        metric = "Adjusted odds ratio",
        value = median_adjusted_odds_ratio
    )]
))
plot_data[, metric := factor(
    metric,
    levels = c(
        "Retained anchor fraction", "Delta ROC AUC vs TP73",
        "Adjusted odds ratio"
    )
)]
reference_data <- data.table(
    metric = factor(
        c("Retained anchor fraction", "Delta ROC AUC vs TP73",
          "Adjusted odds ratio"),
        levels = levels(plot_data$metric)
    ),
    reference = c(1, 0, 1)
)
plot <- ggplot(plot_data, aes(x = threshold, y = value, color = motif_id)) +
    geom_hline(
        data = reference_data,
        aes(yintercept = reference),
        linewidth = 0.35,
        color = "#6B7280"
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_grid(metric ~ motif_id, scales = "free_y") +
    scale_color_manual(values = setNames(
        rep(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"),
            length.out = length(motifs)),
        motifs
    )) +
    scale_x_continuous(breaks = sort(unique(threshold_metrics$threshold))) +
    labs(
        title = "Cofactor score-floor sensitivity around TP73 motifs",
        subtitle = paste0(
            "Chromosome 1; five held-out contiguous regions; ",
            "interval edge-distance <= ", unique(maxima$context_flank_bp), " bp"
        ),
        x = "Minimum retained cofactor score",
        y = NULL,
        color = "Motif",
        caption = paste0(
            "Adjusted odds ratio < 1 indicates an inhibitory association. ",
            if (values$comparison_mode == "fixed-negative-reference") {
                paste0(
                    "Negative maxima < ",
                    values$negative_reference_threshold,
                    "; intermediate maxima excluded. "
                )
            } else {
                "Historical threshold-complement comparison. "
            },
            "Exploratory threshold selection; no independent chromosome test."
        )
    ) +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none"
    )
ggsave(
    paste0(output_prefix, "_threshold_sweep.png"), plot,
    width = 11, height = 8, dpi = 160, bg = "white"
)

message("I: wrote TP73 cofactor threshold sweep: ", output_prefix)
