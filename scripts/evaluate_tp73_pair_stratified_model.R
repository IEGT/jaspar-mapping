#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: evaluate_tp73_pair_stratified_model.R --pair-parquet FILE",
        "       --patz1-parquet FILE --tfap2c-parquet FILE",
        "       --pou2f2-parquet FILE --e2f1-parquet FILE",
        "       --output-prefix PATH [options]",
        "",
        "Evaluate whether TP73 tandem-motif architecture and local cofactor",
        "scores predict strict anti-p73 CUT&RUN immersion relative to matched",
        "control. The pair input is tp73_pair_feature Parquet from",
        "build_motif_context.py. Cofactor inputs are exact feature Parquets from",
        "cutandrun_score_calibration --feature-parquet.",
        "",
        "The target uses discordant anti-p73/control pairs: anti-p73-only is 1",
        "and control-only is 0. Shared models are fitted across pair classes,",
        "then evaluated within each class. Genomic folds keep every sample copy",
        "of an anchor in the same held-out fold.",
        "",
        "Options:",
        "  --pair-parquet FILE    TP73 pair-feature Parquet (required)",
        "  --patz1-parquet FILE   TP73/PATZ1 feature Parquet (required)",
        "  --tfap2c-parquet FILE  TP73/TFAP2C feature Parquet (required)",
        "  --pou2f2-parquet FILE  TP73/POU2F2 feature Parquet (required)",
        "  --e2f1-parquet FILE    TP73/E2F1 feature Parquet (required)",
        "  --output-prefix PATH   Basename for result tables and plot (required)",
        "  --context-floor X      Retained cofactor/partner threshold (default: -1)",
        "  --duckdb COMMAND       DuckDB CLI executable (default: duckdb)",
        "  --folds N              Number of block folds (default: 5)",
        "  --fold-mode MODE       interleaved or contiguous (default: interleaved)",
        "  --block-size BP        Genomic block size (default: 5000000)",
        "  --spline-df N          Natural-spline degrees of freedom (default: 4)",
        "  -h, --help             Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) usage()

values <- list(duckdb = "duckdb", folds = 5L, fold_mode = "interleaved",
               block_size = 5000000L, spline_df = 4L, context_floor = -1)
known <- c("--pair-parquet", "--patz1-parquet", "--tfap2c-parquet",
           "--pou2f2-parquet", "--e2f1-parquet", "--output-prefix",
           "--context-floor", "--duckdb", "--folds", "--fold-mode",
           "--block-size", "--spline-df")
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

required <- c("pair_parquet", "patz1_parquet", "tfap2c_parquet",
              "pou2f2_parquet", "e2f1_parquet", "output_prefix")
missing <- required[vapply(required, function(name) is.null(values[[name]]),
                           logical(1))]
if (length(missing) > 0L) {
    writeLines(paste("E: missing required options:",
                     paste(paste0("--", gsub("_", "-", missing)),
                           collapse = ", ")), con = stderr())
    usage(2L)
}
for (name in c("folds", "block_size", "spline_df")) {
    parsed <- suppressWarnings(as.integer(values[[name]]))
    if (is.na(parsed) || parsed <= 0L) {
        stop("--", gsub("_", "-", name), " must be a positive integer",
             call. = FALSE)
    }
    values[[name]] <- parsed
}
if (values$folds < 2L || values$spline_df < 2L) {
    stop("--folds and --spline-df must both be at least 2", call. = FALSE)
}
if (!values$fold_mode %in% c("interleaved", "contiguous")) {
    stop("--fold-mode must be interleaved or contiguous", call. = FALSE)
}
context_floor <- suppressWarnings(as.numeric(values$context_floor))
if (length(context_floor) != 1L || !is.finite(context_floor)) {
    stop("--context-floor must be a finite number", call. = FALSE)
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

path_options <- c("pair_parquet", "patz1_parquet", "tfap2c_parquet",
                  "pou2f2_parquet", "e2f1_parquet")
for (name in path_options) {
    values[[name]] <- normalizePath(values[[name]], mustWork = TRUE)
}
duckdb_path <- Sys.which(values$duckdb)
if (!nzchar(duckdb_path)) stop("DuckDB CLI not found: ", values$duckdb,
                               call. = FALSE)
output_prefix <- values$output_prefix
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

sql_string <- function(value) {
    paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
}

duckdb_query <- function(query, label) {
    sql_path <- tempfile("tp73-pair-model-", fileext = ".sql")
    on.exit(unlink(sql_path), add = TRUE)
    writeLines(query, sql_path, useBytes = TRUE)
    command <- paste(shQuote(duckdb_path), "-csv :memory: <", shQuote(sql_path))
    message("I: streaming ", label, " from DuckDB")
    # DuckDB CSV represents SQL NULL as an empty field. Declaring that here is
    # essential for nullable singleton-only pair fields such as partner gap.
    fread(cmd = command, showProgress = FALSE,
          na.strings = c("", "NA", "NULL"))
}

samples <- data.table(
    sample_id = c("saos2_TA", "saos2_DN", "skmel29_1_TA",
                  "skmel29_1_DN", "skmel29_2_TA", "skmel29_2_DN"),
    anti = c("supported_anti_saos2_TA", "supported_anti_saos2_DN",
             "supported_anti_skmel29_1_TA",
             "supported_anti_skmel29_1_DN",
             "supported_anti_skmel29_2_TA",
             "supported_anti_skmel29_2_DN"),
    control = c("supported_control_saos2_TA", "supported_control_saos2_DN",
                "supported_control_skmel29_1_TA",
                "supported_control_skmel29_1_DN",
                "supported_control_skmel29_2_TA",
                "supported_control_skmel29_2_DN")
)

# Pair features are joined to the orientation selected by the calibrator. This
# prevents a near-palindromic dual-strand alignment from duplicating labels.
base_cte <- paste0(
    "WITH patz1 AS (SELECT * FROM read_parquet(",
    sql_string(values$patz1_parquet), ")),\n",
    "tfap2c AS (SELECT chrom, anchor_start, anchor_score, context_score ",
    "FROM read_parquet(", sql_string(values$tfap2c_parquet), ")),\n",
    "pou2f2 AS (SELECT chrom, anchor_start, anchor_score, context_score ",
    "FROM read_parquet(", sql_string(values$pou2f2_parquet), ")),\n",
    "e2f1 AS (SELECT chrom, anchor_start, anchor_score, context_score ",
    "FROM read_parquet(", sql_string(values$e2f1_parquet), ")),\n",
    "pair_feature AS (SELECT * FROM read_parquet(",
    sql_string(values$pair_parquet), ", hive_partitioning=1)),\n",
    "joined AS (\n",
    "  SELECT p.*, p.anchor_score AS tp73_score, ",
    "p.context_score AS patz1_score,\n",
    "         t.context_score AS tfap2c_score, ",
    "o.context_score AS pou2f2_score,\n",
    "         e.context_score AS e2f1_score, q.pair_class,\n",
    "         q.n_tandem_tp73_partner_loci, ",
    "q.has_multiple_tandem_partner_loci,\n",
    "         q.nearest_tandem_inter_motif_gap_bp, q.best_partner_score,\n",
    "         q.best_pair_min_score, q.best_pair_sum_score, q.tandem_flank_bp\n",
    "  FROM patz1 p\n",
    "  JOIN tfap2c t ON CAST(t.chrom AS VARCHAR) = CAST(p.chrom AS VARCHAR)",
    " AND t.anchor_start = p.anchor_start\n",
    "  JOIN pou2f2 o ON CAST(o.chrom AS VARCHAR) = CAST(p.chrom AS VARCHAR)",
    " AND o.anchor_start = p.anchor_start\n",
    "  JOIN e2f1 e ON CAST(e.chrom AS VARCHAR) = CAST(p.chrom AS VARCHAR)",
    " AND e.anchor_start = p.anchor_start\n",
    "  JOIN pair_feature q\n",
    "    ON CAST(q.chrom AS VARCHAR) = CAST(p.chrom AS VARCHAR)\n",
    "   AND q.start = p.anchor_start AND q.\"end\" = p.anchor_end\n",
    "   AND q.strand = CASE WHEN p.anchor_strand = 1 THEN '+' ELSE '-' END\n",
    "   AND abs(q.score - p.anchor_score) < 1e-5\n",
    "  WHERE p.anchor_score = t.anchor_score\n",
    "    AND p.anchor_score = o.anchor_score\n",
    "    AND p.anchor_score = e.anchor_score\n",
    ")\n"
)

support_selects <- samples[, sprintf(
    paste0(
        "SELECT '%s' AS sample_id, pair_class,\n",
        "       CASE WHEN tp73_score >= 0 THEN 'tp73_ge_0' ",
        "ELSE 'tp73_minus1_to_0' END AS score_stratum,\n",
        "       count(*) AS anchors,\n",
        "       sum(CASE WHEN %s THEN 1 ELSE 0 END) AS anti_supported,\n",
        "       sum(CASE WHEN %s THEN 1 ELSE 0 END) AS control_supported,\n",
        "       sum(CASE WHEN %s AND NOT %s THEN 1 ELSE 0 END) AS anti_only,\n",
        "       sum(CASE WHEN %s AND NOT %s THEN 1 ELSE 0 END) AS control_only,\n",
        "       sum(CASE WHEN %s AND %s THEN 1 ELSE 0 END) AS both_supported\n",
        "FROM joined GROUP BY pair_class, score_stratum"
    ), sample_id, anti, control, anti, control, control, anti, anti, control
)]
support <- duckdb_query(
    paste0(base_cte, paste(support_selects, collapse = "\nUNION ALL\n"), ";\n"),
    "pair-class support summaries"
)
support[, `:=`(
    anti_support_fraction = anti_supported / anchors,
    control_support_fraction = control_supported / anchors,
    support_difference = (anti_supported - control_supported) / anchors,
    discordant = anti_only + control_only,
    conditional_anti_fraction = anti_only / (anti_only + control_only),
    discordant_log_odds = log((anti_only + 0.5) / (control_only + 0.5))
)]
support_all <- support[, copy(.SD)][, score_stratum := "all"]
support_all <- support_all[, lapply(.SD, sum),
                           by = .(sample_id, pair_class, score_stratum),
                           .SDcols = c("anchors", "anti_supported",
                                       "control_supported", "anti_only",
                                       "control_only", "both_supported")]
support_all[, `:=`(
    anti_support_fraction = anti_supported / anchors,
    control_support_fraction = control_supported / anchors,
    support_difference = (anti_supported - control_supported) / anchors,
    discordant = anti_only + control_only,
    conditional_anti_fraction = anti_only / (anti_only + control_only),
    discordant_log_odds = log((anti_only + 0.5) / (control_only + 0.5))
)]
support <- rbindlist(list(support, support_all), use.names = TRUE, fill = TRUE)
fwrite(support, paste0(output_prefix, "_pair_class_sample_support.tsv"),
       sep = "\t")

support_macro <- support[, .(
    anchors_per_sample = median(anchors),
    discordant_observations = sum(discordant),
    anti_support_fraction = mean(anti_support_fraction),
    control_support_fraction = mean(control_support_fraction),
    support_difference = mean(support_difference),
    support_difference_min = min(support_difference),
    support_difference_max = max(support_difference),
    conditional_anti_fraction = mean(conditional_anti_fraction),
    anti_vs_control_discordant_odds = exp(mean(discordant_log_odds))
), by = .(pair_class, score_stratum)]
fwrite(support_macro, paste0(output_prefix, "_pair_class_support.tsv"),
       sep = "\t")

evidence_selects <- samples[, sprintf(
    paste0(
        "SELECT CAST(chrom AS VARCHAR) AS chrom, anchor_start, anchor_end, ",
        "anchor_strand, tp73_score, ",
        "patz1_score, tfap2c_score, pou2f2_score, e2f1_score, pair_class, ",
        "n_tandem_tp73_partner_loci, has_multiple_tandem_partner_loci, ",
        "nearest_tandem_inter_motif_gap_bp, best_partner_score, ",
        "best_pair_min_score, best_pair_sum_score, tandem_flank_bp, ",
        "'%s' AS sample_id, CASE WHEN %s THEN 1 ELSE 0 END AS outcome ",
        "FROM joined WHERE %s <> %s"
    ), sample_id, anti, anti, control
)]
evidence <- duckdb_query(
    paste0(base_cte, paste(evidence_selects, collapse = "\nUNION ALL\n"),
           "\nORDER BY chrom, anchor_start, sample_id;\n"),
    "pair-stratified discordant evidence"
)

# fread preserves DuckDB BIGINT columns as bit64 integer64 values. Genomic
# coordinates at this scale and the small count fields remain exact as ordinary
# doubles, which also gives base is.finite/arithmetic semantics.
integer64_columns <- intersect(
    c("anchor_start", "anchor_end", "n_tandem_tp73_partner_loci",
      "nearest_tandem_inter_motif_gap_bp", "tandem_flank_bp"),
    names(evidence)
)
for (column in integer64_columns) {
    original_missing <- sum(is.na(evidence[[column]]))
    converted <- suppressWarnings(as.numeric(evidence[[column]]))
    introduced <- sum(is.na(converted)) - original_missing
    if (introduced > 0L) {
        stop("numeric conversion introduced ", introduced,
             " missing values in ", column, " (source class: ",
             paste(class(evidence[[column]]), collapse = "/"), ")",
             call. = FALSE)
    }
    set(evidence, j = column, value = converted)
}

score_columns <- c("tp73_score", "patz1_score", "tfap2c_score",
                   "pou2f2_score", "e2f1_score")
expected_columns <- c("chrom", "anchor_start", "anchor_end", "anchor_strand",
                      score_columns, "pair_class", "sample_id", "outcome")
missing_columns <- setdiff(expected_columns, names(evidence))
if (length(missing_columns) > 0L) {
    stop("DuckDB result lacks columns: ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
}
if (nrow(evidence) == 0L || evidence[, any(!outcome %in% 0:1)] ||
    evidence[, any(!is.finite(unlist(.SD))), .SDcols = score_columns]) {
    stop("model evidence is empty or contains invalid values", call. = FALSE)
}
duplicate_evidence <- evidence[, .N, by = .(
    chrom, anchor_start, anchor_strand, sample_id
)][N > 1L]
if (nrow(duplicate_evidence) > 0L) {
    stop("feature joins duplicated chromosome/anchor/sample rows", call. = FALSE)
}
rm(duplicate_evidence)

pair_levels <- c("singleton", "tandem_same_orientation",
                 "tandem_opposite_orientation",
                 "tandem_orientation_ambiguous", "tandem_mixed_orientation")
unknown_pair_classes <- setdiff(unique(evidence$pair_class), pair_levels)
if (length(unknown_pair_classes) > 0L) {
    stop("unknown pair classes: ", paste(unknown_pair_classes, collapse = ", "),
         call. = FALSE)
}
evidence[, pair_class := factor(pair_class, levels = pair_levels)]
evidence[, sample_id := factor(sample_id, levels = samples$sample_id)]
if (evidence[, any(is.na(pair_class) | is.na(sample_id))]) {
    stop("pair class or sample factor conversion failed", call. = FALSE)
}
if (values$fold_mode == "interleaved") {
    evidence[, fold := as.integer((anchor_start %/% values$block_size) %%
                                  values$folds) + 1L]
    fold_definition <- paste0(
        "cyclic_", values$block_size, "bp_blocks_within_chromosome"
    )
} else {
    chromosome_bounds <- unique(evidence[, .(chrom, anchor_start)])[, .(
        fold_min = min(anchor_start),
        fold_span = max(anchor_start) - min(anchor_start) + 1
    ), by = chrom]
    evidence[chromosome_bounds, on = "chrom", `:=`(
        fold_min = i.fold_min,
        fold_span = i.fold_span
    )]
    evidence[, fold := pmin(
        values$folds,
        as.integer(floor(
            (anchor_start - fold_min) * values$folds / fold_span
        )) + 1L
    )]
    evidence[, c("fold_min", "fold_span") := NULL]
    fold_definition <- paste0(
        values$folds, "_equal_width_contiguous_spans_within_each_chromosome"
    )
}
fold_loci <- unique(
    evidence[, .(chrom, fold, anchor_start, anchor_strand)]
)
fold_manifest <- fold_loci[, .(
    min_anchor_start = min(anchor_start),
    max_anchor_start = max(anchor_start),
    genomic_span_bp = max(anchor_start) - min(anchor_start) + 1,
    anchor_orientation_loci = .N
), by = .(chrom, fold)]
fold_manifest[, `:=`(
    fold_mode = values$fold_mode,
    fold_definition = fold_definition
)]
setorder(fold_manifest, chrom, fold)
fwrite(fold_manifest, paste0(output_prefix, "_fold_manifest.tsv"), sep = "\t")

evidence[, has_tandem := as.integer(pair_class != "singleton")]
if (evidence[has_tandem == 1L,
             any(!is.finite(best_partner_score) |
                 !is.finite(nearest_tandem_inter_motif_gap_bp))]) {
    stop("a tandem pair lacks partner score or gap", call. = FALSE)
}
evidence[, pair_partner_strength := fifelse(
    has_tandem == 1L,
    log1p(pmax(best_partner_score - context_floor, 0)), 0
)]
evidence[, pair_closeness := fifelse(
    has_tandem == 1L & tandem_flank_bp > 0,
    pmax(tandem_flank_bp - nearest_tandem_inter_motif_gap_bp, 0) /
        tandem_flank_bp,
    as.numeric(has_tandem)
)]
evidence[, pair_multiple_strength := log1p(pmax(
    n_tandem_tp73_partner_loci - 1, 0
))]

cofactor_columns <- c("patz1_score", "tfap2c_score", "pou2f2_score",
                      "e2f1_score")
anchor_features <- unique(
    evidence[, c("chrom", "anchor_start", "anchor_strand", score_columns,
                 "pair_class", "n_tandem_tp73_partner_loci",
                 "best_partner_score", "nearest_tandem_inter_motif_gap_bp"),
             with = FALSE],
    by = c("chrom", "anchor_start", "anchor_strand")
)

format_r_number <- function(value) {
    format(value, digits = 17L, scientific = TRUE, trim = TRUE)
}
context_terms <- setNames(vector("list", length(cofactor_columns)),
                          cofactor_columns)
context_basis_rows <- vector("list", length(cofactor_columns))
for (column in cofactor_columns) {
    retained_column <- paste0(column, "_retained")
    excess_column <- paste0(column, "_excess")
    strength_column <- paste0(column, "_strength")
    evidence[, (retained_column) := as.integer(get(column) >= context_floor)]
    evidence[, (excess_column) := pmax(get(column) - context_floor, 0)]
    evidence[, (strength_column) := log1p(get(excess_column))]

    retained_values <- anchor_features[[column]][
        anchor_features[[column]] >= context_floor
    ]
    excess_values <- retained_values - context_floor
    boundary_max <- max(excess_values)
    if (!is.finite(boundary_max) || boundary_max <= 0) {
        stop("cofactor has no score variation above --context-floor: ", column,
             call. = FALSE)
    }
    probabilities <- seq_len(values$spline_df - 1L) / values$spline_df
    knots <- unique(as.numeric(quantile(
        excess_values, probs = probabilities, names = FALSE, type = 7
    )))
    knots <- knots[knots > 0 & knots < boundary_max]
    if (length(knots) == 0L) {
        context_terms[[column]] <- paste(retained_column, excess_column,
                                          sep = " + ")
        encoding <- "retained_indicator_plus_linear_score_excess"
    } else {
        context_terms[[column]] <- paste0(
            retained_column, " + splines::ns(", excess_column,
            ", knots = c(", paste(format_r_number(knots), collapse = ", "),
            "), Boundary.knots = c(0, ", format_r_number(boundary_max), "))"
        )
        encoding <- "retained_indicator_plus_score_excess_spline"
    }
    context_basis_rows[[column]] <- data.table(
        motif = sub("_score$", "", column), context_floor = context_floor,
        anchors_total = nrow(anchor_features),
        anchors_retained = length(retained_values),
        retained_fraction = length(retained_values) / nrow(anchor_features),
        boundary_max_excess = boundary_max,
        interior_knots_excess = paste(format_r_number(knots), collapse = ","),
        encoding = encoding
    )
}
context_basis <- rbindlist(context_basis_rows)
fwrite(context_basis, paste0(output_prefix, "_context_basis.tsv"), sep = "\t")

tp73_term <- sprintf("splines::ns(tp73_score, df = %d)", values$spline_df)
pair_terms <- paste("pair_class + pair_partner_strength + pair_closeness +",
                    "pair_multiple_strength")
cofactor_terms <- paste(unlist(context_terms), collapse = " + ")
interaction_terms <- paste(unlist(lapply(cofactor_columns, function(column) {
    c(paste0("pair_class:", column, "_retained"),
      paste0("pair_class:", column, "_strength"))
})), collapse = " + ")

model_specs <- list(
    "Sample only" = "outcome ~ sample_id",
    "TP73 score" = paste("outcome ~ sample_id +", tp73_term),
    "Pair architecture only" = paste("outcome ~ sample_id +", pair_terms),
    "TP73 + pair architecture" = paste(
        "outcome ~ sample_id +", tp73_term, "+", pair_terms
    ),
    "Cofactors only" = paste("outcome ~ sample_id +", cofactor_terms),
    "TP73 + cofactors" = paste(
        "outcome ~ sample_id +", tp73_term, "+", cofactor_terms
    ),
    "TP73 + pair + cofactors" = paste(
        "outcome ~ sample_id +", tp73_term, "+", pair_terms, "+",
        cofactor_terms
    ),
    "TP73 + pair x cofactors" = paste(
        "outcome ~ sample_id +", tp73_term, "+", pair_terms, "+",
        cofactor_terms, "+", interaction_terms
    )
)
fwrite(
    data.table(model = names(model_specs), formula = unlist(model_specs,
                                                            use.names = FALSE)),
    paste0(output_prefix, "_model_specifications.tsv"), sep = "\t"
)

roc_auc <- function(outcome, prediction) {
    positives <- as.double(sum(outcome == 1L))
    negatives <- as.double(sum(outcome == 0L))
    if (positives == 0 || negatives == 0) return(NA_real_)
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
    grouped[, `:=`(cumulative_positive = cumsum(positive),
                   cumulative_total = cumsum(total))]
    sum((grouped$positive / positives) *
        (grouped$cumulative_positive / grouped$cumulative_total))
}

metric_row <- function(model, scope, outcome, prediction) {
    if (length(outcome) == 0L) {
        return(data.table(model = model, scope = scope, n = 0L,
                          anti_only_fraction = NA_real_, roc_auc = NA_real_,
                          average_precision = NA_real_, log_loss = NA_real_,
                          brier_score = NA_real_))
    }
    epsilon <- 1e-15
    bounded <- pmin(pmax(prediction, epsilon), 1 - epsilon)
    data.table(
        model = model, scope = scope, n = length(outcome),
        anti_only_fraction = mean(outcome == 1L),
        roc_auc = roc_auc(outcome, prediction),
        average_precision = average_precision(outcome, prediction),
        log_loss = -mean(outcome * log(bounded) +
                         (1 - outcome) * log(1 - bounded)),
        brier_score = mean((outcome - prediction)^2)
    )
}

fold_check <- evidence[, .(classes = uniqueN(outcome),
                           samples = uniqueN(sample_id)), by = fold]
if (fold_check[, any(classes != 2L | samples != nrow(samples))]) {
    stop("one or more folds lacks both outcomes or a sample", call. = FALSE)
}

predictions <- vector("list", length(model_specs))
names(predictions) <- names(model_specs)
fold_metrics <- list()
for (model_name in names(model_specs)) {
    message("I: fitting ", model_name)
    formula <- as.formula(model_specs[[model_name]])
    out_of_fold <- rep(NA_real_, nrow(evidence))
    for (held_out in seq_len(values$folds)) {
        training <- evidence[fold != held_out]
        testing_indexes <- which(evidence$fold == held_out)
        fit <- glm(formula, data = training, family = binomial(),
                   control = glm.control(maxit = 50L, epsilon = 1e-8))
        if (!isTRUE(fit$converged)) {
            stop("model did not converge: ", model_name, ", fold ", held_out,
                 call. = FALSE)
        }
        predicted <- predict(fit, newdata = evidence[testing_indexes],
                             type = "response")
        if (any(!is.finite(predicted))) {
            stop("non-finite predictions: ", model_name, ", fold ", held_out,
                 call. = FALSE)
        }
        out_of_fold[testing_indexes] <- predicted
        fold_metrics[[length(fold_metrics) + 1L]] <- metric_row(
            model_name, paste0("fold_", held_out),
            evidence$outcome[testing_indexes], predicted
        )[, fold := held_out]
        rm(fit, training, predicted)
        gc(verbose = FALSE)
    }
    predictions[[model_name]] <- out_of_fold
}

sample_fold_metrics <- rbindlist(lapply(names(predictions), function(model_name) {
    prediction <- predictions[[model_name]]
    rbindlist(lapply(seq_len(values$folds), function(held_out) {
        rbindlist(lapply(samples$sample_id, function(sample) {
            selected <- evidence$fold == held_out & evidence$sample_id == sample
            metric_row(model_name, sample, evidence$outcome[selected],
                       prediction[selected])[, fold := held_out]
        }))
    }))
}))

metric_means <- function(rows) {
    columns <- c("anti_only_fraction", "roc_auc", "average_precision",
                 "log_loss", "brier_score")
    result <- lapply(rows[, ..columns], function(column) mean(column, na.rm = TRUE))
    result <- lapply(result, function(value) if (is.nan(value)) NA_real_ else value)
    as.data.table(result)
}

performance <- rbindlist(lapply(names(predictions), function(model_name) {
    prediction <- predictions[[model_name]]
    pooled <- metric_row(model_name, "pooled", evidence$outcome, prediction)
    cells <- sample_fold_metrics[model == model_name]
    macro <- metric_means(cells)
    macro[, `:=`(model = model_name, scope = "sample_fold_macro_mean",
                 n = sum(cells$n))]
    setcolorder(macro, names(pooled))
    rbindlist(list(pooled, macro))
}))
fold_metrics <- rbindlist(fold_metrics, fill = TRUE)
fwrite(performance, paste0(output_prefix, "_model_metrics.tsv"), sep = "\t")
fwrite(fold_metrics, paste0(output_prefix, "_fold_metrics.tsv"), sep = "\t")
fwrite(sample_fold_metrics,
       paste0(output_prefix, "_sample_fold_metrics.tsv"), sep = "\t")

pair_sample_fold_metrics <- rbindlist(lapply(names(predictions),
                                              function(model_name) {
    prediction <- predictions[[model_name]]
    rbindlist(lapply(pair_levels, function(pair) {
        rbindlist(lapply(seq_len(values$folds), function(held_out) {
            rbindlist(lapply(samples$sample_id, function(sample) {
                selected <- evidence$pair_class == pair &
                    evidence$fold == held_out & evidence$sample_id == sample
                metric_row(model_name, sample, evidence$outcome[selected],
                           prediction[selected])[
                    , `:=`(pair_class = pair, fold = held_out)
                ]
            }))
        }))
    }))
}))
pair_metrics <- rbindlist(lapply(names(predictions), function(model_name) {
    prediction <- predictions[[model_name]]
    rbindlist(lapply(pair_levels, function(pair) {
        selected <- evidence$pair_class == pair
        pooled <- metric_row(model_name, "pooled", evidence$outcome[selected],
                             prediction[selected])[, pair_class := pair]
        cells <- pair_sample_fold_metrics[
            model == model_name & pair_class == pair & n > 0
        ]
        macro <- metric_means(cells)
        macro[, `:=`(model = model_name, scope = "sample_fold_macro_mean",
                     n = sum(cells$n), pair_class = pair)]
        setcolorder(macro, names(pooled))
        rbindlist(list(pooled, macro))
    }))
}))
fwrite(pair_metrics, paste0(output_prefix, "_pair_class_model_metrics.tsv"),
       sep = "\t")
fwrite(pair_sample_fold_metrics,
       paste0(output_prefix, "_pair_class_sample_fold_metrics.tsv"), sep = "\t")

cofactor_names <- c("PATZ1", "TFAP2C", "POU2F2", "E2F1")
retained_columns <- paste0(tolower(cofactor_names), "_score_retained")
combination_mask <- Reduce(`+`, Map(
    function(column, bit) evidence[[column]] * bit,
    retained_columns, 2 ^ (seq_along(retained_columns) - 1L)
))
combination_labels <- vapply(0:15, function(mask) {
    present <- cofactor_names[bitwAnd(mask, 2 ^ (0:3)) != 0]
    if (length(present) == 0L) "none" else paste(present, collapse = "+")
}, character(1))
evidence[, cofactor_combination := factor(
    combination_labels[combination_mask + 1L], levels = combination_labels
)]

combination_candidates <- performance[
    scope == "sample_fold_macro_mean" &
    model %in% c("TP73 + pair + cofactors", "TP73 + pair x cofactors")
]
best_model <- combination_candidates[which.max(roc_auc), as.character(model)]
evidence[, best_prediction := predictions[[best_model]]]
combination_cells <- evidence[, .(
    n = .N, anti_only = sum(outcome == 1L),
    observed_anti_fraction = mean(outcome == 1L),
    mean_held_out_prediction = mean(best_prediction)
), by = .(pair_class, cofactor_combination, sample_id, fold)]
combination_summary <- combination_cells[, .(
    n = sum(n), anti_only = sum(anti_only), sample_fold_cells = .N,
    observed_anti_fraction = sum(anti_only) / sum(n),
    observed_anti_fraction_macro = mean(observed_anti_fraction),
    mean_held_out_prediction = weighted.mean(mean_held_out_prediction, n),
    mean_held_out_prediction_macro = mean(mean_held_out_prediction)
), by = .(pair_class, cofactor_combination)]
fwrite(combination_summary,
       paste0(output_prefix, "_pair_cofactor_combinations.tsv"), sep = "\t")

run_config <- data.table(
    pair_parquet = values$pair_parquet,
    patz1_parquet = values$patz1_parquet,
    tfap2c_parquet = values$tfap2c_parquet,
    pou2f2_parquet = values$pou2f2_parquet,
    e2f1_parquet = values$e2f1_parquet,
    context_floor = context_floor,
    model_anchors = nrow(anchor_features),
    discordant_observations = nrow(evidence),
    folds = values$folds,
    fold_mode = values$fold_mode,
    fold_definition = fold_definition,
    block_size = values$block_size,
    block_size_applies = values$fold_mode == "interleaved",
    spline_df = values$spline_df,
    target = "anti_p73_only_vs_matched_control_only_strict_immersion",
    model_scope = "shared_model_with_pair_stratified_evaluation",
    cofactor_combination_prediction_model = best_model
)
fwrite(run_config, paste0(output_prefix, "_run_config.tsv"), sep = "\t")

pair_labels <- c(
    singleton = "Singleton",
    tandem_same_orientation = "Tandem: same",
    tandem_opposite_orientation = "Tandem: opposite",
    tandem_orientation_ambiguous = "Tandem: ambiguous",
    tandem_mixed_orientation = "Tandem: mixed"
)
plot_support <- support_macro[score_stratum == "all"]
plot_support[, pair_label := factor(pair_labels[pair_class],
                                    levels = unname(pair_labels))]
p_support <- ggplot(plot_support,
                    aes(x = pair_label, y = support_difference)) +
    geom_hline(yintercept = 0, colour = "#777777", linewidth = 0.4) +
    geom_errorbar(aes(ymin = support_difference_min,
                      ymax = support_difference_max),
                  width = 0.16, colour = "#0072B2") +
    geom_point(size = 2.5, colour = "#0072B2") +
    scale_y_continuous(labels = function(x) sprintf("%+.2f%%", 100 * x)) +
    labs(title = "TP73-specific support by motif architecture",
         subtitle = "Mean anti-p73 minus control immersion; bars span six samples",
         x = NULL, y = "Support difference") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 28, hjust = 1),
          plot.title = element_text(face = "bold"))

model_order <- names(model_specs)
plot_models <- performance[scope == "sample_fold_macro_mean"]
plot_models[, model := factor(model, levels = rev(model_order))]
plot_models <- melt(
    plot_models, id.vars = "model",
    measure.vars = c("roc_auc", "average_precision"),
    variable.name = "metric", value.name = "value"
)
plot_models[, metric := factor(metric,
    levels = c("roc_auc", "average_precision"),
    labels = c("ROC AUC", "Average precision"))]
p_models <- ggplot(plot_models, aes(x = value, y = model, fill = metric)) +
    geom_col(width = 0.72, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.3f", value)), hjust = -0.08, size = 3) +
    facet_wrap(~metric, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = c("#009E73", "#E69F00")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Held-out prediction",
         subtitle = paste0(
             "Macro-mean across six samples and ", values$folds, " ",
             values$fold_mode, " genomic folds"
         ),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))

plot_combinations <- combination_summary[n >= 100]
plot_combinations[, pair_label := factor(pair_labels[as.character(pair_class)],
                                         levels = rev(unname(pair_labels)))]
p_combinations <- ggplot(
    plot_combinations,
    aes(x = cofactor_combination, y = pair_label,
        fill = mean_held_out_prediction_macro)
) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_viridis_c(option = "C", labels = function(x) sprintf("%.0f%%", 100*x),
                         na.value = "#EEEEEE") +
    labs(title = "Pair architecture with retained local cofactors",
         subtitle = paste0(
             "Out-of-fold predicted anti-p73-only fraction; cofactors retained at score >= ",
             format(context_floor, trim = TRUE), "; cells with at least 100 observations"
         ),
         x = NULL, y = NULL, fill = "Predicted") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
          panel.grid = element_blank(), plot.title = element_text(face = "bold"))

png(paste0(output_prefix, "_overview.png"), width = 2800, height = 2000,
    res = 180, bg = "white")
grid::grid.newpage()
layout <- grid::grid.layout(2, 2, heights = grid::unit(c(1, 1.15), "null"))
grid::pushViewport(grid::viewport(layout = layout))
print(p_support, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_models, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_combinations,
      vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2))
grid::popViewport()
dev.off()

message("I: wrote pair-stratified results for ", nrow(anchor_features),
        " anchors and ", nrow(evidence), " discordant observations to ",
        output_prefix, "_*")
