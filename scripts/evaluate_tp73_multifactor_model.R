#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: evaluate_tp73_multifactor_model.R --patz1-parquet FILE",
        "       --tfap2c-parquet FILE --pou2f2-parquet FILE",
        "       --output-prefix PATH [options]",
        "",
        "Compare held-out predictors of TP73-specific CUT&RUN immersion. Each input",
        "must contain the same TP73 anchors and the exact best local cofactor score",
        "from cutandrun_score_calibration --feature-parquet. The PATZ1 table must",
        "also contain six anti-p73 and six matched-control support columns.",
        "",
        "The binary target uses discordant pairs only: anti-p73-only immersion is 1",
        "and matched-control-only immersion is 0. Sites supported by both or neither",
        "are excluded. Five chromosome-block folds keep every copy of a genomic site",
        "in one fold. Models are nonlinear additive logistic regressions with sample",
        "identity as a nuisance term.",
        "When --context-floor is set, each cofactor is encoded as a retained-hit",
        "indicator plus a nonlinear score excess above that floor. Scores below",
        "the floor are therefore indistinguishable, as they would be after storage.",
        "",
        "Options:",
        "  --patz1-parquet FILE   TP73/PATZ1 exact feature Parquet (required)",
        "  --tfap2c-parquet FILE  TP73/TFAP2C exact feature Parquet (required)",
        "  --pou2f2-parquet FILE  TP73/POU2F2 exact feature Parquet (required)",
        "  --output-prefix PATH   Basename for result tables and plot (required)",
        "  --duckdb COMMAND       DuckDB CLI executable (default: duckdb)",
        "  --folds N              Number of chromosome-block folds (default: 5)",
        "  --block-size BP        Contiguous genomic block size (default: 5000000)",
        "  --spline-df N          Natural-spline degrees of freedom (default: 4)",
        "  --context-floor X      Censor cofactor maxima below X; use 'none' for",
        "                         exact scores (default: none)",
        "  -h, --help             Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) {
    usage()
}

values <- list(duckdb = "duckdb", folds = 5L, block_size = 5000000L,
               spline_df = 4L, context_floor = "none")
known <- c("--patz1-parquet", "--tfap2c-parquet", "--pou2f2-parquet",
           "--output-prefix", "--duckdb", "--folds", "--block-size",
           "--spline-df", "--context-floor")
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

required <- c("patz1_parquet", "tfap2c_parquet", "pou2f2_parquet",
              "output_prefix")
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
if (values$folds < 2L) {
    stop("--folds must be at least 2", call. = FALSE)
}
if (values$spline_df < 2L) {
    stop("--spline-df must be at least 2", call. = FALSE)
}
if (tolower(values$context_floor) == "none") {
    context_floor <- NULL
} else {
    context_floor <- suppressWarnings(as.numeric(values$context_floor))
    if (length(context_floor) != 1L || !is.finite(context_floor)) {
        stop("--context-floor must be a finite number or 'none'", call. = FALSE)
    }
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

for (name in c("patz1_parquet", "tfap2c_parquet", "pou2f2_parquet")) {
    values[[name]] <- normalizePath(values[[name]], mustWork = TRUE)
}
duckdb_path <- Sys.which(values$duckdb)
if (!nzchar(duckdb_path)) {
    stop("DuckDB CLI not found: ", values$duckdb, call. = FALSE)
}
output_prefix <- values$output_prefix
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

sql_string <- function(value) {
    paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
}

samples <- data.table(
    sample_id = c("saos2_TA", "saos2_DN", "skmel29_1_TA",
                  "skmel29_1_DN", "skmel29_2_TA", "skmel29_2_DN"),
    anti = c("supported_anti_saos2_TA", "supported_anti_saos2_DN",
             "supported_anti_skmel29_1_TA", "supported_anti_skmel29_1_DN",
             "supported_anti_skmel29_2_TA", "supported_anti_skmel29_2_DN"),
    control = c("supported_control_saos2_TA", "supported_control_saos2_DN",
                "supported_control_skmel29_1_TA",
                "supported_control_skmel29_1_DN",
                "supported_control_skmel29_2_TA",
                "supported_control_skmel29_2_DN")
)

pair_queries <- samples[, sprintf(
    paste(
        "SELECT anchor_start, tp73_score, patz1_score, tfap2c_score,",
        "       pou2f2_score, '%s' AS sample_id,",
        "       CASE WHEN %s THEN 1 ELSE 0 END AS outcome",
        "FROM joined WHERE %s <> %s"
    ), sample_id, anti, anti, control
)]
query <- paste0(
    "WITH patz1 AS (SELECT * FROM read_parquet(",
    sql_string(values$patz1_parquet), ")),\n",
    "tfap2c AS (SELECT anchor_start, anchor_score, context_score ",
    "FROM read_parquet(", sql_string(values$tfap2c_parquet), ")),\n",
    "pou2f2 AS (SELECT anchor_start, anchor_score, context_score ",
    "FROM read_parquet(", sql_string(values$pou2f2_parquet), ")),\n",
    "joined AS (\n",
    "  SELECT p.*, p.anchor_score AS tp73_score, ",
    "p.context_score AS patz1_score,\n",
    "         t.context_score AS tfap2c_score, ",
    "o.context_score AS pou2f2_score\n",
    "  FROM patz1 p JOIN tfap2c t USING (anchor_start) ",
    "JOIN pou2f2 o USING (anchor_start)\n",
    "  WHERE p.anchor_score = t.anchor_score ",
    "AND p.anchor_score = o.anchor_score\n",
    ")\n",
    paste(pair_queries, collapse = "\nUNION ALL\n"),
    "\nORDER BY anchor_start, sample_id;\n"
)
sql_path <- tempfile("tp73-multifactor-", fileext = ".sql")
on.exit(unlink(sql_path), add = TRUE)
writeLines(query, sql_path, useBytes = TRUE)
duckdb_command <- paste(shQuote(duckdb_path), "-csv :memory: <",
                        shQuote(sql_path))
message("I: streaming joined discordant evidence from DuckDB")
evidence <- fread(cmd = duckdb_command, showProgress = FALSE)

expected_columns <- c("anchor_start", "tp73_score", "patz1_score",
                      "tfap2c_score", "pou2f2_score", "sample_id", "outcome")
missing_columns <- setdiff(expected_columns, names(evidence))
if (length(missing_columns) > 0L) {
    stop("DuckDB result lacks columns: ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
}
if (nrow(evidence) == 0L || evidence[, any(!outcome %in% 0:1)]) {
    stop("discordant evidence is empty or has an invalid outcome", call. = FALSE)
}
score_columns <- c("tp73_score", "patz1_score", "tfap2c_score",
                   "pou2f2_score")
if (evidence[, any(!is.finite(unlist(.SD))), .SDcols = score_columns]) {
    stop("model score columns contain non-finite values", call. = FALSE)
}
context_columns <- c("patz1_score", "tfap2c_score", "pou2f2_score")
anchor_features <- unique(evidence[, c("anchor_start", score_columns), with = FALSE],
                          by = "anchor_start")
anchor_strata <- list(
    all = rep(TRUE, nrow(anchor_features)),
    tp73_ge_0 = anchor_features$tp73_score >= 0,
    tp73_minus1_to_0 = anchor_features$tp73_score >= -1 &
                       anchor_features$tp73_score < 0
)
context_availability <- rbindlist(lapply(names(anchor_strata), function(stratum) {
    selected <- anchor_strata[[stratum]]
    if (!any(selected)) return(NULL)
    rbindlist(lapply(c(-1, 0), function(threshold) {
        rbindlist(lapply(context_columns, function(column) {
            values_at_anchors <- anchor_features[[column]][selected]
            data.table(
                stratum = stratum,
                threshold = threshold,
                motif = sub("_score$", "", column),
                anchors_total = length(values_at_anchors),
                anchors_with_retained_context = sum(values_at_anchors >= threshold),
                retained_fraction = mean(values_at_anchors >= threshold)
            )
        }))
    }))
}))
fwrite(context_availability,
       paste0(output_prefix, "_context_availability.tsv"), sep = "\t")
evidence[, sample_id := factor(sample_id, levels = samples$sample_id)]
if (evidence[, any(is.na(sample_id))]) {
    stop("model evidence contains an unknown sample ID", call. = FALSE)
}
evidence[, fold := as.integer((anchor_start %/% values$block_size) %%
                              values$folds) + 1L]

label_counts <- evidence[, .(
    n = .N,
    anti_only = sum(outcome == 1L),
    control_only = sum(outcome == 0L),
    anti_only_fraction = mean(outcome == 1L)
), by = sample_id]
fwrite(label_counts, paste0(output_prefix, "_label_counts.tsv"), sep = "\t")

run_config <- data.table(
    patz1_parquet = values$patz1_parquet,
    tfap2c_parquet = values$tfap2c_parquet,
    pou2f2_parquet = values$pou2f2_parquet,
    context_floor = if (is.null(context_floor)) "none" else
                    format(context_floor, scientific = FALSE, trim = TRUE),
    anchor_score_min = min(anchor_features$tp73_score),
    anchor_score_max = max(anchor_features$tp73_score),
    model_anchors = nrow(anchor_features),
    discordant_observations = nrow(evidence),
    folds = values$folds,
    block_size = values$block_size,
    spline_df = values$spline_df,
    context_encoding = if (is.null(context_floor)) "exact_score_spline" else
                       "retained_indicator_plus_score_excess_spline"
)
fwrite(run_config, paste0(output_prefix, "_run_config.tsv"), sep = "\t")

fold_check <- evidence[, .(classes = uniqueN(outcome), samples = uniqueN(sample_id)),
                       by = fold]
if (fold_check[, any(classes != 2L | samples != nrow(samples))]) {
    stop("one or more chromosome folds lack both outcomes or a sample", call. = FALSE)
}

df <- values$spline_df
ns_term <- function(column) {
    sprintf("splines::ns(%s, df = %d)", column, df)
}
format_r_number <- function(value) {
    format(value, digits = 17L, scientific = TRUE, trim = TRUE)
}
context_terms <- setNames(vector("list", length(context_columns)),
                          context_columns)
if (is.null(context_floor)) {
    context_terms <- setNames(lapply(context_columns, ns_term), context_columns)
} else {
    context_basis_rows <- vector("list", length(context_columns))
    names(context_basis_rows) <- context_columns
    for (column in context_columns) {
        retained_column <- paste0(column, "_retained")
        excess_column <- paste0(column, "_excess")
        evidence[, (retained_column) := as.integer(get(column) >= context_floor)]
        evidence[, (excess_column) := pmax(get(column) - context_floor, 0)]

        anchor_values <- anchor_features[[column]]
        retained_values <- anchor_values[anchor_values >= context_floor]
        excess_values <- retained_values - context_floor
        boundary_max <- max(excess_values)
        if (!is.finite(boundary_max) || boundary_max <= 0) {
            stop("cofactor has no score variation above --context-floor: ", column,
                 call. = FALSE)
        }
        probabilities <- seq_len(df - 1L) / df
        knots <- unique(as.numeric(quantile(
            excess_values, probs = probabilities, names = FALSE, type = 7
        )))
        knots <- knots[knots > 0 & knots < boundary_max]
        if (length(knots) == 0L) {
            context_terms[[column]] <- paste(retained_column, excess_column,
                                              sep = " + ")
            encoding <- "retained_indicator_plus_linear_score_excess"
        } else {
            knot_expression <- paste(format_r_number(knots), collapse = ", ")
            context_terms[[column]] <- paste0(
                retained_column, " + splines::ns(", excess_column,
                ", knots = c(", knot_expression,
                "), Boundary.knots = c(0, ", format_r_number(boundary_max), "))"
            )
            encoding <- "retained_indicator_plus_score_excess_spline"
        }
        context_basis_rows[[column]] <- data.table(
            motif = sub("_score$", "", column),
            context_floor = context_floor,
            anchors_total = length(anchor_values),
            anchors_retained = length(retained_values),
            retained_fraction = length(retained_values) / length(anchor_values),
            boundary_max_excess = boundary_max,
            interior_knots_excess = paste(format_r_number(knots), collapse = ","),
            encoding = encoding
        )
    }
    context_basis <- rbindlist(context_basis_rows)
    fwrite(context_basis, paste0(output_prefix, "_context_basis.tsv"), sep = "\t")
}
tp73 <- ns_term("tp73_score")
patz1 <- context_terms[["patz1_score"]]
tfap2c <- context_terms[["tfap2c_score"]]
pou2f2 <- context_terms[["pou2f2_score"]]
model_specs <- list(
    "Sample only" = "outcome ~ sample_id",
    "TP73" = paste("outcome ~ sample_id +", tp73),
    "Three cofactors" = paste("outcome ~ sample_id +", patz1, "+", tfap2c,
                              "+", pou2f2),
    "TP73 + PATZ1" = paste("outcome ~ sample_id +", tp73, "+", patz1),
    "TP73 + TFAP2C" = paste("outcome ~ sample_id +", tp73, "+", tfap2c),
    "TP73 + POU2F2" = paste("outcome ~ sample_id +", tp73, "+", pou2f2),
    "TP73 + all three" = paste("outcome ~ sample_id +", tp73, "+", patz1,
                               "+", tfap2c, "+", pou2f2)
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
    grouped[, `:=`(cumulative_positive = cumsum(positive),
                   cumulative_total = cumsum(total))]
    sum((grouped$positive / positives) *
        (grouped$cumulative_positive / grouped$cumulative_total))
}

metric_row <- function(model, scope, outcome, prediction) {
    epsilon <- 1e-15
    bounded <- pmin(pmax(prediction, epsilon), 1 - epsilon)
    data.table(
        model = model,
        scope = scope,
        n = length(outcome),
        anti_only_fraction = mean(outcome == 1L),
        roc_auc = roc_auc(outcome, prediction),
        average_precision = average_precision(outcome, prediction),
        log_loss = -mean(outcome * log(bounded) +
                         (1 - outcome) * log(1 - bounded)),
        brier_score = mean((outcome - prediction)^2)
    )
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
            stop("model produced non-finite predictions: ", model_name,
                 ", fold ", held_out, call. = FALSE)
        }
        out_of_fold[testing_indexes] <- predicted
        fold_metrics[[length(fold_metrics) + 1L]] <- metric_row(
            model_name, paste0("fold_", held_out),
            evidence$outcome[testing_indexes], predicted
        )[, fold := held_out]
        rm(fit, training, predicted)
        gc(verbose = FALSE)
    }
    if (any(!is.finite(out_of_fold))) {
        stop("out-of-fold prediction vector is incomplete for ", model_name,
             call. = FALSE)
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
setcolorder(sample_fold_metrics, c("model", "fold", "scope", "n",
                                   "anti_only_fraction", "roc_auc",
                                   "average_precision", "log_loss",
                                   "brier_score"))

performance <- rbindlist(lapply(names(predictions), function(model_name) {
    prediction <- predictions[[model_name]]
    rows <- list(metric_row(model_name, "pooled", evidence$outcome, prediction))
    for (sample in samples$sample_id) {
        selected <- evidence$sample_id == sample
        rows[[length(rows) + 1L]] <- metric_row(
            model_name, sample, evidence$outcome[selected], prediction[selected]
        )
    }
    sample_rows <- rbindlist(rows[-1L])
    rows[[length(rows) + 1L]] <- data.table(
        model = model_name,
        scope = "sample_macro_mean",
        n = sum(sample_rows$n),
        anti_only_fraction = mean(sample_rows$anti_only_fraction),
        roc_auc = mean(sample_rows$roc_auc),
        average_precision = mean(sample_rows$average_precision),
        log_loss = mean(sample_rows$log_loss),
        brier_score = mean(sample_rows$brier_score)
    )
    sample_fold_rows <- sample_fold_metrics[model == model_name]
    rows[[length(rows) + 1L]] <- data.table(
        model = model_name,
        scope = "sample_fold_macro_mean",
        n = sum(sample_fold_rows$n),
        anti_only_fraction = mean(sample_fold_rows$anti_only_fraction),
        roc_auc = mean(sample_fold_rows$roc_auc),
        average_precision = mean(sample_fold_rows$average_precision),
        log_loss = mean(sample_fold_rows$log_loss),
        brier_score = mean(sample_fold_rows$brier_score)
    )
    rbindlist(rows)
}))
fold_metrics <- rbindlist(fold_metrics, fill = TRUE)
setcolorder(fold_metrics, c("model", "fold", "scope", "n",
                            "anti_only_fraction", "roc_auc",
                            "average_precision", "log_loss", "brier_score"))
fwrite(performance, paste0(output_prefix, "_model_metrics.tsv"), sep = "\t")
fwrite(fold_metrics, paste0(output_prefix, "_fold_metrics.tsv"), sep = "\t")
fwrite(sample_fold_metrics,
       paste0(output_prefix, "_sample_fold_metrics.tsv"), sep = "\t")

evaluation_strata <- list(
    all = rep(TRUE, nrow(evidence)),
    tp73_ge_0 = evidence$tp73_score >= 0,
    tp73_minus1_to_0 = evidence$tp73_score >= -1 & evidence$tp73_score < 0
)
evaluation_strata <- evaluation_strata[vapply(evaluation_strata, any, logical(1))]
stratum_sample_fold_metrics <- rbindlist(lapply(names(predictions),
                                                function(model_name) {
    prediction <- predictions[[model_name]]
    rbindlist(lapply(names(evaluation_strata), function(stratum_name) {
        in_stratum <- evaluation_strata[[stratum_name]]
        rbindlist(lapply(seq_len(values$folds), function(held_out) {
            rbindlist(lapply(samples$sample_id, function(sample) {
                selected <- in_stratum & evidence$fold == held_out &
                            evidence$sample_id == sample
                metric_row(model_name, sample, evidence$outcome[selected],
                           prediction[selected])[
                    , `:=`(stratum = stratum_name, fold = held_out)
                ]
            }))
        }))
    }))
}))
setcolorder(stratum_sample_fold_metrics,
            c("model", "stratum", "fold", "scope", "n",
              "anti_only_fraction", "roc_auc", "average_precision",
              "log_loss", "brier_score"))
stratum_metrics <- rbindlist(lapply(names(predictions), function(model_name) {
    prediction <- predictions[[model_name]]
    rbindlist(lapply(names(evaluation_strata), function(stratum_name) {
        selected <- evaluation_strata[[stratum_name]]
        pooled <- metric_row(model_name, "pooled", evidence$outcome[selected],
                             prediction[selected])[, stratum := stratum_name]
        cells <- stratum_sample_fold_metrics[
            model == model_name & stratum == stratum_name
        ]
        macro <- data.table(
            model = model_name,
            scope = "sample_fold_macro_mean",
            n = sum(cells$n),
            anti_only_fraction = mean(cells$anti_only_fraction),
            roc_auc = mean(cells$roc_auc),
            average_precision = mean(cells$average_precision),
            log_loss = mean(cells$log_loss),
            brier_score = mean(cells$brier_score),
            stratum = stratum_name
        )
        rbindlist(list(pooled, macro), use.names = TRUE)
    }))
}))
setcolorder(stratum_metrics,
            c("model", "stratum", "scope", "n", "anti_only_fraction",
              "roc_auc", "average_precision", "log_loss", "brier_score"))
fwrite(stratum_metrics,
       paste0(output_prefix, "_anchor_stratum_metrics.tsv"), sep = "\t")
fwrite(stratum_sample_fold_metrics,
       paste0(output_prefix, "_anchor_stratum_sample_fold_metrics.tsv"),
       sep = "\t")

primary <- performance[scope == "sample_fold_macro_mean"]
model_order <- names(model_specs)
primary[, model := factor(model, levels = rev(model_order))]
plot_data <- melt(
    primary,
    id.vars = "model",
    measure.vars = c("roc_auc", "average_precision", "log_loss"),
    variable.name = "metric", value.name = "value"
)
plot_data[, metric := factor(
    metric,
    levels = c("roc_auc", "average_precision", "log_loss"),
    labels = c("ROC AUC", "Average precision", "Log loss (lower is better)")
)]
plot <- ggplot(plot_data, aes(x = value, y = model, fill = model)) +
    geom_col(width = 0.72, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.4f", value)), hjust = -0.08, size = 3.5) +
    facet_wrap(~metric, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = c(
        "#6B7280", "#0072B2", "#CC79A7", "#009E73",
        "#E69F00", "#56B4E9", "#D55E00"
    )) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
        title = "Within-sample prediction of TP73-specific CUT&RUN immersion",
        subtitle = paste0(
            "Macro-average across six samples and ", values$folds,
            " chromosome-block folds; discordant anti-p73/control pairs"
        ),
        x = NULL,
        y = NULL,
        caption = paste0(
            if (is.null(context_floor)) {
                "Exact local motif maxima"
            } else {
                paste0("Local motif maxima censored below ", context_floor)
            },
            " within +/-150 bp; exploratory chromosome-1 analysis"
        )
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold"))
ggsave(paste0(output_prefix, "_comparison.png"), plot,
       width = 14, height = 7, dpi = 180, bg = "white")

message("I: wrote model results with ", nrow(anchor_features), " anchors, ",
        nrow(evidence), " discordant observations, and context floor ",
        if (is.null(context_floor)) "none" else context_floor, " to ",
        output_prefix, "_*")
