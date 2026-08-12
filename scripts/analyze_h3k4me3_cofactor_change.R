#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: analyze_h3k4me3_cofactor_change.R (--signal FILE.parquet |",
        "       --change FILE.parquet)",
        "       --tp73-evidence FILE.parquet --cofactor-maxima FILE.parquet",
        "       --thresholds FILE.tsv --output-prefix PATH [options]",
        "",
        "Estimate, separately by experimental series, whether a cofactor class",
        "predicts GFP-referenced H3K4me3 change at TP73 motif anchors. The primary",
        "outcome is log2((H3K4me3 area + alpha)/(input area + alpha)) in TA or DN",
        "minus the matching GFP value. TP73 confirmation and its cofactor",
        "interaction are secondary, descriptive post-treatment analyses.",
        "",
        "Options:",
        "  --signal FILE             Windowed H3K4me3/input Parquet; repeatable",
        "  --change FILE             Precomputed GFP-referenced change; repeatable",
        "  --tp73-evidence FILE      Strict TP73/negative-control evidence Parquet",
        "  --cofactor-maxima FILE    Schema-v4 cofactor maxima; repeatable",
        "  --annotation FILE         Schema-9 tp73_context_anchor; repeatable",
        "  --thresholds FILE         motif_id and positive/recommended threshold TSV",
        "  --window NAME             Signal window (required with --signal)",
        "  --output-prefix PATH      Prefix for result TSV files",
        "  --series ID               Required series; repeat (default: all in signal)",
        "  --negative-references CSV Strict negative limits (default: -1,0)",
        "  --pseudocount VALUE       Integrated-signal pseudocount (default: 1)",
        "  --block-size BP           Genomic SE cluster width (default: 5000000)",
        "  --spline-df N             TP73-score spline degrees (default: 3)",
        "  --minimum-class-fraction F  Per-class minimum (default: 0.005)",
        "  --minimum-class-count N   Per-class minimum (default: 100)",
        "  --minimum-interaction-cell-count N  Four-cell minimum (default: 100)",
        "  --duckdb PATH             DuckDB CLI (default: duckdb)",
        "  --analysis-role TEXT      Development/validation role for provenance",
        "                            (default: unspecified)",
        "  -h, --help                Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) usage()

values <- list(
    signal = character(),
    change = character(),
    cofactor_maxima = character(),
    annotation = character(),
    series = character(),
    negative_references = "-1,0",
    pseudocount = 1,
    block_size = 5000000,
    spline_df = 3L,
    minimum_class_fraction = 0.005,
    minimum_class_count = 100L,
    minimum_interaction_cell_count = 100L,
    duckdb = "duckdb",
    analysis_role = "unspecified"
)
value_options <- c(
    "--signal", "--change", "--tp73-evidence", "--cofactor-maxima", "--thresholds",
    "--annotation",
    "--window", "--output-prefix", "--series", "--negative-references",
    "--pseudocount", "--block-size", "--spline-df",
    "--minimum-class-fraction", "--minimum-class-count",
    "--minimum-interaction-cell-count", "--duckdb", "--analysis-role"
)
index <- 1L
while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (!option %in% value_options) {
        writeLines(paste("E: unknown argument:", option), con = stderr())
        usage(2L)
    }
    index <- index + 1L
    if (index > length(arguments)) {
        writeLines(paste("E:", option, "requires a value"), con = stderr())
        usage(2L)
    }
    key <- gsub("-", "_", substring(option, 3L), fixed = TRUE)
    if (option %in% c(
        "--signal", "--change", "--cofactor-maxima", "--annotation", "--series"
    )) {
        values[[key]] <- c(values[[key]], arguments[[index]])
    } else {
        values[[key]] <- arguments[[index]]
    }
    index <- index + 1L
}

required <- c("tp73_evidence", "thresholds", "output_prefix")
missing <- required[!vapply(required, function(name) {
    !is.null(values[[name]]) && nzchar(values[[name]])
}, logical(1))]
if (length(missing) > 0L) {
    writeLines(paste("E: missing options:", paste(missing, collapse = ", ")),
               con = stderr())
    usage(2L)
}
if (length(values$cofactor_maxima) == 0L) {
    stop("provide at least one --cofactor-maxima", call. = FALSE)
}
has_signal <- length(values$signal) > 0L
has_change <- length(values$change) > 0L
if (has_signal == has_change) {
    stop("provide exactly one of --signal and --change", call. = FALSE)
}
if (has_signal && (is.null(values$window) || !nzchar(values$window))) {
    stop("--signal requires --window", call. = FALSE)
}

suppressPackageStartupMessages(library(data.table))

vector_input_paths <- c(
    if (has_signal) "signal" else "change", "cofactor_maxima",
    if (length(values$annotation) > 0L) "annotation"
)
for (name in vector_input_paths) {
    values[[name]] <- normalizePath(values[[name]], mustWork = TRUE)
    if (anyDuplicated(values[[name]])) {
        stop("duplicate --", gsub("_", "-", name), " input", call. = FALSE)
    }
}
for (name in c("tp73_evidence", "thresholds")) {
    values[[name]] <- normalizePath(values[[name]], mustWork = TRUE)
}
values$pseudocount <- suppressWarnings(as.numeric(values$pseudocount))
values$block_size <- suppressWarnings(as.numeric(values$block_size))
values$spline_df <- suppressWarnings(as.integer(values$spline_df))
values$minimum_class_fraction <- suppressWarnings(
    as.numeric(values$minimum_class_fraction)
)
values$minimum_class_count <- suppressWarnings(as.integer(values$minimum_class_count))
values$minimum_interaction_cell_count <- suppressWarnings(
    as.integer(values$minimum_interaction_cell_count)
)
if (!is.finite(values$pseudocount) || values$pseudocount <= 0 ||
    !is.finite(values$block_size) || values$block_size <= 300 ||
    is.na(values$spline_df) || values$spline_df < 1L ||
    !is.finite(values$minimum_class_fraction) ||
    values$minimum_class_fraction <= 0 || values$minimum_class_fraction >= 0.5 ||
    is.na(values$minimum_class_count) || values$minimum_class_count < 1L ||
    is.na(values$minimum_interaction_cell_count) ||
    values$minimum_interaction_cell_count < 1L) {
    stop("invalid numeric analysis option", call. = FALSE)
}
negative_references <- sort(unique(suppressWarnings(as.numeric(strsplit(
    values$negative_references, ",", fixed = TRUE
)[[1L]]))))
if (length(negative_references) == 0L || any(!is.finite(negative_references))) {
    stop("--negative-references must contain finite numbers", call. = FALSE)
}
if (any(!nzchar(values$series)) || anyDuplicated(values$series)) {
    stop("--series values must be unique and non-empty", call. = FALSE)
}
if (!nzchar(values$analysis_role) || grepl("[\r\n]", values$analysis_role)) {
    stop("--analysis-role must be non-empty and single-line", call. = FALSE)
}
duckdb_path <- Sys.which(values$duckdb)
if (!nzchar(duckdb_path)) stop("DuckDB CLI not found: ", values$duckdb, call. = FALSE)
dir.create(dirname(values$output_prefix), recursive = TRUE, showWarnings = FALSE)

sql_string <- function(value) paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
sql_identifier <- function(value) paste0('"', gsub('"', '""', value, fixed = TRUE), '"')
sql_parquet_paths <- function(paths) {
    quoted <- vapply(paths, sql_string, character(1))
    if (length(quoted) == 1L) quoted else paste0("[", paste(quoted, collapse = ","), "]")
}

duckdb_fread <- function(query) {
    sql_path <- tempfile("h3k4me3-change-", fileext = ".sql")
    on.exit(unlink(sql_path), add = TRUE)
    writeLines(query, sql_path, useBytes = TRUE)
    command <- paste(
        shQuote(duckdb_path), "-light-mode -csv -nullvalue NA :memory: <",
        shQuote(sql_path)
    )
    fread(cmd = command, na.strings = "NA", showProgress = FALSE)
}

thresholds <- fread(values$thresholds, sep = "\t", na.strings = "NA")
if (!"motif_id" %in% names(thresholds)) {
    stop("threshold table lacks motif_id", call. = FALSE)
}
threshold_column <- if ("positive_threshold" %in% names(thresholds)) {
    "positive_threshold"
} else if ("recommended_threshold" %in% names(thresholds)) {
    "recommended_threshold"
} else {
    stop("threshold table lacks positive_threshold/recommended_threshold", call. = FALSE)
}
for (optional in c("factor_name", "positive_threshold_source", "selection_semantics")) {
    if (!optional %in% names(thresholds)) thresholds[, (optional) := NA_character_]
}
thresholds <- thresholds[, .(
    motif_id = as.character(motif_id),
    factor_name = as.character(factor_name),
    positive_threshold = as.numeric(get(threshold_column)),
    positive_threshold_source = as.character(positive_threshold_source),
    selection_semantics = as.character(selection_semantics)
)]
if (nrow(thresholds) == 0L || any(!nzchar(thresholds$motif_id)) ||
    any(!is.finite(thresholds$positive_threshold)) || anyDuplicated(thresholds$motif_id)) {
    stop("threshold table is empty, duplicated, or invalid", call. = FALSE)
}
setorder(thresholds, motif_id)

if (has_signal) {
    input_mode <- "windowed_signal"
    signal <- duckdb_fread(paste0(
        "SELECT CAST(chrom AS VARCHAR) AS chrom, anchor_start, anchor_end, ",
        "anchor_score, series_id, cell_line, condition, replicate, window_name, ",
        "effective_window_bp, h3k4me3_area, h3k4me3_positive_bp, input_area, ",
        "input_positive_bp FROM read_parquet(", sql_parquet_paths(values$signal), ") ",
        "WHERE window_name = ", sql_string(values$window), " ORDER BY series_id, ",
        "condition, chrom, anchor_start, anchor_end;"
    ))
    required_signal <- c(
        "chrom", "anchor_start", "anchor_end", "anchor_score", "series_id",
        "cell_line", "condition", "replicate", "window_name",
        "effective_window_bp", "h3k4me3_area", "h3k4me3_positive_bp",
        "input_area", "input_positive_bp"
    )
    if (nrow(signal) == 0L ||
        length(setdiff(required_signal, names(signal))) > 0L) {
        stop("selected H3K4me3 signal window is empty or incomplete", call. = FALSE)
    }
    available_series <- sort(unique(signal$series_id))
    analysis_series <- if (length(values$series) > 0L) {
        values$series
    } else available_series
    if (!setequal(analysis_series, intersect(analysis_series, available_series))) {
        stop("a required --series is absent from signal", call. = FALSE)
    }
    signal <- signal[series_id %in% analysis_series]
    if (!setequal(unique(signal$condition), c("GFP", "TA", "DN"))) {
        stop("signal does not contain the GFP/TA/DN factorial", call. = FALSE)
    }
    if (signal[, anyDuplicated(paste(
        chrom, anchor_start, anchor_end, series_id, condition
    ))] != 0L) {
        stop("signal has duplicate anchor/series/condition rows", call. = FALSE)
    }
    if (signal[, any(
        !is.finite(anchor_score) | !is.finite(h3k4me3_area) |
        !is.finite(input_area) | h3k4me3_area < 0 | input_area < 0 |
        effective_window_bp <= 0
    )]) stop("signal contains invalid values", call. = FALSE)

    signal[, normalized_mark := log2(
        (h3k4me3_area + values$pseudocount) /
        (input_area + values$pseudocount)
    )]
    signal[, mark_above_input := h3k4me3_area > input_area]
    signal[, any_observed_signal := h3k4me3_area > 0 | input_area > 0]
    wide <- dcast(
        signal,
        chrom + anchor_start + anchor_end + anchor_score + series_id + cell_line +
            replicate + effective_window_bp ~ condition,
        value.var = c(
            "normalized_mark", "mark_above_input", "any_observed_signal",
            "h3k4me3_area", "input_area"
        )
    )
    expected_wide <- c(
        "normalized_mark_GFP", "normalized_mark_TA", "normalized_mark_DN",
        "mark_above_input_GFP", "mark_above_input_TA", "mark_above_input_DN"
    )
    if (length(setdiff(expected_wide, names(wide))) > 0L ||
        nrow(wide) * 3L != nrow(signal)) {
        stop("failed to reshape the condition-factorial signal", call. = FALSE)
    }
    wide[, delta_TA := normalized_mark_TA - normalized_mark_GFP]
    wide[, delta_DN := normalized_mark_DN - normalized_mark_GFP]
} else {
    input_mode <- "precomputed_change"
    window_filter <- if (!is.null(values$window) && nzchar(values$window)) {
        paste0(" WHERE window_name = ", sql_string(values$window))
    } else ""
    wide <- duckdb_fread(paste0(
        "SELECT CAST(chrom AS VARCHAR) AS chrom, anchor_start, anchor_end, ",
        "anchor_score, series_id, cell_line, replicate, window_name, ",
        "effective_window_bp, ",
        "gfp_log2_h3k4me3_input_ratio AS normalized_mark_GFP, ",
        "ta_log2_h3k4me3_input_ratio AS normalized_mark_TA, ",
        "dn_log2_h3k4me3_input_ratio AS normalized_mark_DN, ",
        "gfp_h3k4me3_area AS h3k4me3_area_GFP, ",
        "ta_h3k4me3_area AS h3k4me3_area_TA, ",
        "dn_h3k4me3_area AS h3k4me3_area_DN, ",
        "gfp_input_area AS input_area_GFP, ta_input_area AS input_area_TA, ",
        "dn_input_area AS input_area_DN, delta_ta_vs_gfp AS delta_TA, ",
        "delta_dn_vs_gfp AS delta_DN ",
        "FROM read_parquet(", sql_parquet_paths(values$change), ")", window_filter,
        " ORDER BY series_id, chrom, anchor_start, anchor_end;"
    ))
    required_change <- c(
        "chrom", "anchor_start", "anchor_end", "anchor_score", "series_id",
        "cell_line", "replicate", "window_name", "effective_window_bp",
        "normalized_mark_GFP", "normalized_mark_TA", "normalized_mark_DN",
        "h3k4me3_area_GFP", "h3k4me3_area_TA", "h3k4me3_area_DN",
        "input_area_GFP", "input_area_TA", "input_area_DN", "delta_TA", "delta_DN"
    )
    if (nrow(wide) == 0L || length(setdiff(required_change, names(wide))) > 0L ||
        uniqueN(wide$window_name) != 1L ||
        wide[, anyDuplicated(paste(
            chrom, anchor_start, anchor_end, series_id
        ))] != 0L) {
        stop("precomputed H3K4me3 change is empty or violates its key", call. = FALSE)
    }
    values$window <- unique(wide$window_name)
    available_series <- sort(unique(wide$series_id))
    analysis_series <- if (length(values$series) > 0L) {
        values$series
    } else available_series
    if (!setequal(analysis_series, intersect(analysis_series, available_series))) {
        stop("a required --series is absent from change input", call. = FALSE)
    }
    wide <- wide[series_id %in% analysis_series]
    if (wide[, any(
        !is.finite(anchor_score) | !is.finite(delta_TA) | !is.finite(delta_DN) |
        h3k4me3_area_GFP < 0 | h3k4me3_area_TA < 0 | h3k4me3_area_DN < 0 |
        input_area_GFP < 0 | input_area_TA < 0 | input_area_DN < 0 |
        effective_window_bp <= 0
    )]) stop("precomputed change contains invalid values", call. = FALSE)
    wide[, mark_above_input_GFP := h3k4me3_area_GFP > input_area_GFP]
    wide[, mark_above_input_TA := h3k4me3_area_TA > input_area_TA]
    wide[, mark_above_input_DN := h3k4me3_area_DN > input_area_DN]
    wide[, any_observed_signal_GFP :=
             h3k4me3_area_GFP > 0 | input_area_GFP > 0]
    wide[, any_observed_signal_TA :=
             h3k4me3_area_TA > 0 | input_area_TA > 0]
    wide[, any_observed_signal_DN :=
             h3k4me3_area_DN > 0 | input_area_DN > 0]
}
wide[, occurrence_TA := fifelse(
    mark_above_input_TA & !mark_above_input_GFP, "gained",
    fifelse(!mark_above_input_TA & mark_above_input_GFP, "lost",
            fifelse(mark_above_input_TA, "shared_positive", "shared_nonpositive"))
)]
wide[, occurrence_DN := fifelse(
    mark_above_input_DN & !mark_above_input_GFP, "gained",
    fifelse(!mark_above_input_DN & mark_above_input_GFP, "lost",
            fifelse(mark_above_input_DN, "shared_positive", "shared_nonpositive"))
)]
wide[, h3k4me3_informative :=
         h3k4me3_area_GFP > 0 | h3k4me3_area_TA > 0 | h3k4me3_area_DN > 0]

evidence_description <- duckdb_fread(paste0(
    "DESCRIBE SELECT * FROM read_parquet(", sql_string(values$tp73_evidence), ");"
))
evidence_columns <- as.character(evidence_description$column_name)
base_evidence <- c("chrom", "anchor_start", "anchor_end", "anchor_score")
binding_columns <- unlist(lapply(analysis_series, function(series_id) {
    unlist(lapply(c("GFP", "TA", "DN"), function(condition) c(
        paste0("supported_tp73_", series_id, "_", condition),
        paste0("depth_tp73_", series_id, "_", condition),
        paste0("supported_negative_control_", series_id, "_", condition),
        paste0("depth_negative_control_", series_id, "_", condition)
    )))
}))
missing_evidence <- setdiff(c(base_evidence, binding_columns), evidence_columns)
if (length(missing_evidence) > 0L) {
    stop("TP73 evidence lacks columns: ", paste(missing_evidence, collapse = ", "),
         call. = FALSE)
}
evidence_select <- c(
    "CAST(chrom AS VARCHAR) AS chrom", "anchor_start", "anchor_end",
    vapply(binding_columns, sql_identifier, character(1))
)
evidence <- duckdb_fread(paste0(
    "SELECT ", paste(evidence_select, collapse = ", "), " FROM read_parquet(",
    sql_string(values$tp73_evidence), ") ORDER BY chrom, anchor_start, anchor_end;"
))
if (evidence[, anyDuplicated(paste(chrom, anchor_start, anchor_end))] != 0L) {
    stop("TP73 evidence has duplicate anchor spans", call. = FALSE)
}

binding_rows <- lapply(analysis_series, function(series_id) {
    result <- evidence[, .(chrom, anchor_start, anchor_end)]
    result[, series_id := series_id]
    for (condition in c("GFP", "TA", "DN")) {
        tp73_support <- as.logical(evidence[[paste0(
            "supported_tp73_", series_id, "_", condition
        )]])
        negative_support <- as.logical(evidence[[paste0(
            "supported_negative_control_", series_id, "_", condition
        )]])
        tp73_depth <- as.numeric(evidence[[paste0(
            "depth_tp73_", series_id, "_", condition
        )]])
        negative_depth <- as.numeric(evidence[[paste0(
            "depth_negative_control_", series_id, "_", condition
        )]])
        if (anyNA(tp73_support) || anyNA(negative_support) ||
            any(!is.finite(tp73_depth)) || any(!is.finite(negative_depth))) {
            stop("invalid TP73 evidence for ", series_id, " ", condition,
                 call. = FALSE)
        }
        result[, (paste0("confirmed_", condition)) :=
                   tp73_support & !negative_support]
        result[, (paste0("tp73_depth_", condition)) := tp73_depth]
        result[, (paste0("negative_depth_", condition)) := negative_depth]
    }
    result[, binding_state_TA := fifelse(
        confirmed_TA & !confirmed_GFP, "gained_in_TA",
        fifelse(confirmed_TA & confirmed_GFP, "shared_TA_GFP",
                fifelse(!confirmed_TA & confirmed_GFP, "lost_in_TA", "unsupported"))
    )]
    result[, binding_state_DN := fifelse(
        confirmed_DN & !confirmed_GFP, "gained_in_DN",
        fifelse(confirmed_DN & confirmed_GFP, "shared_DN_GFP",
                fifelse(!confirmed_DN & confirmed_GFP, "lost_in_DN", "unsupported"))
    )]
    result
})
binding <- rbindlist(binding_rows)
setkey(binding, chrom, anchor_start, anchor_end, series_id)
setkey(wide, chrom, anchor_start, anchor_end, series_id)
wide <- binding[wide]
if (nrow(wide) == 0L || wide[, anyNA(confirmed_GFP)]) {
    stop("TP73 evidence did not join completely to H3K4me3 signal", call. = FALSE)
}

annotation_available <- length(values$annotation) > 0L
if (annotation_available) {
    annotation <- duckdb_fread(paste0(
        "SELECT CAST(annotation.chrom AS VARCHAR) AS chrom, ",
        "annotation.start AS anchor_start, annotation.\"end\" AS anchor_end, ",
        "min(primary_genomic_context) AS ",
        "primary_genomic_context, count(DISTINCT primary_genomic_context) ",
        "AS context_values, bool_or(strict_intergenic) AS strict_intergenic, ",
        "count(DISTINCT strict_intergenic) AS intergenic_values, ",
        "bool_or(overlaps_any_promoter) AS overlaps_any_promoter, ",
        "bool_or(overlaps_any_downstream_region) AS ",
        "overlaps_any_downstream_region, ",
        "min(gene_relation_class) AS gene_relation_class, ",
        "count(DISTINCT gene_relation_class) AS gene_relation_values, ",
        "bool_or(in_any_transcript) AS in_any_transcript, ",
        "bool_or(in_any_exon) AS in_any_exon, bool_or(in_any_cds) AS in_any_cds, ",
        "min(nearest_tss_distance_bp) AS nearest_tss_distance_bp, ",
        "count(DISTINCT nearest_tss_distance_bp) AS nearest_tss_distance_values, ",
        "min(nearest_tss_genomic_distance_bp) AS nearest_tss_genomic_distance_bp, ",
        "count(DISTINCT nearest_tss_genomic_distance_bp) AS tss_genomic_values, ",
        "bool_or(nearest_tss_has_mixed_strands) AS nearest_tss_has_mixed_strands, ",
        "count(DISTINCT nearest_tss_has_mixed_strands) AS tss_mixed_values, ",
        "min(nearest_tss_relation) AS nearest_tss_relation, ",
        "count(DISTINCT nearest_tss_relation) AS tss_relation_values, ",
        "min(nearest_cds_distance_bp) AS nearest_cds_distance_bp, ",
        "count(DISTINCT nearest_cds_distance_bp) AS nearest_cds_distance_values, ",
        "min(nearest_cds_genomic_distance_bp) AS nearest_cds_genomic_distance_bp, ",
        "count(DISTINCT nearest_cds_genomic_distance_bp) AS cds_genomic_values, ",
        "bool_or(nearest_cds_has_mixed_strands) AS nearest_cds_has_mixed_strands, ",
        "count(DISTINCT nearest_cds_has_mixed_strands) AS cds_mixed_values, ",
        "min(nearest_cds_relation) AS nearest_cds_relation, ",
        "count(DISTINCT nearest_cds_relation) AS cds_relation_values ",
        "FROM read_parquet(", sql_parquet_paths(values$annotation), ", ",
        "hive_partitioning=true) ",
        "AS annotation GROUP BY CAST(annotation.chrom AS VARCHAR), ",
        "annotation.start, annotation.\"end\" ORDER BY chrom, anchor_start, ",
        "anchor_end;"
    ))
    if (nrow(annotation) == 0L || annotation[, any(
        context_values != 1L | intergenic_values != 1L |
        gene_relation_values != 1L |
        nearest_tss_distance_values > 1L | nearest_cds_distance_values > 1L |
        tss_genomic_values > 1L | cds_genomic_values > 1L |
        tss_mixed_values > 1L | cds_mixed_values > 1L |
        tss_relation_values > 1L | cds_relation_values > 1L
    )] || annotation[, anyDuplicated(paste(
        chrom, anchor_start, anchor_end
    ))] != 0L) {
        stop("schema-9 annotation is empty or differs across anchor orientations",
             call. = FALSE)
    }
    annotation[, analysis_context_class := fcase(
        overlaps_any_promoter, "promoter",
        overlaps_any_downstream_region, "downstream",
        strict_intergenic, "strict_intergenic",
        primary_genomic_context %in% c("cds", "cds_boundary"), "cds",
        primary_genomic_context %in% c("exonic_non_cds", "exon_boundary"),
            "exonic_non_cds",
        primary_genomic_context %in% c("intron", "intron_boundary"), "intron",
        default = "other_transcribed"
    )]
    annotation[, expected_gene_relation_class := fcase(
        overlaps_any_promoter, "promoter",
        overlaps_any_downstream_region, "downstream",
        in_any_transcript, "gene_body",
        default = "intergenic"
    )]
    if (annotation[, any(
        gene_relation_class != expected_gene_relation_class |
        strict_intergenic != (gene_relation_class == "intergenic")
    )]) {
        stop("schema-9 gene-relation precedence is inconsistent", call. = FALSE)
    }
    annotation[, expected_gene_relation_class := NULL]
    annotation[, nearest_tss_direction_class := fcase(
        is.na(nearest_tss_distance_bp), "unavailable",
        nearest_tss_has_mixed_strands %in% TRUE, "mixed_strand_tie",
        nearest_tss_relation %in% c("spans_tss", "coincident_center"),
            "overlap_or_coincident",
        nearest_tss_distance_bp < 0, "upstream",
        nearest_tss_distance_bp > 0, "downstream",
        default = "coincident"
    )]
    annotation[, nearest_cds_direction_class := fcase(
        is.na(nearest_cds_distance_bp), "unavailable",
        nearest_cds_has_mixed_strands %in% TRUE, "mixed_strand_tie",
        nearest_cds_relation %in% "overlaps_cds", "overlap",
        nearest_cds_distance_bp < 0, "upstream",
        nearest_cds_distance_bp > 0, "downstream",
        default = "coincident"
    )]
    annotation[, log_nearest_tss_genomic_distance :=
                   log1p(nearest_tss_genomic_distance_bp)]
    annotation[, log_nearest_cds_genomic_distance :=
                   log1p(nearest_cds_genomic_distance_bp)]
    annotation[, c(
        "context_values", "intergenic_values", "gene_relation_values",
        "nearest_tss_distance_values",
        "nearest_cds_distance_values", "tss_genomic_values",
        "cds_genomic_values", "tss_mixed_values", "cds_mixed_values",
        "tss_relation_values", "cds_relation_values"
    ) := NULL]
    setkey(annotation, chrom, anchor_start, anchor_end)
    setkey(wide, chrom, anchor_start, anchor_end)
    wide <- annotation[wide]
    if (nrow(wide) == 0L || wide[, anyNA(primary_genomic_context)]) {
        stop("schema-9 annotation did not join every H3K4me3 anchor", call. = FALSE)
    }
} else {
    wide[, `:=`(
        primary_genomic_context = "not_assessed",
        strict_intergenic = NA,
        overlaps_any_promoter = NA,
        overlaps_any_downstream_region = NA,
        in_any_transcript = NA,
        in_any_exon = NA,
        in_any_cds = NA,
        nearest_tss_distance_bp = NA_real_,
        nearest_cds_distance_bp = NA_real_,
        nearest_tss_genomic_distance_bp = NA_real_,
        nearest_cds_genomic_distance_bp = NA_real_,
        nearest_tss_has_mixed_strands = NA,
        nearest_cds_has_mixed_strands = NA,
        nearest_tss_relation = NA_character_,
        nearest_cds_relation = NA_character_,
        analysis_context_class = "not_assessed",
        gene_relation_class = "not_assessed",
        nearest_tss_direction_class = "not_assessed",
        nearest_cds_direction_class = "not_assessed",
        log_nearest_tss_genomic_distance = NA_real_,
        log_nearest_cds_genomic_distance = NA_real_
    )]
}

motif_list <- paste(vapply(thresholds$motif_id, sql_string, character(1)), collapse = ",")
maxima <- duckdb_fread(paste0(
    "SELECT CAST(chrom AS VARCHAR) AS chrom, anchor_start, anchor_end, motif_id, ",
    "context_score, source_score_floor, context_flank_bp, context_distance_metric ",
    "FROM read_parquet(", sql_parquet_paths(values$cofactor_maxima),
    ") WHERE motif_id IN (",
    motif_list, ") ORDER BY motif_id, chrom, anchor_start, anchor_end;"
))
if (nrow(maxima) == 0L || !setequal(unique(maxima$motif_id), thresholds$motif_id)) {
    stop("cofactor maxima and threshold motifs differ", call. = FALSE)
}
anchor_count <- uniqueN(wide, by = c("chrom", "anchor_start", "anchor_end"))
if (nrow(maxima) != anchor_count * nrow(thresholds) ||
    maxima[, anyDuplicated(paste(motif_id, chrom, anchor_start, anchor_end))] != 0L) {
    stop("cofactor maxima are not rectangular over the signal anchors", call. = FALSE)
}
if (maxima[, any(context_distance_metric != "signed_interval_edge_distance")]) {
    stop("cofactor maxima do not use schema-v4 interval geometry", call. = FALSE)
}
maxima_provenance <- maxima[, .(
    source_floor_values = uniqueN(source_score_floor),
    source_score_floor = first(source_score_floor),
    context_flank_values = uniqueN(context_flank_bp),
    context_flank_bp = first(context_flank_bp),
    distance_metric_values = uniqueN(context_distance_metric)
), by = motif_id]
if (maxima_provenance[, any(
    source_floor_values != 1L | context_flank_values != 1L |
    distance_metric_values != 1L | !is.finite(source_score_floor) |
    is.na(context_flank_bp) | context_flank_bp <= 0
)]) {
    stop("cofactor maxima have inconsistent or invalid per-motif provenance",
         call. = FALSE)
}
threshold_provenance <- merge(
    thresholds[, .(motif_id, positive_threshold)],
    maxima_provenance[, .(motif_id, source_score_floor)],
    by = "motif_id", all.x = TRUE, sort = FALSE
)
if (threshold_provenance[, any(
    is.na(source_score_floor) | positive_threshold < source_score_floor
)]) {
    stop("a positive threshold is below its cofactor scan floor", call. = FALSE)
}

clustered_lm <- function(data, formula, contrasts) {
    fit <- tryCatch(lm(formula, data = data, x = TRUE, y = TRUE),
                    error = function(condition) NULL)
    empty <- function(status, note) data.table(
        estimate = NA_real_, standard_error = NA_real_,
        confidence_interval_95_lower = NA_real_,
        confidence_interval_95_upper = NA_real_, p_value = NA_real_,
        genomic_blocks = uniqueN(data$genomic_block),
        observations = nrow(data), evaluation_status = status,
        evaluation_note = note
    )
    if (is.null(fit) || any(!is.finite(coef(fit)))) {
        return(lapply(contrasts, function(x) empty("not_estimable", "linear model failed")))
    }
    x <- fit$x
    inverse <- tryCatch(solve(crossprod(x)), error = function(condition) NULL)
    if (is.null(inverse) || any(!is.finite(inverse))) {
        return(lapply(contrasts, function(x) empty("not_estimable", "singular model matrix")))
    }
    cluster_scores <- rowsum(x * residuals(fit), data$genomic_block, reorder = FALSE)
    clusters <- nrow(cluster_scores)
    observations <- nrow(x)
    parameters <- ncol(x)
    if (clusters < 2L || observations <= parameters) {
        return(lapply(contrasts, function(x) empty(
            "not_estimable", "too few genomic blocks for clustered uncertainty"
        )))
    }
    correction <- (clusters / (clusters - 1)) *
        ((observations - 1) / (observations - parameters))
    covariance <- correction * inverse %*% crossprod(cluster_scores) %*% inverse
    lapply(contrasts, function(contrast) {
        vector <- setNames(rep(0, length(coef(fit))), names(coef(fit)))
        if (!all(names(contrast) %in% names(vector))) {
            return(empty("not_estimable", "requested coefficient is absent"))
        }
        vector[names(contrast)] <- contrast
        estimate <- sum(vector * coef(fit))
        standard_error <- sqrt(drop(t(vector) %*% covariance %*% vector))
        if (!is.finite(standard_error) || standard_error <= 0) {
            return(empty("not_estimable", "clustered standard error is invalid"))
        }
        z <- estimate / standard_error
        data.table(
            estimate = estimate, standard_error = standard_error,
            confidence_interval_95_lower = estimate - 1.96 * standard_error,
            confidence_interval_95_upper = estimate + 1.96 * standard_error,
            p_value = 2 * pnorm(abs(z), lower.tail = FALSE),
            genomic_blocks = clusters, observations = observations,
            evaluation_status = "ok",
            evaluation_note = paste("SE clusters", values$block_size, "bp genomic blocks")
        )
    })
}

status_result <- function(status, note, observations = 0L, blocks = 0L) data.table(
    estimate = NA_real_, standard_error = NA_real_,
    confidence_interval_95_lower = NA_real_, confidence_interval_95_upper = NA_real_,
    p_value = NA_real_, genomic_blocks = blocks, observations = observations,
    evaluation_status = status, evaluation_note = note
)

model_formula <- function(data, interaction = FALSE,
                          adjust_context_class = annotation_available,
                          exposure = "cofactor_present") {
    terms <- paste0(
        "splines::ns(anchor_score, df = ", values$spline_df, ")"
    )
    if (uniqueN(data$chrom) > 1L) {
        terms <- c(terms, "factor(chrom)")
    }
    if (annotation_available && adjust_context_class &&
        uniqueN(data$analysis_context_class) > 1L) {
        terms <- c(terms, "factor(analysis_context_class)")
    }
    if (annotation_available &&
        all(is.finite(data$log_nearest_tss_genomic_distance)) &&
        uniqueN(data$log_nearest_tss_genomic_distance) > values$spline_df) {
        terms <- c(terms, paste0(
            "splines::ns(log_nearest_tss_genomic_distance, df = ",
            values$spline_df, ")"
        ))
    }
    if (annotation_available &&
        all(is.finite(data$log_nearest_cds_genomic_distance)) &&
        uniqueN(data$log_nearest_cds_genomic_distance) > values$spline_df) {
        terms <- c(terms, paste0(
            "splines::ns(log_nearest_cds_genomic_distance, df = ",
            values$spline_df, ")"
        ))
    }
    if (annotation_available && uniqueN(data$nearest_tss_direction_class) > 1L) {
        terms <- c(terms, "factor(nearest_tss_direction_class)")
    }
    if (annotation_available && uniqueN(data$nearest_cds_direction_class) > 1L) {
        terms <- c(terms, "factor(nearest_cds_direction_class)")
    }
    exposure <- if (interaction) {
        "cofactor_present * tp73_confirmed"
    } else exposure
    as.formula(paste("outcome ~", paste(c(terms, exposure), collapse = " + ")))
}

intensity_rows <- list()
interaction_rows <- list()
state_rows <- list()
occurrence_rows <- list()
context_intensity_rows <- list()
gene_relation_intensity_rows <- list()
score_effect_rows <- list()
intensity_index <- interaction_index <- state_index <- occurrence_index <- 1L
context_intensity_index <- gene_relation_intensity_index <- score_effect_index <- 1L

for (motif_index in seq_len(nrow(thresholds))) {
    threshold <- thresholds[motif_index]
    motif <- maxima[motif_id == threshold$motif_id]
    setkey(motif, chrom, anchor_start, anchor_end)
    joined <- motif[wide]
    if (nrow(joined) != nrow(wide)) {
        stop("cofactor maxima join changed anchor-series cardinality for ",
             threshold$motif_id, call. = FALSE)
    }
    motif_source_floor <- unique(motif$source_score_floor)
    if (length(motif_source_floor) != 1L || !is.finite(motif_source_floor)) {
        stop("cofactor maxima have no unique source floor for ",
             threshold$motif_id, call. = FALSE)
    }
    for (negative_reference in negative_references) {
        overlap_invalid <- threshold$positive_threshold < negative_reference
        negative_reference_observable <- motif_source_floor <= negative_reference
        joined[, cofactor_positive := !is.na(context_score) &
                   context_score >= threshold$positive_threshold]
        joined[, cofactor_negative := negative_reference_observable &
                   (is.na(context_score) | context_score < negative_reference)]
        joined[, cofactor_intermediate :=
                   !cofactor_positive & !cofactor_negative]
        for (isoform in c("TA", "DN")) {
            outcome_column <- paste0("delta_", isoform)
            confirmed_column <- paste0("confirmed_", isoform)
            state_column <- paste0("binding_state_", isoform)
            occurrence_column <- paste0("occurrence_", isoform)
            for (series_id in analysis_series) {
                current_series <- series_id
                all_current <- joined[series_id == current_series]
                current <- all_current[h3k4me3_informative == TRUE]
                eligible <- current$cofactor_positive | current$cofactor_negative
                positive_count <- sum(current$cofactor_positive)
                negative_count <- sum(current$cofactor_negative)
                eligible_count <- sum(eligible)
                class_supported <- negative_reference_observable &&
                    !overlap_invalid && eligible_count > 0L &&
                    positive_count >= values$minimum_class_count &&
                    negative_count >= values$minimum_class_count &&
                    positive_count / eligible_count >= values$minimum_class_fraction &&
                    negative_count / eligible_count >= values$minimum_class_fraction
                comparison_status <- if (!negative_reference_observable) {
                    "negative_reference_below_source_floor"
                } else if (overlap_invalid) {
                    "overlapping_score_classes"
                } else if (!class_supported) {
                    "underpowered_anchor_class"
                } else {
                    "ok"
                }
                comparison_note <- if (!negative_reference_observable) {
                    paste("cofactor source floor", motif_source_floor,
                          "is above the requested negative reference")
                } else if (overlap_invalid) {
                    "positive threshold is below the negative reference"
                } else if (!class_supported) {
                    "positive/negative cofactor class is below the declared support"
                } else {
                    "positive and negative cofactor classes are disjoint and supported"
                }
                model_data <- current[eligible, .(
                    chrom, anchor_start, anchor_score,
                    outcome = get(outcome_column),
                    cofactor_present = as.integer(cofactor_positive),
                    tp73_confirmed = as.integer(get(confirmed_column)),
                    binding_state = get(state_column),
                    occurrence = get(occurrence_column),
                    analysis_context_class,
                    log_nearest_tss_genomic_distance,
                    log_nearest_cds_genomic_distance,
                    nearest_tss_direction_class,
                    nearest_cds_direction_class
                )]
                model_data[, genomic_block := paste(
                    chrom, floor(anchor_start / values$block_size), sep = ":"
                )]
                if (!negative_reference_observable) {
                    total <- status_result(
                        "negative_reference_below_source_floor",
                        paste("cofactor source floor", motif_source_floor,
                              "is above the requested negative reference")
                    )
                } else if (overlap_invalid) {
                    total <- status_result(
                        "overlapping_score_classes",
                        "positive threshold is below the negative reference"
                    )
                } else if (!class_supported) {
                    total <- status_result(
                        "underpowered_anchor_class",
                        "positive/negative cofactor class is below the declared support"
                    )
                } else {
                    total <- clustered_lm(
                        model_data,
                        model_formula(model_data),
                        list(c(cofactor_present = 1))
                    )[[1L]]
                }
                intensity_rows[[intensity_index]] <- cbind(data.table(
                    motif_id = threshold$motif_id,
                    factor_name = threshold$factor_name,
                    positive_threshold = threshold$positive_threshold,
                    source_score_floor = motif_source_floor,
                    negative_reference_threshold = negative_reference,
                    negative_reference_observable = negative_reference_observable,
                    isoform = isoform, series_id = series_id,
                    anchors_positive = positive_count,
                    anchors_negative = negative_count,
                    anchors_intermediate = sum(current$cofactor_intermediate),
                    anchors_uninformative_h3k4me3 =
                        sum(!all_current$h3k4me3_informative),
                    mean_change_positive = if (positive_count) {
                        mean(current[[outcome_column]][current$cofactor_positive])
                    } else NA_real_,
                    mean_change_negative = if (negative_count) {
                        mean(current[[outcome_column]][current$cofactor_negative])
                    } else NA_real_,
                    raw_change_difference = if (positive_count && negative_count) {
                        mean(current[[outcome_column]][current$cofactor_positive]) -
                        mean(current[[outcome_column]][current$cofactor_negative])
                    } else NA_real_
                ), total)
                intensity_index <- intensity_index + 1L

                if (annotation_available) {
                    for (context_class in sort(unique(
                        current$analysis_context_class
                    ))) {
                        context <- current[
                            analysis_context_class == context_class
                        ]
                        context_eligible <- context$cofactor_positive |
                            context$cofactor_negative
                        context_positive <- sum(context$cofactor_positive)
                        context_negative <- sum(context$cofactor_negative)
                        context_eligible_count <- sum(context_eligible)
                        context_supported <- negative_reference_observable &&
                            !overlap_invalid && context_eligible_count > 0L &&
                            context_positive >= values$minimum_class_count &&
                            context_negative >= values$minimum_class_count &&
                            context_positive / context_eligible_count >=
                                values$minimum_class_fraction &&
                            context_negative / context_eligible_count >=
                                values$minimum_class_fraction
                        context_data <- context[context_eligible, .(
                            chrom, anchor_start, anchor_score,
                            outcome = get(outcome_column),
                            cofactor_present = as.integer(cofactor_positive),
                            tp73_confirmed = as.integer(get(confirmed_column)),
                            analysis_context_class,
                            log_nearest_tss_genomic_distance,
                            log_nearest_cds_genomic_distance,
                            nearest_tss_direction_class,
                            nearest_cds_direction_class
                        )]
                        context_data[, genomic_block := paste(
                            chrom, floor(anchor_start / values$block_size), sep = ":"
                        )]
                        context_effect <- if (!negative_reference_observable) {
                            status_result(
                                "negative_reference_below_source_floor",
                                paste("cofactor source floor", motif_source_floor,
                                      "is above the requested negative reference")
                            )
                        } else if (overlap_invalid) {
                            status_result(
                                "overlapping_score_classes",
                                "positive threshold is below the negative reference"
                            )
                        } else if (!context_supported) {
                            status_result(
                                "underpowered_context_anchor_class",
                                "context-specific positive/negative classes lack support"
                            )
                        } else {
                            clustered_lm(
                                context_data,
                                model_formula(
                                    context_data, adjust_context_class = FALSE
                                ),
                                list(c(cofactor_present = 1))
                            )[[1L]]
                        }
                        context_intensity_rows[[context_intensity_index]] <- cbind(
                            data.table(
                                motif_id = threshold$motif_id,
                                factor_name = threshold$factor_name,
                                positive_threshold = threshold$positive_threshold,
                                source_score_floor = motif_source_floor,
                                negative_reference_threshold = negative_reference,
                                negative_reference_observable =
                                    negative_reference_observable,
                                isoform = isoform,
                                series_id = series_id,
                                genomic_context_class = context_class,
                                anchors_positive = context_positive,
                                anchors_negative = context_negative,
                                anchors_intermediate =
                                    sum(context$cofactor_intermediate),
                                mean_change_positive = if (context_positive) {
                                    mean(context[[outcome_column]][
                                        context$cofactor_positive
                                    ])
                                } else NA_real_,
                                mean_change_negative = if (context_negative) {
                                    mean(context[[outcome_column]][
                                        context$cofactor_negative
                                    ])
                                } else NA_real_
                            ),
                            context_effect
                        )
                        context_intensity_index <- context_intensity_index + 1L
                    }

                    # The four-way gene relation is the stable reporting
                    # projection. It keeps promoter/downstream precedence
                    # explicit while the richer context class remains a model
                    # covariate inside the heterogeneous gene-body stratum.
                    for (relation_class in c(
                        "promoter", "downstream", "gene_body", "intergenic"
                    )) {
                        relation <- current[
                            gene_relation_class == relation_class
                        ]
                        relation_eligible <- relation$cofactor_positive |
                            relation$cofactor_negative
                        relation_positive <- sum(relation$cofactor_positive)
                        relation_negative <- sum(relation$cofactor_negative)
                        relation_eligible_count <- sum(relation_eligible)
                        relation_supported <- negative_reference_observable &&
                            !overlap_invalid && relation_eligible_count > 0L &&
                            relation_positive >= values$minimum_class_count &&
                            relation_negative >= values$minimum_class_count &&
                            relation_positive / relation_eligible_count >=
                                values$minimum_class_fraction &&
                            relation_negative / relation_eligible_count >=
                                values$minimum_class_fraction
                        relation_data <- relation[relation_eligible, .(
                            chrom, anchor_start, anchor_score,
                            outcome = get(outcome_column),
                            cofactor_present = as.integer(cofactor_positive),
                            tp73_confirmed = as.integer(get(confirmed_column)),
                            analysis_context_class,
                            log_nearest_tss_genomic_distance,
                            log_nearest_cds_genomic_distance,
                            nearest_tss_direction_class,
                            nearest_cds_direction_class
                        )]
                        relation_data[, genomic_block := paste(
                            chrom, floor(anchor_start / values$block_size), sep = ":"
                        )]
                        relation_effect <- if (!negative_reference_observable) {
                            status_result(
                                "negative_reference_below_source_floor",
                                paste("cofactor source floor", motif_source_floor,
                                      "is above the requested negative reference")
                            )
                        } else if (overlap_invalid) {
                            status_result(
                                "overlapping_score_classes",
                                "positive threshold is below the negative reference"
                            )
                        } else if (!relation_supported) {
                            status_result(
                                "underpowered_gene_relation_anchor_class",
                                paste0(
                                    "gene-relation-specific positive/negative ",
                                    "classes lack support"
                                )
                            )
                        } else {
                            clustered_lm(
                                relation_data,
                                model_formula(relation_data),
                                list(c(cofactor_present = 1))
                            )[[1L]]
                        }
                        gene_relation_intensity_rows[[
                            gene_relation_intensity_index
                        ]] <- cbind(
                            data.table(
                                motif_id = threshold$motif_id,
                                factor_name = threshold$factor_name,
                                positive_threshold = threshold$positive_threshold,
                                source_score_floor = motif_source_floor,
                                negative_reference_threshold = negative_reference,
                                negative_reference_observable =
                                    negative_reference_observable,
                                isoform = isoform,
                                series_id = series_id,
                                gene_relation_class = relation_class,
                                anchors_positive = relation_positive,
                                anchors_negative = relation_negative,
                                anchors_intermediate =
                                    sum(relation$cofactor_intermediate),
                                mean_change_positive = if (relation_positive) {
                                    mean(relation[[outcome_column]][
                                        relation$cofactor_positive
                                    ])
                                } else NA_real_,
                                mean_change_negative = if (relation_negative) {
                                    mean(relation[[outcome_column]][
                                        relation$cofactor_negative
                                    ])
                                } else NA_real_
                            ),
                            relation_effect
                        )
                        gene_relation_intensity_index <-
                            gene_relation_intensity_index + 1L
                    }
                }

                score_data <- current[, .(
                    chrom, anchor_start, anchor_score,
                    outcome = get(outcome_column), context_score,
                    analysis_context_class,
                    log_nearest_tss_genomic_distance,
                    log_nearest_cds_genomic_distance,
                    nearest_tss_direction_class,
                    nearest_cds_direction_class
                )]
                score_data[, cofactor_score_clamped := fifelse(
                    is.na(context_score), negative_reference,
                    pmax(context_score, negative_reference)
                )]
                score_center <- mean(score_data$cofactor_score_clamped)
                score_scale <- sd(score_data$cofactor_score_clamped)
                score_data[, cofactor_score_standardized := if (
                    is.finite(score_scale) && score_scale > 0
                ) {
                    (cofactor_score_clamped - score_center) / score_scale
                } else 0]
                score_data[, genomic_block := paste(
                    chrom, floor(anchor_start / values$block_size), sep = ":"
                )]
                score_effect <- if (!negative_reference_observable) {
                    status_result(
                        "negative_reference_below_source_floor",
                        paste("cofactor source floor", motif_source_floor,
                              "is above the requested score-clamp reference")
                    )
                } else if (!is.finite(score_scale) || score_scale <= 0 ||
                           nrow(score_data) < 2L * values$minimum_class_count) {
                    status_result(
                        "underpowered_score_gradient",
                        "clamped cofactor score lacks variation or observations"
                    )
                } else {
                    clustered_lm(
                        score_data,
                        model_formula(
                            score_data,
                            exposure = "cofactor_score_standardized"
                        ),
                        list(c(cofactor_score_standardized = 1))
                    )[[1L]]
                }
                score_effect_rows[[score_effect_index]] <- cbind(data.table(
                    motif_id = threshold$motif_id,
                    factor_name = threshold$factor_name,
                    positive_threshold = threshold$positive_threshold,
                    source_score_floor = motif_source_floor,
                    score_clamp_reference = negative_reference,
                    score_clamp_reference_observable =
                        negative_reference_observable,
                    isoform = isoform,
                    series_id = series_id,
                    anchors = nrow(score_data),
                    anchors_with_observed_context_score =
                        sum(!is.na(score_data$context_score)),
                    cofactor_score_center = score_center,
                    cofactor_score_scale = score_scale,
                    estimate_unit = "one_SD_increase_in_clamped_cofactor_score"
                ), score_effect)
                score_effect_index <- score_effect_index + 1L

                cell_counts <- model_data[, .N, by = .(cofactor_present, tp73_confirmed)]
                four_cells <- nrow(cell_counts) == 4L &&
                    min(cell_counts$N) >= values$minimum_interaction_cell_count
                if (!class_supported) {
                    interaction <- replicate(3L, total, simplify = FALSE)
                } else if (!four_cells) {
                    interaction <- replicate(3L, status_result(
                        "underpowered_interaction",
                        "one of cofactor+/- by TP73-confirmed+/- cells is too small",
                        nrow(model_data), uniqueN(model_data$genomic_block)
                    ), simplify = FALSE)
                } else {
                    interaction <- clustered_lm(
                        model_data,
                        model_formula(model_data, interaction = TRUE),
                        list(
                            c(cofactor_present = 1),
                            c(cofactor_present = 1,
                              "cofactor_present:tp73_confirmed" = 1),
                            c("cofactor_present:tp73_confirmed" = 1)
                        )
                    )
                }
                contrast_names <- c(
                    "cofactor_effect_unconfirmed", "cofactor_effect_confirmed",
                    "cofactor_by_tp73_confirmation_interaction"
                )
                for (contrast_index in seq_along(contrast_names)) {
                    interaction_rows[[interaction_index]] <- cbind(data.table(
                        motif_id = threshold$motif_id,
                        factor_name = threshold$factor_name,
                        positive_threshold = threshold$positive_threshold,
                        source_score_floor = motif_source_floor,
                        negative_reference_threshold = negative_reference,
                        negative_reference_observable =
                            negative_reference_observable,
                        isoform = isoform, series_id = series_id,
                        contrast = contrast_names[[contrast_index]],
                        cell_cofactor0_tp730 = cell_counts[
                            cofactor_present == 0 & tp73_confirmed == 0, sum(N)
                        ],
                        cell_cofactor0_tp731 = cell_counts[
                            cofactor_present == 0 & tp73_confirmed == 1, sum(N)
                        ],
                        cell_cofactor1_tp730 = cell_counts[
                            cofactor_present == 1 & tp73_confirmed == 0, sum(N)
                        ],
                        cell_cofactor1_tp731 = cell_counts[
                            cofactor_present == 1 & tp73_confirmed == 1, sum(N)
                        ]
                    ), interaction[[contrast_index]])
                    interaction_index <- interaction_index + 1L
                }

                if (!negative_reference_observable || overlap_invalid) {
                    state_summary <- data.table(
                        cofactor_present = NA_integer_, binding_state = NA_character_,
                        anchors = NA_integer_, mean_change = NA_real_,
                        median_change = NA_real_
                    )
                } else {
                    state_summary <- model_data[, .(
                        anchors = .N,
                        mean_change = mean(outcome),
                        median_change = median(outcome)
                    ), by = .(cofactor_present, binding_state)]
                }
                state_summary[, `:=`(
                    motif_id = threshold$motif_id,
                    factor_name = threshold$factor_name,
                    positive_threshold = threshold$positive_threshold,
                    source_score_floor = motif_source_floor,
                    negative_reference_threshold = negative_reference,
                    negative_reference_observable = negative_reference_observable,
                    isoform = isoform, series_id = series_id,
                    comparison_status = comparison_status,
                    comparison_note = comparison_note
                )]
                state_rows[[state_index]] <- state_summary
                state_index <- state_index + 1L

                if (!negative_reference_observable || overlap_invalid) {
                    occurrence_summary <- data.table(
                        cofactor_present = NA_integer_, tp73_confirmed = NA_integer_,
                        occurrence = NA_character_, N = NA_integer_,
                        fraction_within_class = NA_real_
                    )
                } else {
                    occurrence_data <- all_current[
                        cofactor_positive | cofactor_negative,
                        .(
                            cofactor_present = as.integer(cofactor_positive),
                            tp73_confirmed = as.integer(get(confirmed_column)),
                            occurrence = get(occurrence_column)
                        )
                    ]
                    occurrence_summary <- occurrence_data[, .N,
                        by = .(cofactor_present, tp73_confirmed, occurrence)]
                    occurrence_summary[, fraction_within_class := N / sum(N),
                                       by = .(cofactor_present, tp73_confirmed)]
                }
                occurrence_summary[, `:=`(
                    motif_id = threshold$motif_id,
                    factor_name = threshold$factor_name,
                    positive_threshold = threshold$positive_threshold,
                    source_score_floor = motif_source_floor,
                    negative_reference_threshold = negative_reference,
                    negative_reference_observable = negative_reference_observable,
                    isoform = isoform, series_id = series_id,
                    comparison_status = comparison_status,
                    comparison_note = comparison_note
                )]
                occurrence_rows[[occurrence_index]] <- occurrence_summary
                occurrence_index <- occurrence_index + 1L
            }
        }
    }
    message("I: analyzed motif ", threshold$motif_id, " (", motif_index, "/",
            nrow(thresholds), ")")
}

intensity <- rbindlist(intensity_rows, fill = TRUE)
interaction <- rbindlist(interaction_rows, fill = TRUE)
binding_state <- rbindlist(state_rows, fill = TRUE)
occurrence <- rbindlist(occurrence_rows, fill = TRUE)
context_intensity <- if (annotation_available) {
    rbindlist(context_intensity_rows, fill = TRUE)
} else NULL
gene_relation_intensity <- if (annotation_available) {
    rbindlist(gene_relation_intensity_rows, fill = TRUE)
} else NULL
score_effect <- rbindlist(score_effect_rows, fill = TRUE)

intensity[, q_value_bh := p.adjust(p_value, method = "BH"),
          by = .(series_id, isoform, negative_reference_threshold)]
interaction[, q_value_bh := p.adjust(p_value, method = "BH"),
            by = .(series_id, isoform, negative_reference_threshold, contrast)]
if (annotation_available) {
    context_intensity[, q_value_bh := p.adjust(p_value, method = "BH"),
        by = .(series_id, isoform, negative_reference_threshold,
               genomic_context_class)]
    gene_relation_intensity[, q_value_bh := p.adjust(p_value, method = "BH"),
        by = .(series_id, isoform, negative_reference_threshold,
               gene_relation_class)]
}
score_effect[, q_value_bh := p.adjust(p_value, method = "BH"),
    by = .(series_id, isoform, score_clamp_reference)]
series_summary <- intensity[, .(
    series_estimable = sum(evaluation_status == "ok"),
    series_total = .N,
    mean_adjusted_change_difference = mean(estimate[evaluation_status == "ok"],
                                           na.rm = TRUE),
    minimum_adjusted_change_difference = suppressWarnings(
        min(estimate[evaluation_status == "ok"], na.rm = TRUE)
    ),
    maximum_adjusted_change_difference = suppressWarnings(
        max(estimate[evaluation_status == "ok"], na.rm = TRUE)
    ),
    series_positive = sum(estimate > 0 & evaluation_status == "ok"),
    series_negative = sum(estimate < 0 & evaluation_status == "ok"),
    direction_consistency = if (sum(evaluation_status == "ok") == .N &&
        all(estimate[evaluation_status == "ok"] > 0)) {
        "all_series_positive"
    } else if (sum(evaluation_status == "ok") == .N &&
               all(estimate[evaluation_status == "ok"] < 0)) {
        "all_series_negative"
    } else if (sum(evaluation_status == "ok") == .N) {
        "mixed_direction"
    } else {
        "not_all_series_estimable"
    }
), by = .(
    motif_id, factor_name, positive_threshold, source_score_floor,
    negative_reference_threshold, negative_reference_observable, isoform
)]
for (column in c("mean_adjusted_change_difference", "minimum_adjusted_change_difference",
                 "maximum_adjusted_change_difference")) {
    set(series_summary, which(!is.finite(series_summary[[column]])), column, NA_real_)
}

setorder(intensity, motif_id, negative_reference_threshold, isoform, series_id)
setorder(interaction, motif_id, negative_reference_threshold, isoform, series_id,
         contrast)
setorder(series_summary, motif_id, negative_reference_threshold, isoform)
setorder(binding_state, motif_id, negative_reference_threshold, isoform, series_id,
         cofactor_present, binding_state)
setorder(occurrence, motif_id, negative_reference_threshold, isoform, series_id,
         cofactor_present, tp73_confirmed, occurrence)
if (annotation_available) {
    setorder(context_intensity, motif_id, negative_reference_threshold, isoform,
             series_id, genomic_context_class)
    setorder(gene_relation_intensity, motif_id, negative_reference_threshold,
             isoform, series_id, gene_relation_class)
}
setorder(score_effect, motif_id, score_clamp_reference, isoform, series_id)

fwrite(intensity, paste0(values$output_prefix, "_intensity_effect.tsv"),
       sep = "\t", na = "NA", quote = FALSE)
fwrite(interaction, paste0(values$output_prefix, "_tp73_interaction.tsv"),
       sep = "\t", na = "NA", quote = FALSE)
fwrite(series_summary, paste0(values$output_prefix, "_series_summary.tsv"),
       sep = "\t", na = "NA", quote = FALSE)
fwrite(binding_state, paste0(values$output_prefix, "_binding_state_summary.tsv"),
       sep = "\t", na = "NA", quote = FALSE)
fwrite(occurrence, paste0(values$output_prefix, "_occurrence_summary.tsv"),
       sep = "\t", na = "NA", quote = FALSE)
if (annotation_available) {
    fwrite(
        context_intensity,
        paste0(values$output_prefix, "_context_stratified_intensity_effect.tsv"),
        sep = "\t", na = "NA", quote = FALSE
    )
    fwrite(
        gene_relation_intensity,
        paste0(
            values$output_prefix,
            "_gene_relation_stratified_intensity_effect.tsv"
        ),
        sep = "\t", na = "NA", quote = FALSE
    )
}
fwrite(score_effect, paste0(values$output_prefix, "_score_gradient.tsv"),
       sep = "\t", na = "NA", quote = FALSE)

run_config <- data.table(
    schema_version = 3L,
    analysis = "gfp_referenced_h3k4me3_cofactor_change",
    analysis_role = values$analysis_role,
    chromosomes = paste(sort(unique(wide$chrom)), collapse = ","),
    input_mode = input_mode,
    signal = if (has_signal) paste(values$signal, collapse = ";") else NA_character_,
    change = if (has_change) paste(values$change, collapse = ";") else NA_character_,
    tp73_evidence = values$tp73_evidence,
    cofactor_maxima = paste(values$cofactor_maxima, collapse = ";"),
    annotation = if (annotation_available) {
        paste(values$annotation, collapse = ";")
    } else NA_character_,
    annotation_schema = if (annotation_available) {
        "schema9_tp73_context_anchor_collapsed_to_physical_span"
    } else "not_supplied",
    thresholds = values$thresholds,
    window_name = values$window,
    pseudocount = values$pseudocount,
    signal_formula = "log2((h3k4me3_area + alpha)/(input_area + alpha))",
    primary_outcome = "condition_minus_GFP_normalized_H3K4me3",
    primary_estimand = "cofactor_positive_minus_negative_difference_in_change",
    secondary_score_estimand =
        "change_per_SD_of_cofactor_score_clamped_at_observable_negative_reference",
    primary_adjustment = if (annotation_available) {
        paste0(
            "TP73_score_spline + chromosome + genomic_context + ",
            "nearest_TSS/CDS_unsigned_distance + direction_class"
        )
    } else "TP73_score_spline + chromosome_when_multichromosome",
    context_stratified_output = annotation_available,
    gene_relation_stratified_output = annotation_available,
    gene_relation_precedence = if (annotation_available) {
        "promoter_then_downstream_then_gene_body_then_intergenic"
    } else NA_character_,
    strict_intergenic_definition = if (annotation_available) {
        "no_transcript_promoter_or_downstream_region_overlap"
    } else NA_character_,
    tp73_confirmation = "strict TP73 immersion AND no strict negative-control immersion",
    tp73_interaction_interpretation = "secondary_descriptive_post_treatment",
    all_zero_anchor_policy =
        "excluded_from_intensity_and_score_gradient_retained_for_occurrence",
    series = paste(analysis_series, collapse = ","),
    block_size_bp = values$block_size,
    spline_df = values$spline_df,
    negative_references = paste(negative_references, collapse = ","),
    minimum_class_fraction = values$minimum_class_fraction,
    minimum_class_count = values$minimum_class_count,
    minimum_interaction_cell_count = values$minimum_interaction_cell_count,
    motifs_in_threshold_input = nrow(thresholds),
    multiple_testing_scope =
        "motifs_in_threshold_input_within_declared_series_isoform_reference_family"
)
fwrite(run_config, paste0(values$output_prefix, "_run_config.tsv"),
       sep = "\t", na = "NA", quote = FALSE)

message("I: wrote GFP-referenced H3K4me3 cofactor results: ", values$output_prefix)
