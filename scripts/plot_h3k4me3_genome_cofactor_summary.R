#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: plot_h3k4me3_genome_cofactor_summary.R --joint JOINT.tsv",
        "       --context CONTEXT.tsv --output-effect EFFECT.png",
        "       --output-context CONTEXT.png --context-table CONTEXT.tsv",
        "       [--label-motifs ID[,ID...]]",
        "",
        "Plot the genome-wide relationship between TP73 CUT&RUN enrichment and",
        "GFP-referenced H3K4me3 cofactor effects. The context plot and table",
        "summarize that relationship separately by genomic annotation class.",
        "Only motifs with an estimable strict score < -1 reference enter these",
        "figures; scan-floor-censored motifs remain in motif_coverage instead.",
        "",
        "Options:",
        "  --joint FILE          joint_primary_motif TSV",
        "  --context FILE        context_primary_effect TSV",
        "  --output-effect FILE  Four-panel enrichment/effect PNG",
        "  --output-context FILE Genomic-context correlation PNG",
        "  --context-table FILE  Derived context-correlation TSV",
        "  --label-motifs IDS    Comma-separated motif IDs to label",
        "  -h, --help            Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) usage()
values <- list(label_motifs = paste(
    c("MA0024.3", "MA0507.3", "MA0138.3", "MA0660.1", "MA0773.1",
      "MA2494.1", "MA1575.2", "MA0730.1", "MA1538.1"),
    collapse = ","
))
index <- 1L
valid <- c(
    "--joint", "--context", "--output-effect", "--output-context",
    "--context-table", "--label-motifs"
)
while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (!option %in% valid) {
        writeLines(paste("E: unknown argument:", option), con = stderr())
        usage(2L)
    }
    index <- index + 1L
    if (index > length(arguments)) usage(2L)
    key <- gsub("-", "_", substring(option, 3L), fixed = TRUE)
    values[[key]] <- arguments[[index]]
    index <- index + 1L
}
required_arguments <- c(
    "joint", "context", "output_effect", "output_context", "context_table"
)
if (any(vapply(required_arguments, function(key) is.null(values[[key]]),
               logical(1)))) {
    usage(2L)
}
if (!file.exists(values$joint) || !file.exists(values$context)) {
    stop("joint or context input not found", call. = FALSE)
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

joint <- fread(values$joint, sep = "\t", na.strings = c("NA", ""))
required_joint <- c(
    "motif_id", "factor_name", "tp73_adjusted_odds_ratio", "tp73_q",
    "tp73_association_direction", "ta_saos", "ta_skmel", "dn_saos",
    "dn_skmel"
)
if (nrow(joint) == 0L || length(setdiff(required_joint, names(joint))) > 0L) {
    stop("joint table is empty or incomplete", call. = FALSE)
}
if (joint[, any(!is.finite(tp73_adjusted_odds_ratio) |
                tp73_adjusted_odds_ratio <= 0)]) {
    stop("joint table contains an invalid TP73 odds ratio", call. = FALSE)
}

effects <- melt(
    joint,
    id.vars = c(
        "motif_id", "factor_name", "tp73_adjusted_odds_ratio", "tp73_q",
        "tp73_association_direction"
    ),
    measure.vars = c("ta_saos", "ta_skmel", "dn_saos", "dn_skmel"),
    variable.name = "panel", value.name = "h3_effect"
)
effects[, `:=`(
    tp73_log2_odds_ratio = log2(tp73_adjusted_odds_ratio),
    isoform = fifelse(grepl("^ta_", panel), "TA", "DN"),
    series = fifelse(grepl("_saos$", panel), "SaOS-2", "SK-Mel-29 series 2"),
    enrichment_class = fcase(
        tp73_q > 0.05, "TP73 association q > 0.05",
        tp73_association_direction == "anti_p73_enriched", "TP73 enriched",
        tp73_association_direction == "anti_p73_depleted", "TP73 depleted",
        default = "Other"
    )
)]
effects[, enrichment_class := factor(
    enrichment_class,
    levels = c("TP73 depleted", "TP73 association q > 0.05", "TP73 enriched")
)]
effects[, isoform := factor(isoform, levels = c("TA", "DN"))]
effects[, series := factor(series, levels = c("SaOS-2", "SK-Mel-29 series 2"))]

labels <- trimws(strsplit(values$label_motifs, ",", fixed = TRUE)[[1L]])
labels <- labels[nzchar(labels)]
label_data <- effects[motif_id %in% labels & isoform == "TA"]
label_data[, label := paste0(factor_name, "\n", motif_id)]

palette <- c(
    "TP73 depleted" = "#b2182b",
    "TP73 association q > 0.05" = "#8a8a8a",
    "TP73 enriched" = "#2166ac"
)
effect_figure <- ggplot(
    effects,
    aes(tp73_log2_odds_ratio, h3_effect, colour = enrichment_class)
) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey65") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey65") +
    geom_point(alpha = 0.55, size = 1.25) +
    geom_smooth(
        data = effects[tp73_q <= 0.05], method = "lm", formula = y ~ x,
        se = FALSE, linewidth = 0.65, colour = "black"
    ) +
    geom_text(
        data = label_data, aes(label = label), colour = "black", size = 2.5,
        check_overlap = TRUE, vjust = -0.65, show.legend = FALSE
    ) +
    facet_grid(rows = vars(series), cols = vars(isoform), scales = "free_y") +
    scale_colour_manual(values = palette, drop = FALSE) +
    labs(
        x = expression(log[2] * " adjusted TP73 CUT&RUN odds ratio"),
        y = expression(Delta * " H3K4me3 cofactor effect vs GFP"),
        colour = NULL,
        title = "TP73 occupancy association and H3K4me3 change",
        subtitle = paste0(
            "Strict cofactor-negative reference: score < -1; line uses motifs ",
            "with TP73 association q <= 0.05"
        )
    ) +
    theme_bw(base_size = 10.5) +
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = "grey70"),
        legend.position = "bottom"
    )

context <- fread(values$context, sep = "\t", na.strings = c("NA", ""))
required_context <- c(
    "motif_id", "tp73_adjusted_odds_ratio", "tp73_q", "series_id",
    "isoform", "genomic_context_class", "evaluation_status", "estimate",
    "q_value_bh_all_motifs"
)
if (nrow(context) == 0L ||
    length(setdiff(required_context, names(context))) > 0L) {
    stop("context table is empty or incomplete", call. = FALSE)
}
context_correlations <- context[
    evaluation_status == "ok" & tp73_q <= 0.05 &
        is.finite(tp73_adjusted_odds_ratio) & tp73_adjusted_odds_ratio > 0 &
        is.finite(estimate),
    .(
        motifs = uniqueN(motif_id),
        pearson_r = cor(log2(tp73_adjusted_odds_ratio), estimate),
        h3_q05 = sum(q_value_bh_all_motifs <= 0.05, na.rm = TRUE)
    ),
    by = .(series_id, isoform, genomic_context_class)
]
context_correlations[, `:=`(
    series = fifelse(series_id == "saos2", "SaOS-2", "SK-Mel-29 series 2"),
    context_label = factor(
        genomic_context_class,
        levels = c("strict_intergenic", "promoter", "intron", "exonic_non_cds", "cds"),
        labels = c("Strict intergenic", "Promoter", "Intron", "Exonic, non-CDS", "CDS")
    )
)]
setorder(context_correlations, series_id, isoform, genomic_context_class)

context_figure <- ggplot(
    context_correlations,
    aes(isoform, context_label, fill = pearson_r)
) +
    geom_tile(colour = "white", linewidth = 0.7) +
    geom_text(aes(label = sprintf("%.2f", pearson_r)), size = 3.2) +
    facet_wrap(vars(series), ncol = 2) +
    scale_fill_gradient2(
        low = "#b2182b", mid = "white", high = "#2166ac",
        midpoint = 0, limits = c(-1, 1)
    ) +
    labs(
        x = "p73 isoform", y = NULL, fill = "Pearson r",
        title = "Association persists within genomic annotation classes",
        subtitle = paste0(
            "Correlation of TP73 log2 odds ratio with H3K4me3 cofactor effect; ",
            "TP73 association q <= 0.05"
        )
    ) +
    theme_bw(base_size = 10.5) +
    theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = "grey70")
    )

for (path in c(values$output_effect, values$output_context, values$context_table)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
ggsave(values$output_effect, effect_figure, width = 10, height = 7, dpi = 180)
ggsave(values$output_context, context_figure, width = 8.5, height = 4.8, dpi = 180)
fwrite(context_correlations, values$context_table, sep = "\t")
message("I: wrote enrichment/effect figure: ", values$output_effect)
message("I: wrote genomic-context figure: ", values$output_context)
