#!/usr/bin/env Rscript

usage <- function(status = 0L) {
    stream <- if (status == 0L) stdout() else stderr()
    writeLines(c(
        "Usage: plot_tp73_distance_frequency_enrichment.R --input TABLE.tsv",
        "       --output-adjusted FIGURE.png",
        "       --output-frequency-ratio FIGURE.png --output-table TABLE.tsv",
        "       [--human-source-only] [--label-motifs ID[,ID...]]",
        "",
        "Plot motif frequency near TP73 anchors against CUT&RUN enrichment or",
        "depletion, separately by p73 isoform and exclusive distance band.",
        "The adjusted figure uses the primary matched log2 odds ratio. The",
        "frequency-ratio figure is a descriptive compatibility view analogous",
        "to the historical frequency-versus-log-ratio plot.",
        "",
        "Options:",
        "  --input FILE                    Schema-5 frequency-enrichment TSV",
        "  --output-adjusted FILE          Primary adjusted-odds scatter plot",
        "  --output-frequency-ratio FILE   Descriptive log-ratio scatter plot",
        "  --output-table FILE             Exact plotted rows and derived fields",
        "  --human-source-only             Keep matrices sourced from human data",
        "  --label-motifs IDS              Comma-separated motif IDs to label",
        "  -h, --help                      Show this help"
    ), con = stream)
    quit(status = status)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (any(arguments %in% c("-h", "--help"))) usage()
values <- list(
    human_source_only = FALSE,
    label_motifs = paste(
        c("MA0024.3", "MA0507.3", "MA0079.5", "MA1961.2", "MA0138.3"),
        collapse = ","
    )
)
value_options <- c(
    "--input", "--output-adjusted", "--output-frequency-ratio",
    "--output-table", "--label-motifs"
)
index <- 1L
while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (option == "--human-source-only") {
        values$human_source_only <- TRUE
        index <- index + 1L
        next
    }
    if (!option %in% value_options) {
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
    "input", "output_adjusted", "output_frequency_ratio", "output_table"
)
if (any(vapply(required_arguments, function(key) is.null(values[[key]]),
               logical(1)))) {
    usage(2L)
}
if (!file.exists(values$input)) {
    stop("frequency-enrichment input not found", call. = FALSE)
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

source <- fread(values$input, sep = "\t", na.strings = c("NA", ""))
required_columns <- c(
    "motif_id", "motif_name", "isoform", "distance_band",
    "distance_band_order", "positive_threshold",
    "all_tp73_anchor_vicinity_frequency",
    "anti_supported_positive_anchor_fraction_discordant",
    "control_supported_positive_anchor_fraction_discordant",
    "anti_to_control_positive_anchor_log2_ratio_discordant",
    "adjusted_log2_odds", "q_value_bh_tax_group", "evaluation_status",
    "includes_homo_sapiens"
)
missing_columns <- setdiff(required_columns, names(source))
if (nrow(source) == 0L || length(missing_columns) > 0L) {
    stop(
        paste("frequency-enrichment table is empty or incomplete:",
              paste(missing_columns, collapse = ", ")),
        call. = FALSE
    )
}
if (values$human_source_only) source <- source[includes_homo_sapiens == TRUE]

distance_levels <- source[
    order(distance_band_order), unique(as.character(distance_band))
]
plot_data <- source[
    evaluation_status == "ok" & is.finite(adjusted_log2_odds) &
        is.finite(anti_supported_positive_anchor_fraction_discordant)
]
if (nrow(plot_data) == 0L) {
    stop("no estimable frequency/enrichment rows remain", call. = FALSE)
}
plot_data[, `:=`(
    isoform = factor(isoform, levels = c("TA", "DN")),
    distance_band = factor(distance_band, levels = distance_levels),
    anti_supported_frequency_percent =
        100 * anti_supported_positive_anchor_fraction_discordant,
    control_supported_frequency_percent =
        100 * control_supported_positive_anchor_fraction_discordant,
    all_anchor_vicinity_frequency_percent =
        100 * all_tp73_anchor_vicinity_frequency,
    association_class = fcase(
        q_value_bh_tax_group > 0.05, "q > 0.05",
        adjusted_log2_odds > 0, "anti-p73 enriched",
        adjusted_log2_odds < 0, "anti-p73 depleted",
        default = "neutral"
    )
)]
plot_data[, association_class := factor(
    association_class,
    levels = c("anti-p73 depleted", "q > 0.05", "anti-p73 enriched", "neutral")
)]
labels <- trimws(strsplit(values$label_motifs, ",", fixed = TRUE)[[1L]])
labels <- labels[nzchar(labels)]
label_data <- plot_data[motif_id %in% labels]
label_data[, label := paste0(motif_name, "\n", motif_id)]

palette <- c(
    "anti-p73 depleted" = "#b2182b",
    "q > 0.05" = "#8a8a8a",
    "anti-p73 enriched" = "#2166ac",
    "neutral" = "#4d4d4d"
)
common_layers <- list(
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55"),
    geom_point(alpha = 0.48, size = 1.15),
    geom_text(
        data = label_data, aes(label = label), colour = "black", size = 2.1,
        check_overlap = TRUE, vjust = -0.65, show.legend = FALSE
    ),
    facet_grid(rows = vars(isoform), cols = vars(distance_band), scales = "free_y"),
    scale_colour_manual(values = palette, drop = FALSE),
    theme_bw(base_size = 9.5),
    theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = "grey70"),
        legend.position = "bottom"
    )
)

adjusted_figure <- ggplot(
    plot_data,
    aes(adjusted_log2_odds, anti_supported_frequency_percent,
        colour = association_class)
) +
    common_layers +
    labs(
        x = expression(log[2] * " adjusted matched TP73 CUT&RUN odds ratio"),
        y = "Motif-positive anti-p73-supported anchors (%)",
        colour = NULL,
        title = "Cofactor frequency and TP73 CUT&RUN association",
        subtitle = paste0(
            "Frequency is calculated within matched-discordant anti-p73-supported ",
            "anchors; one presence per anchor and distance band"
        )
    )

compatibility_data <- plot_data[
    is.finite(anti_to_control_positive_anchor_log2_ratio_discordant)
]
if (nrow(compatibility_data) == 0L) {
    stop("no finite descriptive frequency ratios remain", call. = FALSE)
}
compatibility_labels <- label_data[
    is.finite(anti_to_control_positive_anchor_log2_ratio_discordant)
]
compatibility_layers <- common_layers
compatibility_layers[[3L]] <- geom_text(
    data = compatibility_labels, aes(label = label), colour = "black",
    size = 2.1, check_overlap = TRUE, vjust = -0.65, show.legend = FALSE
)
compatibility_figure <- ggplot(
    compatibility_data,
    aes(anti_to_control_positive_anchor_log2_ratio_discordant,
        anti_supported_frequency_percent, colour = association_class)
) +
    compatibility_layers +
    labs(
        x = expression(log[2] * " motif-frequency ratio (anti-p73 / control)"),
        y = "Motif-positive anti-p73-supported anchors (%)",
        colour = NULL,
        title = "Historical-style frequency versus enrichment/depletion",
        subtitle = paste0(
            "Descriptive matched-discordant compatibility view; primary inference ",
            "uses the adjusted odds-ratio figure"
        )
    )

for (path in c(
    values$output_adjusted, values$output_frequency_ratio, values$output_table
)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
ggsave(values$output_adjusted, adjusted_figure, width = 16, height = 7.5, dpi = 180)
ggsave(
    values$output_frequency_ratio, compatibility_figure,
    width = 16, height = 7.5, dpi = 180
)
setorder(plot_data, distance_band_order, isoform, motif_id)
fwrite(plot_data, values$output_table, sep = "\t", na = "NA")
message("I: wrote adjusted frequency/enrichment figure: ", values$output_adjusted)
message("I: wrote descriptive frequency-ratio figure: ",
        values$output_frequency_ratio)
message("I: wrote plotted frequency table: ", values$output_table)
