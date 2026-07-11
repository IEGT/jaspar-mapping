#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function(output = stdout()) {
  cat(
    paste0(
      "Usage: plot_tp73_score_distributions.R [LOG_RISK.tsv] [LOG_ODDS.tsv] ",
      "[OUTPUT.png] [TAIL_SUMMARY.tsv]\n\n",
      "Plot TP73/MA0861.2 score-density and cumulative-tail curves from two ",
      "pssm_scan score-distribution tables. Positional arguments override, ",
      "in order, the default log2-relative-risk input, log-odds input, PNG ",
      "plot, and tail-summary TSV paths.\n\n",
      "Options:\n",
      "  -h, --help  Show this help and exit\n"
    ),
    file = output
  )
}

if (length(args) == 1L && args[[1L]] %in% c("-h", "--help")) {
  usage()
  quit(status = 0L)
}
if (length(args) > 4L) {
  cat("E: Expected no more than four positional arguments.\n", file = stderr())
  usage(stderr())
  quit(status = 2L)
}

log_risk_file <- if (length(args) >= 1) args[[1]] else
  "score_distributions_log2_relative_risk_bins_adaptive_JASPAR2026_chr1_pseudocount_1/TP73_MA0861.2_score_distribution_log2_relative_risk_bins_adaptive_both_1_pseudocount_1.tsv"
log_odds_file <- if (length(args) >= 2) args[[2]] else
  "score_distributions_log_odds_bins_adaptive_JASPAR2026_chr1_pseudocount_1/TP73_MA0861.2_score_distribution_log_odds_bins_adaptive_both_1_pseudocount_1.tsv"
out_png <- if (length(args) >= 3) args[[3]] else
  "images/tp73_chr1_score_distribution_density_cumulative.png"
out_tsv <- if (length(args) >= 4) args[[4]] else
  "score_distributions_log_odds_bins_adaptive_JASPAR2026_chr1_pseudocount_1/TP73_MA0861.2_chr1_tail_summary.tsv"

read_distribution <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing score distribution file: ", path)
  }

  data <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  starts <- suppressWarnings(as.numeric(data$ScoreBinStart))
  ends <- suppressWarnings(as.numeric(data$ScoreBinEnd))
  finite_bins <- is.finite(starts) & is.finite(ends)
  data <- data[finite_bins, , drop = FALSE]
  data$ScoreBinStart <- starts[finite_bins]
  data$ScoreBinEnd <- ends[finite_bins]
  data$BinCount <- as.numeric(data$BinCount)
  data$ValidWindows <- as.numeric(data$ValidWindows)
  data$BinWidthComputed <- data$ScoreBinEnd - data$ScoreBinStart
  data$ScoreBinMidpoint <- (data$ScoreBinStart + data$ScoreBinEnd) / 2
  data$BinDensity <- data$BinCount / data$BinWidthComputed
  data <- data[order(data$ScoreBinStart), , drop = FALSE]
  data$CumulativeCountGE <- rev(cumsum(rev(data$BinCount)))
  data$CumulativeFractionGE <- data$CumulativeCountGE / data$ValidWindows[1]
  data$ScoreModeLabel <- label
  data
}

tail_summary <- function(data, thresholds) {
  rows <- lapply(thresholds, function(threshold) {
    cumulative_count <- sum(data$BinCount[data$ScoreBinEnd > threshold])
    data.frame(
      ScoreMode = data$ScoreMode[1],
      MotifID = data$MotifID[1],
      MotifName = data$MotifName[1],
      Chromosome = data$Chromosome[1],
      Strand = data$Strand[1],
      Pseudocount = data$Pseudocount[1],
      Threshold = threshold,
      CumulativeCountGE = cumulative_count,
      CumulativeFractionGE = cumulative_count / data$ValidWindows[1]
    )
  })
  do.call(rbind, rows)
}

log_risk <- read_distribution(log_risk_file, "log2 relative risk")
log_odds <- read_distribution(log_odds_file, "log odds")

summary_table <- rbind(
  tail_summary(log_risk, c(-10, -5, 0, 5, 10, 15, 20)),
  tail_summary(log_odds, c(0, 10, 20, 25, 30, 35, 40, 45))
)

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
write.table(summary_table, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
png(out_png, width = 1800, height = 1200, res = 180)

old_par <- par(no.readonly = TRUE)
on.exit({
  par(old_par)
  dev.off()
}, add = TRUE)

colors <- c("log2 relative risk" = "#1F77B4", "log odds" = "#D62728")
x_limits <- range(c(log_risk$ScoreBinStart, log_risk$ScoreBinEnd,
                    log_odds$ScoreBinStart, log_odds$ScoreBinEnd))

par(mfrow = c(2, 1), mar = c(4.2, 5.0, 3.0, 1.2), oma = c(0, 0, 2.0, 0))

density_values <- c(log_risk$BinDensity, log_odds$BinDensity)
density_positive <- density_values[density_values > 0]
plot(log_risk$ScoreBinMidpoint, log10(log_risk$BinDensity),
     type = "l", lwd = 2.3, col = colors[["log2 relative risk"]],
     xlim = x_limits, ylim = range(log10(density_positive)),
     xlab = "Score", ylab = "log10(bin density)",
     main = "TP73 / MA0861.2 chr1 score density")
lines(log_odds$ScoreBinMidpoint, log10(log_odds$BinDensity),
      lwd = 2.3, col = colors[["log odds"]])
abline(v = -10, col = colors[["log2 relative risk"]], lty = 2)
abline(v = 0, col = colors[["log2 relative risk"]], lty = 3)
abline(v = 35, col = colors[["log odds"]], lty = 2)
legend("topright",
       legend = c("log2 relative risk", "log odds", "risk bin-width switch", "risk zero", "odds +35"),
       col = c(colors[["log2 relative risk"]], colors[["log odds"]],
               colors[["log2 relative risk"]], colors[["log2 relative risk"]], colors[["log odds"]]),
       lwd = c(2.3, 2.3, 1.2, 1.2, 1.2),
       lty = c(1, 1, 2, 3, 2),
       bty = "n", cex = 0.85)
mtext("Bin counts are normalized by bin width; the raw-count step at -10 is a binning boundary.",
      side = 3, line = 0.2, cex = 0.8)

tail_positive <- c(log_risk$CumulativeFractionGE, log_odds$CumulativeFractionGE)
tail_positive <- tail_positive[tail_positive > 0]
plot(log_risk$ScoreBinStart, log_risk$CumulativeFractionGE,
     type = "l", lwd = 2.3, col = colors[["log2 relative risk"]],
     log = "y", xlim = x_limits, ylim = range(tail_positive),
     xlab = "Score threshold", ylab = "Fraction with score >= threshold",
     main = "Cumulative tail fraction")
lines(log_odds$ScoreBinStart, log_odds$CumulativeFractionGE,
      lwd = 2.3, col = colors[["log odds"]])
abline(v = -10, col = colors[["log2 relative risk"]], lty = 2)
abline(v = 0, col = colors[["log2 relative risk"]], lty = 3)
abline(v = 35, col = colors[["log odds"]], lty = 2)
grid(col = "grey85")

mtext("TP73 / MA0861.2 chr1 both strands, pseudocount 1", outer = TRUE, cex = 1.1, font = 2)

message("Wrote plot: ", out_png)
message("Wrote summary table: ", out_tsv)
