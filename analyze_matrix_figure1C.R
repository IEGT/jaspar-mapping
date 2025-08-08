#
# Manuscript Figure 1C
#
svg("fractions_by_score.svg", width=800, height=400)
# Plot the fractions    
ggplot(fractions, aes(x = ScoreBin+bin.size/2)) +
    geom_bar(aes(y = TFBS_Count / max(TFBS_Count), fill = "p73 BS Count"), stat = "identity", alpha = 0.5) +
    geom_line(aes(y = Fraction_DN, color = "DN"), linewidth = 1) +
    geom_line(aes(y = Fraction_TA, color = "TA"), linewidth = 1) +
    scale_y_continuous(
        limits = c(0, 1),  # Ensure the y-axis covers 0 to 1 so 0.15 is visible
        name = "Fraction confirmed by CUT&RUN",
        sec.axis = sec_axis(~ . * max(fractions$TFBS_Count), name = "TFBS Count")
    ) +
    labs(
        title = "Fraction of Entries with Value > 0 by Score and TFBS Count",
        x = "Score",
        color = "Legend",
        fill = "Legend"
    ) +
    theme_minimal()
dev.off()
