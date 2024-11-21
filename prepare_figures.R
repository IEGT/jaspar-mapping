# Load necessary library
library(data.table)

# Define the file path
file_path <- "/home/moeller/GitHub/jaspar-mapping/output_Chr1/TP73_MA0861.1_positive_1.combined.bed"

# Read the large tab-separated data file
data <- fread(file_path, sep = "\t")

# Display the first few rows of the data
head(data)

# Load necessary library for plotting
library(ggplot2)

dot.size <- 0.1
# Create a dot plot of "Score" against "tp73_skmel29_2_DN"
ggplot(data, aes(y = log2(1+tp73_skmel29_2_DN), x = Score)) +
    geom_point(size = dot.size) +
    theme_minimal() +
    labs(title = "Dot Plot of Score vs tp73_skmel29_2_DN",
             y = "tp73_skmel29_2_DN",
             x = "Score")
# Save the last plot as a PNG file
ggsave("scatter_cutnrun_DN_vs_PSSM_score.png", width = 8, height = 6)
ggsave("scatter_cutnrun_DN_vs_PSSM_score.svg", width = 8, height = 6)

# Create a dot plot of "Score" against "tp73_skmel29_2_TA"
ggplot(data, aes(y = log2(1+tp73_skmel29_2_TA), x = Score)) +
    geom_point(size = dot.size) +
    theme_minimal() +
    labs(title = "Dot Plot TAalpha Cut&Run coverage against TP73 predicted binding sites",
             y = "tp73_skmel29_2_TA",
             x = "Score")
# Save the last plot as a PNG file
ggsave("scatter_cutnrun_TA_vs_PSSM_score.png", width = 8, height = 6)
ggsave("scatter_cutnrun_DN_vs_PSSM_score.svg", width = 8, height = 6)

# Create a dot plot of "Score" against "tp73_skmel29_2_TA"-"tp73_skmel29_2_DN"
ggplot(data, aes(y = log2(1+tp73_skmel29_2_TA)-log2(1+tp73_skmel29_2_DN), x = Score)) +
    geom_point(size = dot.size) +
    theme_minimal() +
    labs(title = "Dot Plot TAalpha-DNbeta Cut&Run coverage against TP73 predicted binding sites",
             y = "tp73_skmel29_2_TA-tp73_skmel29_2_DN",
             x = "Score")

# Save the last plot as a PNG file
ggsave("scatter_cutnrun_TA-DN_vs_PSSM_score.png", width = 8, height = 6)
ggsave("scatter_cutnrun_TA-DN_vs_PSSM_score.svg", width = 8, height = 6)

# Create a scatter plot of "tp73_skmel29_2_TA" against "tp73_skmel29_2_DN" with "Score" as color
ggplot(data, aes(x = log2(1+tp73_skmel29_2_DN), y = log2(1+tp73_skmel29_2_TA), color = Score)) +
    geom_point(size = dot.size) +
    theme_minimal() +
    labs(title = "Scatter Plot of tp73_skmel29_2_TA vs tp73_skmel29_2_DN with Score as Color",
                x = "tp73_skmel29_2_DN",
                y = "tp73_skmel29_2_TA",
                color = "Score") +
    scale_color_gradientn(colors = rainbow(7))

score.q <- quantile(data$Score, probs = c(0, 0.25, 0.5, 0.75, 1))

# Define score ranges
score_ranges <- list(
    "0-25%"   = score.q[1:2],
    "25-50%"  = score.q[2:3],
    "50-75%"  = score.q[3:4],
    "75-100%" = score.q[4:5]
)

# Create a list to store plots
plots <- list()

# Generate plots for each score range
for (range_name in names(score_ranges)) {
    range <- score_ranges[[range_name]]
    filtered_data <- data[Score >= range[1] & Score < range[2]]
    filtered_data.cor <- cor.test(filtered_data$tp73_skmel29_2_TA, filtered_data$tp73_skmel29_2_DN,method="spearman",exact=F)
    p <- ggplot(filtered_data, aes(x = log2(1+tp73_skmel29_2_DN), y = log2(1+tp73_skmel29_2_TA), color = Score)) +
        geom_point(size = dot.size) +
        theme_minimal() +
        labs(title = paste(range_name, cor=signif(filtered_data.cor$estimate,3),
                        " p=", signif(filtered_data.cor$p.value,3),
                        " %0(TA)=", signif(sum(filtered_data$tp73_skmel29_2_TA==0)/nrow(filtered_data)*100,3),
                        " %0(DN)=", signif(sum(filtered_data$tp73_skmel29_2_DN==0)/nrow(filtered_data)*100,3),
                        " %0(GFP)=", signif(sum(filtered_data$tp73_skmel29_2_GFP==0)/nrow(filtered_data)*100,3),
                        sep=""
                        ),
                x = "tp73_skmel29_2_DN",
                y = "tp73_skmel29_2_TA",
                color = "Score") +
        scale_color_gradientn(colors = rainbow(7))
    
    plots[[range_name]] <- p
}

require(gridExtra)
png("scatter_cutnrun_TA_vs_DN_score_array.png", width = 1000, height = 1000)
# Arrange and display the plots in a grid
do.call(grid.arrange, c(plots, ncol = 2))
dev.off()

svg("scatter_cutnrun_TA_vs_DN_score_array.svg")
# Arrange and display the plots in a grid
do.call(grid.arrange, c(plots, ncol = 2))
dev.off()


# Save the density plot as a PNG file

#data[, bin := cut(Score, breaks = seq(min(Score), max(Score), by = 1), include.lowest = TRUE)]
# Print the percentage of 0 values in "tp73_skmel29_2_TA" and "tp73_skmel29_2_DN" for each bin
print(percentage_zeros_combined[, .(bin, percentage_zeros_TA, percentage_zeros_DN)])
# Prepare data for plotting
percentage_zeros_combined_long <- melt(percentage_zeros_combined[-1,], id.vars = "bin", variable.name = "Type", value.name = "Percentage")

# Create a line plot of the percentages against the score
ggplot(percentage_zeros_combined_long,
    #aes(x = as.factor(as.character(bin)),
    aes(as.numeric(gsub(pattern=",.*",replacement="",x=gsub(pattern="^[([]",replacement="",x=bin)))
    y = Percentage, color = Type, shape = Type)) +
    geom_line() +
    geom_point(size = dot.size) +
    theme_minimal() +
    labs(title = "Percentage of 0 values in tp73_skmel29_2_TA and tp73_skmel29_2_DN against Score",
         x = "Score",
         y = "Percentage of 0 values") +
    scale_color_manual(values = c("percentage_zeros_TA" = "blue", "percentage_zeros_DN" = "red")) +
    scale_shape_manual(values = c("percentage_zeros_TA" = 16, "percentage_zeros_DN" = 17))

# Save the plot as a PNG file
ggsave("percentage_zeros_line_plot.png", width = 8, height = 6)
ggsave("percentage_zeros_line_plot.svg", width = 8, height = 6)
# Calculate the percentage of 0 values in "tp73_skmel29_2_TA" for each bin
