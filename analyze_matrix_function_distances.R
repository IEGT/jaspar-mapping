#!/usr/bin/R


require(ggplot2)

plot.binding.affinity <- function(m) {
    melted_data <- melt(selected_columns_cutnrun)

    # Create the boxplot
    ggplot(melted_data, aes(x = variable, y = value)) +
    geom_boxplot() +
    labs(title = "CUT&RUN signal distribution for TP73 binding sites",
        x = "Columns",
        y = "Values") +
    theme_minimal()

    # Save the plot to a file
    ggsave("boxplot_selected_columns_cutnrun.png")
}



#####
# Creates distance plots for the selected columns

require(ggplot2)

plotShiftOfBinding <- function(m,
                               chromosome=NA,
                               tf="TP73_MA0861.1",
                               tfbs.selection=NULL,
                               subset.name="unknown",
                               binwidth=2,
                               output.dir="images"

) {
    # Extract relevant columns from the data frame
    selected_column_tf_Score_name <- paste(tf,"Score",sep="_")
    selected_column_tf_Shift_name <- paste(tf,"Shift",sep="_")
    selected_column_tf_SameStrand_name <- paste(tf,"StrandEqual",sep="_")
    #selected_column_tf_Score <- m[tfbs.selection, ..selected_column_tf_Score_name][[1]]
    #selected_column_tf_Shift <- m[tfbs.selection, ..selected_column_tf_Shift_name][[1]]
    selected_column_tf_SameStrand <- m[tfbs.selection, ..selected_column_tf_SameStrand_name][[1]]

    plot_data <- data.frame(
        score=m[tfbs.selection, ..selected_column_tf_Score_name][[1]],
        shift=m[tfbs.selection, ..selected_column_tf_Shift_name][[1]]
    )

    if (!is.numeric(plot_data$shift)) {
        cat("The 'shift' column must be numeric, found:\n")
        print(head(plot_data$shift))
    }

    # Convert selected_column_tf_SameStrand to a factor with labels
    plot_data$selected_column_tf_SameStrand <- factor(selected_column_tf_SameStrand, levels = c(1, 0), labels = c("same", "opposite"))

    plot_data.orig_num_lines <- nrow(plot_data)
    plot_data <- na.omit(plot_data)
    plot_data.contributing_num_lines <- nrow(plot_data)


    # Create the histogram
    #p <- ggplot(plot_data, aes(x = selected_column_tf_Shift, fill = selected_column_tf_SameStrand, group = selected_column_tf_SameStrand)) +
    p <- ggplot(plot_data, aes(x = shift, fill = selected_column_tf_SameStrand, group = selected_column_tf_SameStrand)) +
         geom_histogram(position = "dodge", binwidth = binwidth ) +
         labs(title = paste("Shift of TF ",prettyIdentifierJaspar(tf)," for ",subset.name," grouped by strand",sep=""),
             x = "Shift", y = "Count", fill = "Strand") +
         theme_minimal() +
         annotate("text", x = Inf, y = Inf, 
             label = paste("N / total =", plot_data.contributing_num_lines, " / ",plot_data.orig_num_lines), 
             hjust = 1.1, vjust = 1.1, 
             size = 5, color = "black")  # Add annotation for number of data points
  
    # Save the plot to a file
    dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(paste(output.dir,paste("histogram_of_shift_for_tf_",tf,"_subset_",subset.name,".png",sep=""),sep="/"),
        plot=p, width=8, height=6, dpi=300)
} 