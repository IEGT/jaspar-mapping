#!/usr/bin/R

# Script to interpret TF binding sites for their association with CUT&RUN data.

options(width=180)

# Columns representing results from CUT&RUN data
## The names of the columns that represent the CUT&RUN data for the positions at which
## TP73 is binding.
cols.cutandrun <- c("pos_saos2_DN", "pos_saos2_GFP", "pos_saos2_TA", "pos_skmel29_2_DN",  "pos_skmel29_2_GFP", "pos_skmel29_2_TA",
 "tp73_saos2_DN", "tp73_saos2_GFP", "tp73_saos2_TA", "tp73_skmel29_2_DN", "tp73_skmel29_2_GFP", "tp73_skmel29_2_TA")
# Subset of columns representing results for TP73 binding in CUT&RUN data
cols.cutandrun.tp73 <- c("tp73_saos2_DN", "tp73_saos2_GFP", "tp73_saos2_TA", "tp73_skmel29_2_DN", "tp73_skmel29_2_GFP", "tp73_skmel29_2_TA")
cols.cutandrun.tp73.saos <- c("tp73_saos2_DN", "tp73_saos2_GFP", "tp73_saos2_TA")
cols.cutandrun.tp73.skmel <- c("tp73_skmel29_2_DN", "tp73_skmel29_2_GFP", "tp73_skmel29_2_TA")
cols.cutandrun.tp73.TAa <- c("tp73_saos2_TA", "tp73_skmel29_2_TA")
cols.cutandrun.tp73.DNb <- c("tp73_saos2_DN", "tp73_skmel29_2_DN")
cols.cutandrun.tp73.GFP <- c("tp73_saos2_GFP", "tp73_skmel29_2_GFP")
cols.cutandrun.tp73.TAa.saos <- c("tp73_saos2_TA")
cols.cutandrun.tp73.DNb.saos <- c("tp73_saos2_DN")
cols.cutandrun.tp73.GFP.saos <- c("tp73_saos2_GFP")
cols.cutandrun.tp73.TAa.skmel <- c("tp73_skmel29_2_TA")
cols.cutandrun.tp73.DNb.skmel <- c("tp73_skmel29_2_DN")
cols.cutandrun.tp73.GFP.skmel <- c("tp73_skmel29_2_GFP")
cols.cutandrun.pos <- c("pos_saos2_DN", "pos_saos2_GFP", "pos_saos2_TA", "pos_skmel29_2_DN", "pos_skmel29_2_GFP", "pos_skmel29_2_TA")
cols.cutandrun.pos.saos <- c("pos_saos2_DN", "pos_saos2_GFP", "pos_saos2_TA")
cols.cutandrun.pos.skmel <- c("pos_skmel29_2_DN", "pos_skmel29_2_GFP", "pos_skmel29_2_TA")
cols.cutandrun.pos.TAa <- c("pos_saos2_TA", "pos_skmel29_2_TA")
cols.cutandrun.pos.DNb <- c("pos_saos2_DN", "pos_skmel29_2_DN")
cols.cutandrun.pos.GFP <- c("pos_saos2_GFP","pos_skmel29_2_GFP")
cols.cutandrun.pos.TAa.saos  <- c("pos_saos2_TA")
cols.cutandrun.pos.DNb.saos  <- c("pos_saos2_DN")
cols.cutandrun.pos.GFP.saos  <- c("pos_saos2_GFP")
cols.cutandrun.pos.TAa.skmel <- c("pos_skmel29_2_TA")
cols.cutandrun.pos.DNb.skmel <- c("pos_skmel29_2_DN")
cols.cutandrun.pos.GFP.skmel <- c("pos_skmel29_2_GFP")


# Import function to read and interpret the matrix for individual chromosomes
source("analyze_matrix_function_lists.R")  

## Ugly?!?

library(ggplot2)
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


prettyIdentifierJaspar <- function(X) {
    A<-strsplit(X,"_")
    B<-sapply(A,function(X) paste(X[1]," (",X[2],")",sep=""))
    B
}

#####
# Creates distance plots for the selected columns

require(ggplot2)

plotShiftOfBinding <- function(tf="TP73_MA0861.1",
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

l <- list("TAa"=sum.cutandrun.tp73.TAa>0, "DNb"=sum.cutandrun.tp73.DNb>0, "GFP"=sum.cutandrun.tp73.GFP>0,
          "TAa_without_DNb"=sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.90 & 
                            sum.cutandrun.tp73.DNb<=sum.cutandrun.tp73.DNb.quantile.50,
          "DNb_without_TAa"=sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.90 & 
                            sum.cutandrun.tp73.TAa<=sum.cutandrun.tp73.TAa.quantile.50)

for (l.name in names(l)) {
    for(tf in c("TP73_MA0861.1","TP63_MA0525.2","TP53_MA0106.3",
            "CTCF_MA0139.1","E2F2_MA0864.2", "FOS_MA1951.1","FOSL1--JUN_MA1129.1",
            "FOXO4_MA0848.1","Foxq1_MA0040.1", "KLF14_MA0740.2","Klf15_MA1890.1","REL_MA0101.1",
            "PLAG1_MA0163.1","POU2F1--SOX2_MA1962.1", "PPARG_MA0066.1",
            "RELA_MA0107.1","REST_MA0138.2", "RXRA--VDR_MA0074.1", "SP1_MA0079.5",
            "TBP_MA0108.2","YY1-2_MA1927.1")) {
        cat("I: plotting '",tf,"' (",l.name,") : ",sep="")
        plotShiftOfBinding(tf=tf,tfbs.selection=l[[l.name]],subset.name=l.name)
        cat("\n")
    }
}
