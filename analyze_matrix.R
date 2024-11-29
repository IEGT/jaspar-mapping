library(data.table)
 m <- fread("TP73_datatable.bed.gz")
 print(colnames(m))

 cols.NumInWindow <- grepl("_NumInWindow", colnames(m))
 selected_columns <- m[, cols.NumInWindow]
 quantiles.NumInWindow <- sapply(selected_columns, quantile, probs = c(0,0.25, 0.5, 0.75,1))
 sum.NumInWindow <- sapply(selected_columns, sum)

 cols.Score <- grepl("_Score", colnames(m))
 selected_columns <- m[, cols.Score]
 quantiles.Score <- sapply(selected_columns, quantile, probs = c(0,0.25, 0.5, 0.75,1))
 sum.Score <- sapply(selected_columns, sum, na.rm=T)

 cols.StrandEqual <- grepl("_StrandEqual", colnames(m))
 selected_columns.StrandEqual <- m[, cols.StrandEqual]

cols.cutandrun <- c("pos_saos2_DN", "pos_saos2_GFP", "pos_saos2_TA", "pos_skmel29_2_DN",
  "pos_skmel29_2_GFP", "pos_skmel29_2_TA", "tp73_saos2_DN", "tp73_saos2_GFP",
  "tp73_saos2_TA", "tp73_skmel29_2_DN", "tp73_skmel29_2_GFP", "tp73_skmel29_2_TA")
cols.cutandrun.tp73 <- c("tp73_saos2_DN", "tp73_saos2_GFP", "tp73_saos2_TA", "tp73_skmel29_2_DN",
  "tp73_skmel29_2_GFP", "tp73_skmel29_2_TA")

selected_columns_cutnrun <- m[, cols.cutandrun]
selected_columns_cutnrun_tp73 <- m[, cols.cutandrun.tp73]
quantiles_cutandrun <- sapply(selected_columns_cutnrun, quantile, probs = c(75,90,95,99,99.5,99.9,100)/100)
sum_cutandrun_tp73 <- apply(as.matrix(selected_columns_cutnrun_tp73), 1, sum)
sort(sum.NumInWindow)

sum.NumInWindow.equal.0 <- sapply(selected_columns, function(X) sum(X[sum_cutandrun_tp73 == 0]))
sum.NumInWindow.filtered.0 <- sapply(selected_columns, function(X) sum(X[sum_cutandrun_tp73 > 0]))
sum.NumInWindow.filtered.20 <- sapply(selected_columns, function(X) sum(X[sum_cutandrun_tp73 > 20]))
sort(sum.NumInWindow.filtered.0)
sort(sum.NumInWindow.filtered.20)

ratio.0 <- (sum.NumInWindow.filtered.0/sum(sum_cutandrun_tp73>0)) / (sum.NumInWindow.equal.0/sum(sum_cutandrun_tp73==0))
ratio.20 <- (sum.NumInWindow.filtered.20/sum(sum_cutandrun_tp73>20)) / (sum.NumInWindow.equal.0/sum(sum_cutandrun_tp73==0))

sort(ratio.0)
sort(ratio.20)

pretty.table <- function(x) {
    x.sort <- sort(x)
    d<-data.frame(sapply(strsplit(x=names(x.sort),split="_"),function(X) paste(X[1]," (",X[2],")",sep="")), round(x.sort,3))
    rownames(d) <- NULL
    colnames(d) <- c("TF","Ratio")
    d
}

pretty.table(ratio.0)
pretty.table(ratio.20)


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

plotShiftOfBinding <- function(tf="TP73_MA0861.1",
                               tfbs.selection=(sum_cutandrun_tp73 > 20),
                               binwidth=2

) {
    selected_column_tf_Score <- m[tfbs.selection, paste(tf,"Score",sep="_")]
    selected_column_tf_Shift <- m[tfbs.selection, paste(tf,"Shift",sep="_")]
    selected_column_tf_SameStrand <- m[tfbs.selection, paste(tf,"StrandEqual",sep="_")]

    data <- data.frame(
    shift=selected_column_tf_Shift,
    strand=selected_column_tf_SameStrand
    )
    # Convert selected_column_tf_SameStrand to a factor with labels
    data$selected_column_tf_SameStrand <- factor(data$strand, levels = c(1, 0), labels = c("same", "opposite"))

    # Create the histogram
    ggplot(data, aes(x = selected_column_tf_Shift, fill = selected_column_tf_SameStrand,
                    group = selected_column_tf_SameStrand)) +
    geom_histogram(position = "dodge", binwidth = binwidth ) +
    labs(title = paste("Shift of TF ",prettyIdentifierJaspar(tf)," grouped by strand",sep=""),
        x = "Shift",
        y = "Count",
        fill = "Strand") +
    theme_minimal()
  
    # Save the plot to a file
    ggsave(paste("histogram_of_shift_for_tf_",tf,"_Shift.png",sep=""))
} 

for(tf in c("TP73_MA0861.1","TP63_MA0525.2","TP53_MA0106.3",
            "E2F2_MA0864.2","RELA_MA0107.1","PLAG1_MA0163.1",
            "SP1_MA0079.5","PLAG1_MA0163.1","REST_MA0138.2","CTCF_MA0139.1",
            "KLF14_MA0740.2","Klf15_MA1890.1","REL_MA0101.1","YY1-2_MA1927.1","TBP_MA0108.2",
            "POU2F1--SOX2_MA1962.1")) {
    cat("I: plotting '",tf,"': ",sep="")
    plotShiftOfBinding(tf)
}
