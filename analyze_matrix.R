#!/usr/bin/R

# Script to interpret TF binding sites for their association with CUT&RUN data.

options(width=180)

# Import of multi-GB large compressed data file
library(data.table)
#m <- fread("TP73_datatable.bed.gz")
m <- fread("TP73_datatable.01.bed.gz")
print(colnames(m))

# Identification of columns of particula type
## NumInWIndow - number of binding sites of particular transcription factor
cols.NumInWindow <- grepl("_NumInWindow", colnames(m))
quantiles.NumInWindow <- sapply(m[, ..cols.NumInWindow], quantile, probs = c(0,0.25, 0.5, 0.75,1))
sum.NumInWindow <- sapply(m[, ..cols.NumInWindow], sum)
gc(full=T)

## Score - computed affinity with which the TF is binding
cols.Score <- grepl("_Score", colnames(m))
quantiles.Score <- sapply(m[, ..cols.Score], quantile, probs = c(0,0.25, 0.5, 0.75,1), na.rm=TRUE)
sum.Score <- sapply(m[, ..cols.Score], sum, na.rm=T)
gc(full=T)

## Identification of strand (same/other) at which the motif was found
cols.StrandEqual <- grepl("_StrandEqual", colnames(m))
selected_columns.StrandEqual <- m[, ..cols.StrandEqual]

# Columns representing results from CUT&RUN data
cols.cutandrun <- c("pos_saos2_DN", "pos_saos2_GFP", "pos_saos2_TA", "pos_skmel29_2_DN",
  "pos_skmel29_2_GFP", "pos_skmel29_2_TA", "tp73_saos2_DN", "tp73_saos2_GFP",
  "tp73_saos2_TA", "tp73_skmel29_2_DN", "tp73_skmel29_2_GFP", "tp73_skmel29_2_TA")
# Subset of columns representing results for TP73 binding in CUT&RUN data
cols.cutandrun.tp73 <- c("tp73_saos2_DN", "tp73_saos2_GFP", "tp73_saos2_TA", "tp73_skmel29_2_DN", "tp73_skmel29_2_GFP", "tp73_skmel29_2_TA")
cols.cutandrun.tp73.TAa <- c("tp73_saos2_TA", "tp73_skmel29_2_TA")
cols.cutandrun.tp73.DNb <- c("tp73_saos2_DN", "tp73_skmel29_2_DN")
cols.cutandrun.tp73.GFP <- c("tp73_saos2_GFP", "tp73_skmel29_2_GFP")


# Determination of quantiles 
quantiles_cutandrun <- sapply(m[, ..cols.cutandrun], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
sum_cutandrun_tp73 <- rowSums(m[, ..cols.cutandrun.tp73])
sum_cutandrun_tp73_TAa <- rowSums(m[, ..cols.cutandrun.tp73.TAa])
sum_cutandrun_tp73_DNb <- rowSums(m[, ..cols.cutandrun.tp73.DNb])
sum_cutandrun_tp73_GFP <- rowSums(m[, ..cols.cutandrun.tp73.GFP])
gc()

sort(sum.NumInWindow)
cat("I: Quantiles for individual CUT&RUN experiments:\n")
print(quantiles_cutandrun)
cat("I: Quantiles for sum across all CUT&RUN experiments for p73:\n")
print(quantile(sum_cutandrun_tp73,probs=c(0,25,50,75,90,95,99,99.5,99.9,100)/100) )
cat("I: Quantiles for sum across all CUT&RUN experiments for TAp73alpha:\n")
print(quantile(sum_cutandrun_tp73_TAa,probs=c(0,25,50,75,90,95,99,99.5,99.9,100)/100) )
cat("I: Quantiles for sum across all CUT&RUN experiments for DNp73beta:\n")
print(quantile(sum_cutandrun_tp73_DNb,probs=c(0,25,50,75,90,95,99,99.5,99.9,100)/100) )
cat("I: Quantiles for sum across all CUT&RUN experiments for GFP-only:\n")
print(quantile(sum_cutandrun_tp73_GFP,probs=c(0,25,50,75,90,95,99,99.5,99.9,100)/100) )

# The 99th percentile for the number of reads in TFBS for different transfected isoforms
sum_cutandrun_tp73.50 <- max(quantile(sum_cutandrun_tp73a,probs=0.50),1) # 1
sum_cutandrun_tp73.99 <- quantile(sum_cutandrun_tp73a,probs=0.99) # 20
sum_cutandrun_tp73_TAa.50 <- max(quantile(sum_cutandrun_tp73_TAa,probs=0.50),1) # 1
sum_cutandrun_tp73_DNb.50 <- max(quantile(sum_cutandrun_tp73_DNb,probs=0.50),1) # 1
sum_cutandrun_tp73_GFP.50 <- max(quantile(sum_cutandrun_tp73_GFP,probs=0.50),1) # 1
sum_cutandrun_tp73_TAa.99 <- quantile(sum_cutandrun_tp73_TAa,probs=0.99) # 13
sum_cutandrun_tp73_DNb.99 <- quantile(sum_cutandrun_tp73_DNb,probs=0.99) # 9
sum_cutandrun_tp73_GFP.99 <- quantile(sum_cutandrun_tp73_GFP,probs=0.99) # 2

# Clean up memory
gc(full=T)

# Having encountered frequent crashes, maybe 
# save.image()
# at this point.

sum.NumInWindow.equal.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum_cutandrun_tp73 == 0]))
sum.NumInWindow.filtered.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum_cutandrun_tp73 >= sum_cutandrun_tp73.50]))
sum.NumInWindow.filtered.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum_cutandrun_tp73 >= sum_cutandrun_tp73.99 ]))
cat("head fitered@50\n")
print(head(sort(sum.NumInWindow.filtered.50)))
cat("tail fitered@50\n")
print(tail(sort(sum.NumInWindow.filtered.50)))
cat("head fitered@99\n")
print(head(sort(sum.NumInWindow.filtered.99)))
cat("tail fitered@99\n")
print(tail(sort(sum.NumInWindow.filtered.99)))
gc()

ratio.50 <- (sum.NumInWindow.filtered.50/sum(sum_cutandrun_tp73>=sum_cutandrun_tp73.50)) / (sum.NumInWindow.equal.0/sum(sum_cutandrun_tp73==0))
ratio.99 <- (
                sum.NumInWindow.filtered.99     # Number of predicted occurences of TF in TFBS that are in top 1% of signals of TP73 cut'n'run
                 /
                sum(sum_cutandrun_tp73>=sum_cutandrun_tp73.99) # Number of TFBS that constitute the top 1%
            ) / (
                sum.NumInWindow.equal.0         # Number of predicted occurences of TF in TFBS that have no coverage in TP73 cut'n'run
                /
                sum(sum_cutandrun_tp73==0)      # Number of in TFBS that have no coverage in TP73 cut'n'run
            )
gc()
cat("head ratio@50\n")
print(head(sort(ratio.50)))
cat("tail ratio@50\n")
print(tail(sort(ratio.50)))
cat("head ratio@20\n")
print(head(sort(ratio.99)))
cat("tail ratio@20\n")
print(tail(sort(ratio.99)))
gc()

pretty.table <- function(x) {
    x.sort <- sort(x)
    d<-data.frame(sapply(strsplit(x=names(x.sort),split="_"),function(X) paste(X[1]," (",X[2],")",sep="")), round(x.sort,3))
    rownames(d) <- NULL
    colnames(d) <- c("TF","Ratio")
    d
}

pretty.table(ratio.50)
pretty.table(ratio.99)

# Enrichment 

sum.NumInWindow.equal.TAa.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum_cutandrun_tp73_TAa == 0]))
sum.NumInWindow.filtered.TAa.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum_cutandrun_tp73_TAa >= sum_cutandrun_tp73_TAa.50]))
sum.NumInWindow.filtered.TAa.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum_cutandrun_tp73_TAa >= sum_cutandrun_tp73_TAa.99]))
ratio.TAa.50 <- (sum.NumInWindow.filtered.TAa.50/sum(sum_cutandrun_tp73_TAa>=sum_cutandrun_tp73_TAa.50)) / (sum.NumInWindow.equal.TAa.0/sum(sum_cutandrun_tp73_TAa==0))
ratio.TAa.99 <- (sum.NumInWindow.filtered.TAa.99/sum(sum_cutandrun_tp73_TAa>=sum_cutandrun_tp73_TAa.99)) / (sum.NumInWindow.equal.TAa.0/sum(sum_cutandrun_tp73_TAa==0))
sum.NumInWindow.equal.DNb.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum_cutandrun_tp73_TAa == 0]))
sum.NumInWindow.filtered.DNb.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum_cutandrun_tp73_DNb >= sum_cutandrun_tp73_DNb.99]))
ratio.DNb.99 <- (sum.NumInWindow.filtered.DNb.99/sum(sum_cutandrun_tp73_DNb>=sum_cutandrun_tp73_DNb.99)) / (sum.NumInWindow.equal.DNb.0/sum(sum_cutandrun_tp73_DNb==0))


pretty.table(ratio.TAa.99)
pretty.table(ratio.DNb.99)

ratio.TAa.99.vs.ratio.DNb.99 <- ratio.TAa.99 / ratio.DNb.99
pretty.table(ratio.TAa.99.vs.ratio.DNb.99)



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
