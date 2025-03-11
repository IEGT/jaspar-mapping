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

## Selection of columns representing CUT&RUN data
## This results in a data frame with columns representing the CUT&RUN data
## with columns representing the individual PSSM and the strand of the motif
## found in the genomic sequence.
## The rows represent the individual binding sites of the TF. The value is
## NA if the motif was not found within the window of 100 bp to either side of the
## TP73 binding site. A value of 0 stands for the oposite strand, 1 for the same strand.
selected_columns.StrandEqual <- m[, ..cols.StrandEqual]

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


# Determination of quantiles 
## These quantiles are used to determine the threshold for the number of reads in the CUT&RUN data
## that are considered to be reliable
quantiles.cutandrun <- sapply(m[, ..cols.cutandrun], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.tp73 <- sapply(m[, ..cols.cutandrun.tp73], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.tp73.saos <- sapply(m[, ..cols.cutandrun.tp73.saos], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.tp73.skmel <- sapply(m[, ..cols.cutandrun.tp73.skmel], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.tp73.TAa <- sapply(m[, ..cols.cutandrun.tp73.TAa], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.tp73.DNb <- sapply(m[, ..cols.cutandrun.tp73.DNb], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.tp73.GFP <- sapply(m[, ..cols.cutandrun.tp73.GFP], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.pos <- sapply(m[, ..cols.cutandrun.pos], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.pos.saos <- sapply(m[, ..cols.cutandrun.pos.saos], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.pos.skmel <- sapply(m[, ..cols.cutandrun.pos.skmel], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.pos.TAa <- sapply(m[, ..cols.cutandrun.pos.TAa], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.pos.DNb <- sapply(m[, ..cols.cutandrun.pos.DNb], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)
quantiles.cutandrun.pos.GFP <- sapply(m[, ..cols.cutandrun.pos.GFP], quantile, probs = c(0,25,50,75,90,95,99,99.5,99.9,100)/100)

# Sum of activities in CUT&RUN data
sum.cutandrun.tp73 <- rowSums(m[, ..cols.cutandrun.tp73])
sum.cutandrun.tp73.saos <- rowSums(m[, ..cols.cutandrun.tp73.saos])
sum.cutandrun.tp73.skmel <- rowSums(m[, ..cols.cutandrun.tp73.skmel])
sum.cutandrun.tp73.TAa <- rowSums(m[, ..cols.cutandrun.tp73.TAa])
sum.cutandrun.tp73.DNb <- rowSums(m[, ..cols.cutandrun.tp73.DNb])
sum.cutandrun.tp73.GFP <- rowSums(m[, ..cols.cutandrun.tp73.GFP])
sum.cutandrun.pos <- rowSums(m[, ..cols.cutandrun.pos])
sum.cutandrun.pos.saos <- rowSums(m[, ..cols.cutandrun.pos.saos])
sum.cutandrun.pos.skmel <- rowSums(m[, ..cols.cutandrun.pos.skmel])
sum.cutandrun.pos.TAa <- rowSums(m[, ..cols.cutandrun.pos.TAa])
sum.cutandrun.pos.DNb <- rowSums(m[, ..cols.cutandrun.pos.DNb])
sum.cutandrun.pos.GFP <- rowSums(m[, ..cols.cutandrun.pos.GFP])

gc()

sort(sum.NumInWindow)
cat("I: Quantiles for individual CUT&RUN experiments:\n")
print(quantiles.cutandrun)
cat("I: Quantiles for sum across all CUT&RUN experiments for p73-antibody:\n")
print(quantile(sum.cutandrun.tp73,probs=c(0,25,50,75,90,95,99,99.5,99.9,100)/100) )
cat("I: Quantiles for sum across all CUT&RUN experiments for p73-antibody in Saos cells:\n")
print(quantile(sum.cutandrun.tp73.saos,probs=c(0,25,50,75,90,95,99,99.5,99.9,100)/100) )
cat("I: Quantiles for sum across all CUT&RUN experiments for p73-antibody in SKmel-29 cells:\n")
print(quantile(sum.cutandrun.tp73.skmel,probs=c(0,25,50,75,90,95,99,99.5,99.9,100)/100) )
cat("I: Quantiles for sum across all CUT&RUN experiments for TAp73alpha:\n")
print(quantile(sum.cutandrun.tp73.TAa,probs=c(0,25,50,75,90,95,99,99.5,99.9,100)/100) )
cat("I: Quantiles for sum across all CUT&RUN experiments for DNp73beta:\n")
print(quantile(sum.cutandrun.tp73.DNb,probs=c(0,25,50,75,90,95,99,99.5,99.9,100)/100) )
cat("I: Quantiles for sum across all CUT&RUN experiments for GFP-only:\n")
print(quantile(sum.cutandrun.tp73.GFP,probs=c(0,25,50,75,90,95,99,99.5,99.9,100)/100) )


# The 50th, 75th, 90th and 99th percentile for the number of reads in TFBS for different transfected isoforms
sum.cutandrun.tp73.quantile.50 <- max(quantile(sum.cutandrun.tp73,probs=0.50),0.5) # 1
sum.cutandrun.tp73.quantile.75 <- max(quantile(sum.cutandrun.tp73,probs=0.75),0.5) # 2
sum.cutandrun.tp73.quantile.90 <- quantile(sum.cutandrun.tp73,probs=0.90) # 4
sum.cutandrun.tp73.quantile.95 <- quantile(sum.cutandrun.tp73,probs=0.95) # 6
sum.cutandrun.tp73.quantile.99 <- quantile(sum.cutandrun.tp73,probs=0.99) # 20
sum.cutandrun.tp73.TAa.quantile.50 <- max(quantile(sum.cutandrun.tp73.TAa,probs=0.50),0.5) # 1
sum.cutandrun.tp73.DNb.quantile.50 <- max(quantile(sum.cutandrun.tp73.DNb,probs=0.50),0.5) # 1
sum.cutandrun.tp73.GFP.quantile.50 <- max(quantile(sum.cutandrun.tp73.GFP,probs=0.50),0.5) # 1
sum.cutandrun.tp73.TAa.quantile.75 <- quantile(sum.cutandrun.tp73.TAa,probs=0.75) # 1
sum.cutandrun.tp73.DNb.quantile.75 <- quantile(sum.cutandrun.tp73.DNb,probs=0.75) # 1
sum.cutandrun.tp73.GFP.quantile.75 <- max(quantile(sum.cutandrun.tp73.GFP,probs=0.75), 0.5) # 0
sum.cutandrun.tp73.TAa.quantile.90 <- quantile(sum.cutandrun.tp73.TAa,probs=0.90) # 2
sum.cutandrun.tp73.DNb.quantile.90 <- quantile(sum.cutandrun.tp73.DNb,probs=0.90) # 2
sum.cutandrun.tp73.GFP.quantile.90 <- quantile(sum.cutandrun.tp73.GFP,probs=0.90) # 1
sum.cutandrun.tp73.TAa.quantile.95 <- quantile(sum.cutandrun.tp73.TAa,probs=0.95) # 3
sum.cutandrun.tp73.DNb.quantile.95 <- quantile(sum.cutandrun.tp73.DNb,probs=0.95) # 3
sum.cutandrun.tp73.GFP.quantile.95 <- quantile(sum.cutandrun.tp73.GFP,probs=0.95) # 1
sum.cutandrun.tp73.TAa.quantile.99 <- quantile(sum.cutandrun.tp73.TAa,probs=0.99) # 13
sum.cutandrun.tp73.DNb.quantile.99 <- quantile(sum.cutandrun.tp73.DNb,probs=0.99) # 9
sum.cutandrun.tp73.GFP.quantile.99 <- quantile(sum.cutandrun.tp73.GFP,probs=0.99) # 2
sum.cutandrun.tp73.skmel29_2.TAa.quantile.50 <- quantile(m$"tp73_skmel29_2_TA",probs=0.50)
sum.cutandrun.tp73.skmel29_2.DNb.quantile.50 <- quantile(m$"tp73_skmel29_2_DN",probs=0.50)
sum.cutandrun.tp73.skmel29_2.GFP.quantile.50 <- quantile(m$"tp73_skmel29_2_GFP",probs=0.50)
sum.cutandrun.tp73.skmel29_2.TAa.quantile.75 <- quantile(m$"tp73_skmel29_2_TA",probs=0.75)
sum.cutandrun.tp73.skmel29_2.DNb.quantile.75 <- quantile(m$"tp73_skmel29_2_DN",probs=0.75)
sum.cutandrun.tp73.skmel29_2.GFP.quantile.75 <- quantile(m$"tp73_skmel29_2_GFP",probs=0.75)
sum.cutandrun.tp73.skmel29_2.TAa.quantile.90 <- quantile(m$"tp73_skmel29_2_TA",probs=0.90)
sum.cutandrun.tp73.skmel29_2.DNb.quantile.90 <- quantile(m$"tp73_skmel29_2_DN",probs=0.90)
sum.cutandrun.tp73.skmel29_2.GFP.quantile.90 <- quantile(m$"tp73_skmel29_2_GFP",probs=0.90)
sum.cutandrun.tp73.skmel29_2.TAa.quantile.95 <- quantile(m$"tp73_skmel29_2_TA",probs=0.95)
sum.cutandrun.tp73.skmel29_2.DNb.quantile.95 <- quantile(m$"tp73_skmel29_2_DN",probs=0.95)
sum.cutandrun.tp73.skmel29_2.GFP.quantile.95 <- quantile(m$"tp73_skmel29_2_GFP",probs=0.95)
sum.cutandrun.tp73.skmel29_2.TAa.quantile.99 <- quantile(m$"tp73_skmel29_2_TA",probs=0.99)
sum.cutandrun.tp73.skmel29_2.DNb.quantile.99 <- quantile(m$"tp73_skmel29_2_DN",probs=0.99)
sum.cutandrun.tp73.skmel29_2.GFP.quantile.99 <- quantile(m$"tp73_skmel29_2_GFP",probs=0.99)
sum.cutandrun.tp73.saos2.TAa.quantile.50 <- quantile(m$"tp73_saos2_TA",probs=0.50)
sum.cutandrun.tp73.saos2.DNb.quantile.50 <- quantile(m$"tp73_saos2_DN",probs=0.50)
sum.cutandrun.tp73.saos2.GFP.quantile.50 <- quantile(m$"tp73_saos2_GFP",probs=0.50)
sum.cutandrun.tp73.saos2.TAa.quantile.75 <- quantile(m$"tp73_saos2_TA",probs=0.75)
sum.cutandrun.tp73.saos2.DNb.quantile.75 <- quantile(m$"tp73_saos2_DN",probs=0.75)
sum.cutandrun.tp73.saos2.GFP.quantile.75 <- quantile(m$"tp73_saos2_GFP",probs=0.75)
sum.cutandrun.tp73.saos2.TAa.quantile.90 <- quantile(m$"tp73_saos2_TA",probs=0.90)
sum.cutandrun.tp73.saos2.DNb.quantile.90 <- quantile(m$"tp73_saos2_DN",probs=0.90)
sum.cutandrun.tp73.saos2.GFP.quantile.90 <- quantile(m$"tp73_saos2_GFP",probs=0.90)
sum.cutandrun.tp73.saos2.TAa.quantile.95 <- quantile(m$"tp73_saos2_TA",probs=0.95)
sum.cutandrun.tp73.saos2.DNb.quantile.95 <- quantile(m$"tp73_saos2_DN",probs=0.95)
sum.cutandrun.tp73.saos2.GFP.quantile.95 <- quantile(m$"tp73_saos2_GFP",probs=0.95)
sum.cutandrun.tp73.saos2.TAa.quantile.99 <- quantile(m$"tp73_saos2_TA",probs=0.99)
sum.cutandrun.tp73.saos2.DNb.quantile.99 <- quantile(m$"tp73_saos2_DN",probs=0.99)
sum.cutandrun.tp73.saos2.GFP.quantile.99 <- quantile(m$"tp73_saos2_GFP",probs=0.99)

# analogously for pos
sum.cutandrun.pos.quantile.50 <- max(quantile(sum.cutandrun.pos,probs=0.50),0.5) # 0
sum.cutandrun.pos.quantile.75 <- max(quantile(sum.cutandrun.pos,probs=0.75),0.5) # 1
sum.cutandrun.pos.quantile.90 <- quantile(sum.cutandrun.pos,probs=0.90) # 3
sum.cutandrun.pos.quantile.95 <- quantile(sum.cutandrun.pos,probs=0.95) # 6
sum.cutandrun.pos.quantile.99 <- quantile(sum.cutandrun.pos,probs=0.99) # 80
sum.cutandrun.pos.TAa.quantile.50 <- max(quantile(sum.cutandrun.pos.TAa,probs=0.50),0.5) # 1
sum.cutandrun.pos.DNb.quantile.50 <- max(quantile(sum.cutandrun.pos.DNb,probs=0.50),0.5) # 0
sum.cutandrun.pos.GFP.quantile.50 <- max(quantile(sum.cutandrun.pos.GFP,probs=0.50),0.5) # 1
sum.cutandrun.pos.TAa.quantile.75 <- quantile(sum.cutandrun.pos.TAa,probs=0.75) # 1
sum.cutandrun.pos.DNb.quantile.75 <- quantile(sum.cutandrun.pos.DNb,probs=0.75) # 1
sum.cutandrun.pos.GFP.quantile.75 <- quantile(sum.cutandrun.pos.GFP,probs=0.75) # 1
sum.cutandrun.pos.TAa.quantile.90 <- quantile(sum.cutandrun.pos.TAa,probs=0.90) # 1
sum.cutandrun.pos.DNb.quantile.90 <- quantile(sum.cutandrun.pos.DNb,probs=0.90) # 1
sum.cutandrun.pos.GFP.quantile.90 <- quantile(sum.cutandrun.pos.GFP,probs=0.90) # 1
sum.cutandrun.pos.TAa.quantile.95 <- quantile(sum.cutandrun.pos.TAa,probs=0.95) # 2
sum.cutandrun.pos.DNb.quantile.95 <- quantile(sum.cutandrun.pos.DNb,probs=0.95) # 2
sum.cutandrun.pos.GFP.quantile.95 <- quantile(sum.cutandrun.pos.GFP,probs=0.95) # 2
sum.cutandrun.pos.TAa.quantile.99 <- quantile(sum.cutandrun.pos.TAa,probs=0.99) # 20
sum.cutandrun.pos.DNb.quantile.99 <- quantile(sum.cutandrun.pos.DNb,probs=0.99) # 30
sum.cutandrun.pos.GFP.quantile.99 <- quantile(sum.cutandrun.pos.GFP,probs=0.99) # 28

# Clean up memory
gc(full=T)

# Having encountered frequent crashes, maybe 
# save.image()
# at this point.

# Association of CUT&RUN data with PSSM binding sites
## The number of PSSM matches that are in the top 1%,10%,25%,50% of CUT&RUN data
sum.NumInWindow.equal.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 == 0]))
sum.NumInWindow.filtered.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 >= sum.cutandrun.tp73.quantile.50]))
sum.NumInWindow.filtered.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 >= sum.cutandrun.tp73.quantile.75 ]))
sum.NumInWindow.filtered.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 >= sum.cutandrun.tp73.quantile.90 ]))
sum.NumInWindow.filtered.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 >= sum.cutandrun.tp73.quantile.95 ]))
sum.NumInWindow.filtered.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 >= sum.cutandrun.tp73.quantile.99 ]))
cat("head fitered@50\n")
print(head(sort(sum.NumInWindow.filtered.50),5))
cat("head fitered@75\n")
print(head(sort(sum.NumInWindow.filtered.75),5))
cat("head fitered@90\n")
print(head(sort(sum.NumInWindow.filtered.90),5))
cat("head fitered@95\n")
print(head(sort(sum.NumInWindow.filtered.95),5))
cat("head fitered@99\n")
print(head(sort(sum.NumInWindow.filtered.99),5))

cat("tail fitered@50\n")
print(tail(sort(sum.NumInWindow.filtered.50),5))
cat("tail fitered@75\n")
print(tail(sort(sum.NumInWindow.filtered.75),5))
cat("tail fitered@90\n")
print(tail(sort(sum.NumInWindow.filtered.90),5))
cat("tail fitered@95\n")
print(tail(sort(sum.NumInWindow.filtered.95),5))
cat("tail fitered@99\n")
print(tail(sort(sum.NumInWindow.filtered.99),5))
gc()

ratio.tp73.50 <- (sum.NumInWindow.filtered.50/sum(sum.cutandrun.tp73>=sum.cutandrun.tp73.quantile.50)) / (sum.NumInWindow.equal.0/sum(sum.cutandrun.tp73==0))
ratio.tp73.75 <- (sum.NumInWindow.filtered.75/sum(sum.cutandrun.tp73>=sum.cutandrun.tp73.quantile.75)) / (sum.NumInWindow.equal.0/sum(sum.cutandrun.tp73==0))
ratio.tp73.90 <- (sum.NumInWindow.filtered.90/sum(sum.cutandrun.tp73>=sum.cutandrun.tp73.quantile.90)) / (sum.NumInWindow.equal.0/sum(sum.cutandrun.tp73==0))
ratio.tp73.95 <- (sum.NumInWindow.filtered.95/sum(sum.cutandrun.tp73>=sum.cutandrun.tp73.quantile.95)) / (sum.NumInWindow.equal.0/sum(sum.cutandrun.tp73==0))

ratio.tp73.99 <- (
                sum.NumInWindow.filtered.99     # Number of predicted occurences of TF in TFBS that are in top 1% of signals of TP73 cut'n'run
                 /
                sum(sum.cutandrun.tp73>=sum.cutandrun.tp73.quantile.99) # Number of TFBS that constitute the top 1%
            ) / (
                sum.NumInWindow.equal.0         # Number of predicted occurences of TF in TFBS that have no coverage in TP73 cut'n'run
                /
                sum(sum.cutandrun.tp73==0)      # Number of in TFBS that have no coverage in TP73 cut'n'run
            )

ratio.pos.50 <- (sum.NumInWindow.filtered.50/sum(sum.cutandrun.pos >= sum.cutandrun.pos.quantile.50)) / (sum.NumInWindow.equal.0/sum(sum.cutandrun.pos==0))
ratio.pos.75 <- (sum.NumInWindow.filtered.75/sum(sum.cutandrun.pos >= sum.cutandrun.pos.quantile.75)) / (sum.NumInWindow.equal.0/sum(sum.cutandrun.pos==0))
ratio.pos.90 <- (sum.NumInWindow.filtered.90/sum(sum.cutandrun.pos >= sum.cutandrun.pos.quantile.90)) / (sum.NumInWindow.equal.0/sum(sum.cutandrun.pos==0))
ratio.pos.95 <- (sum.NumInWindow.filtered.95/sum(sum.cutandrun.pos >= sum.cutandrun.pos.quantile.95)) / (sum.NumInWindow.equal.0/sum(sum.cutandrun.pos==0))
ratio.pos.99 <- (sum.NumInWindow.filtered.99/sum(sum.cutandrun.pos >= sum.cutandrun.pos.quantile.99)) / (sum.NumInWindow.equal.0/sum(sum.cutandrun.pos==0))

gc()

cat("head tp73 ratio@50\n")
print(head(sort(ratio.tp73.50),5))
cat("tail tp73 ratio@50\n")
print(tail(sort(ratio.tp73.50),5))
cat("head tp73 ratio@75\n")
print(head(sort(ratio.tp73.75),5))
cat("tail tp73 ratio@75\n")
print(tail(sort(ratio.tp73.75),5))
cat("head tp73 ratio@90\n")
print(head(sort(ratio.tp73.90),5))
cat("tail tp73 ratio@90\n")
print(tail(sort(ratio.tp73.90),5))
cat("head tp73 ratio@99\n")
print(head(sort(ratio.tp73.99),5))
cat("tail tp73 ratio@99\n")
print(tail(sort(ratio.tp73.99),5))
cat("head pos ratio@50\n")
print(head(sort(ratio.pos.50),5))
cat("tail pos ratio@50\n")
print(tail(sort(ratio.pos.50),5))
cat("head pos ratio@75\n")
print(head(sort(ratio.pos.75),5))
cat("tail pos ratio@75\n")
print(tail(sort(ratio.pos.75),5))
cat("head pos ratio@90\n")
print(head(sort(ratio.pos.90),5))
cat("tail pos ratio@90\n")
print(tail(sort(ratio.pos.90),5))
cat("head pos ratio@99\n")
print(head(sort(ratio.pos.99),5))
cat("tail pos ratio@99\n")
print(tail(sort(ratio.pos.99),5))
gc()

pretty.table <- function(x) {
    x.sort <- sort(x)
    d<-data.frame(sapply(strsplit(x=names(x.sort),split="_"),function(X) paste(X[1]," (",X[2],")",sep="")), round(x.sort,3))
    rownames(d) <- NULL
    colnames(d) <- c("TF","Ratio")
    d
}

pretty.table(ratio.tp73.50)
pretty.table(ratio.tp73.75)
pretty.table(ratio.tp73.90)
pretty.table(ratio.tp73.95)
pretty.table(ratio.tp73.99)
pretty.table(ratio.pos.50)
pretty.table(ratio.pos.75)
pretty.table(ratio.pos.90)
pretty.table(ratio.pos.95)
pretty.table(ratio.pos.99)

# Enrichment 

sum.NumInWindow.equal.TAa.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa == 0]))
sum.NumInWindow.filtered.TAa.quantile.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa >= sum.cutandrun.tp73.TAa.quantile.50]))
sum.NumInWindow.filtered.TAa.quantile.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa >= sum.cutandrun.tp73.TAa.quantile.90]))
sum.NumInWindow.filtered.TAa.quantile.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa >= sum.cutandrun.tp73.TAa.quantile.95]))
sum.NumInWindow.filtered.TAa.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa >= sum.cutandrun.tp73.TAa.quantile.99]))
ratio.TAa.quantile.50 <- (sum.NumInWindow.filtered.TAa.quantile.50/sum(sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.50)) / (sum.NumInWindow.equal.TAa.0/sum(sum.cutandrun.tp73.TAa==0))
ratio.TAa.quantile.90 <- (sum.NumInWindow.filtered.TAa.quantile.90/sum(sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.90)) / (sum.NumInWindow.equal.TAa.0/sum(sum.cutandrun.tp73.TAa==0))
ratio.TAa.quantile.95 <- (sum.NumInWindow.filtered.TAa.quantile.95/sum(sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.95)) / (sum.NumInWindow.equal.TAa.0/sum(sum.cutandrun.tp73.TAa==0))
ratio.TAa.quantile.99 <- (sum.NumInWindow.filtered.TAa.quantile.99/sum(sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.99)) / (sum.NumInWindow.equal.TAa.0/sum(sum.cutandrun.tp73.TAa==0))
sum.NumInWindow.equal.DNb.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa == 0]))
sum.NumInWindow.filtered.DNb.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.DNb >= sum.cutandrun.tp73.DNb.quantile.50]))
sum.NumInWindow.filtered.DNb.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.DNb >= sum.cutandrun.tp73.DNb.quantile.90]))
sum.NumInWindow.filtered.DNb.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.DNb >= sum.cutandrun.tp73.DNb.quantile.95]))
sum.NumInWindow.filtered.DNb.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.DNb >= sum.cutandrun.tp73.DNb.quantile.99]))
ratio.DNb.quantile.50 <- (sum.NumInWindow.filtered.DNb.50/sum(sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.50)) / (sum.NumInWindow.equal.DNb.0/sum(sum.cutandrun.tp73.DNb==0))
ratio.DNb.quantile.90 <- (sum.NumInWindow.filtered.DNb.99/sum(sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.90)) / (sum.NumInWindow.equal.DNb.0/sum(sum.cutandrun.tp73.DNb==0))
ratio.DNb.quantile.95 <- (sum.NumInWindow.filtered.DNb.95/sum(sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.95)) / (sum.NumInWindow.equal.DNb.0/sum(sum.cutandrun.tp73.DNb==0))
ratio.DNb.quantile.99 <- (sum.NumInWindow.filtered.DNb.99/sum(sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.99)) / (sum.NumInWindow.equal.DNb.0/sum(sum.cutandrun.tp73.DNb==0))

head(pretty.table(ratio.TAa.quantile.50))
tail(pretty.table(ratio.TAa.quantile.50))
head(pretty.table(ratio.DNb.quantile.50))
tail(pretty.table(ratio.DNb.quantile.50))
head(pretty.table(ratio.TAa.quantile.99))
tail(pretty.table(ratio.TAa.quantile.99))
head(pretty.table(ratio.DNb.quantile.99))
tail(pretty.table(ratio.DNb.quantile.99))

ratio.TAa.quantile.50.vs.ratio.DNb.quantile.50 <- ratio.TAa.quantile.50 / ratio.DNb.quantile.50
head(pretty.table(ratio.TAa.quantile.50.vs.ratio.DNb.quantile.50))
tail(pretty.table(ratio.TAa.quantile.50.vs.ratio.DNb.quantile.50))

ratio.TAa.quantile.99.vs.ratio.DNb.quantile.99 <- ratio.TAa.quantile.99 / ratio.DNb.quantile.99
head(pretty.table(ratio.TAa.quantile.99.vs.ratio.DNb.quantile.99))
# Auffälligkeit p53
tail(pretty.table(ratio.TAa.quantile.99.vs.ratio.DNb.quantile.99))




sum.NumInWindow.equal.saos2.TAa.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_TA" == 0]))
sum.NumInWindow.filtered.saos2.TAa.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_TA" >= sum.cutandrun.tp73.saos2.TAa.quantile.99]))
                 # Number of TFBS for particular TF in top 1% of CUT&RUN confirmed sites
                                                       # Number of TFBS in top 1% of confirmed sites
ratio.saos2.TAa.quantile.99 <- (sum.NumInWindow.filtered.saos2.TAa.quantile.99/sum(m$"tp73_saos2_TA">=sum.cutandrun.tp73.saos2.TAa.quantile.99)) / (sum.NumInWindow.equal.saos2.TAa.0/sum(m$"tp73_saos2_TA"==0))

sum.NumInWindow.equal.saos2.DNb.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_DN" == 0]))
sum.NumInWindow.filtered.saos2.DNb.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_DN" >= sum.cutandrun.tp73.saos2.DNb.quantile.99]))
ratio.saos2.DNb.quantile.99 <- (sum.NumInWindow.filtered.saos2.DNb.quantile.99/sum(m$"tp73_saos2_DN">=sum.cutandrun.tp73.saos2.DNb.quantile.99)) / (sum.NumInWindow.equal.saos2.DNb.0/sum(m$"tp73_saos2_DN"==0))

sum.NumInWindow.equal.saos2.GFP.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_GFP" == 0]))
sum.NumInWindow.filtered.saos2.GFP.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_GFP" >= sum.cutandrun.tp73.saos2.GFP.quantile.99]))
ratio.saos2.GFP.quantile.99 <- (sum.NumInWindow.filtered.saos2.GFP.quantile.99/sum(m$"tp73_saos2_GFP">=sum.cutandrun.tp73.saos2.GFP.quantile.99)) / (sum.NumInWindow.equal.saos2.GFP.0/sum(m$"tp73_saos2_GFP"==0))

head(pretty.table(ratio.saos2.TAa.quantile.99))
tail(pretty.table(ratio.saos2.TAa.quantile.99))
head(pretty.table(ratio.saos2.DNb.quantile.99))
tail(pretty.table(ratio.saos2.DNb.quantile.99))
head(pretty.table(ratio.saos2.GFP.quantile.99))
tail(pretty.table(ratio.saos2.GFP.quantile.99))

ratio.saos2.TAa.99.vs.ratio.saos2.DNb.99 <- ratio.saos2.TAa.99 / ratio.saos2.DNb.99
pretty.table(ratio.saos2.TAa.99.vs.ratio.saos2.DNb.99)
ratio.saos2.TAa.99.vs.ratio.saos2.GFP.99 <- ratio.saos2.TAa.99 / ratio.saos2.GFP.99
pretty.table(ratio.saos2.TAa.99.vs.ratio.saos2.GFP.99)
ratio.saos2.DNb.99.vs.ratio.saos2.GFP.99 <- ratio.saos2.DNb.99 / ratio.saos2.GFP.99
pretty.table(ratio.saos2.DNb.99.vs.ratio.saos2.GFP.99)


sum.NumInWindow.equal.skmel29_2.TAa.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_TA" == 0]))
sum.NumInWindow.filtered.skmel29_2.TAa.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_TA" >= sum.cutandrun.tp73.skmel29_2.TAa.quantile.99]))
                 # Number of TFBS for particular TF in top 1% of CUT&RUN confirmed sites
                                                       # Number of TFBS in top 1% of confirmed sites
ratio.skmel29_2.TAa.quantile.99 <- (sum.NumInWindow.filtered.skmel29_2.TAa.quantile.99/sum(m$"tp73_skmel29_2_TA">=sum.cutandrun.tp73.skmel29_2.TAa.quantile.99)) / (sum.NumInWindow.equal.skmel29_2.TAa.0/sum(m$"tp73_skmel29_2_TA"==0))

sum.NumInWindow.equal.skmel29_2.DNb.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_DN" == 0]))
sum.NumInWindow.filtered.skmel29_2.DNb.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_DN" >= sum.cutandrun.tp73.skmel29_2.DNb.quantile.99]))
ratio.skmel29_2.DNb.quantile.99 <- (sum.NumInWindow.filtered.skmel29_2.DNb.quantile.99/sum(m$"tp73_skmel29_2_DN">=sum.cutandrun.tp73.skmel29_2.DNb.quantile.99)) / (sum.NumInWindow.equal.skmel29_2.DNb.0/sum(m$"tp73_skmel29_2_DN"==0))

sum.NumInWindow.equal.skmel29_2.GFP.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_GFP" == 0]))
sum.NumInWindow.filtered.skmel29_2.GFP.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_GFP" >= sum.cutandrun.tp73.skmel29_2.GFP.quantile.99]))
ratio.skmel29_2.GFP.quantile.99 <- (sum.NumInWindow.filtered.skmel29_2.GFP.quantile.99/sum(m$"tp73_skmel29_2_GFP">=sum.cutandrun.tp73.skmel29_2.GFP.quantile.99)) / (sum.NumInWindow.equal.skmel29_2.GFP.0/sum(m$"tp73_skmel29_2_GFP"==0))

head(pretty.table(ratio.skmel29_2.TAa.quantile.99),10)
tail(pretty.table(ratio.skmel29_2.TAa.quantile.99),10)
head(pretty.table(ratio.skmel29_2.DNb.quantile.99),10)
tail(pretty.table(ratio.skmel29_2.DNb.quantile.99),10)
head(pretty.table(ratio.skmel29_2.GFP.quantile.99),10)
tail(pretty.table(ratio.skmel29_2.GFP.quantile.99),10)

ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.DNb.quantile.99 <- ratio.skmel29_2.TAa.quantile.99 / ratio.skmel29_2.DNb.quantile.99
head(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.DNb.quantile.99),10)
tail(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.DNb.quantile.99),10)
ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99 <- ratio.skmel29_2.TAa.quantile.99 / ratio.skmel29_2.GFP.quantile.99
head(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99),10)
tail(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99),10)
ratio.skmel29_2.DNb.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99 <- ratio.skmel29_2.DNb.quantile.99 / ratio.skmel29_2.GFP.quantile.99
head(pretty.table(ratio.skmel29_2.DNb.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99),10)
tail(pretty.table(ratio.skmel29_2.DNb.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99),10)


ratio.skmel29_2.TAa.quantile.99.vs.ratio.saos2.TAa.quantile.99 <- ratio.skmel29_2.TAa.quantile.99 / ratio.saos2.TAa.quantile.99
head(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.saos2.TAa.quantile.99),10)
tail(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.saos2.TAa.quantile.99),10)
ratio.skmel29_2.DNb.quantile.99.vs.ratio.saos2.DNb.quantile.99 <- ratio.skmel29_2.DNb.quantile.99 / ratio.saos2.DNb.quantile.99
head(pretty.table(ratio.skmel29_2.DNb.quantile.99.vs.ratio.saos2.DNb.quantile.99),10)
tail(pretty.table(ratio.skmel29_2.DNb.quantile.99.vs.ratio.saos2.DNb.quantile.99),10)
ratio.skmel29_2.GFP.quantile.99.vs.ratio.saos2.GFP.quantile.99 <- ratio.skmel29_2.GFP.quantile.99 / ratio.saos2.GFP.quantile.99
head(pretty.table(ratio.skmel29_2.GFP.quantile.99.vs.ratio.saos2.GFP.quantile.99),10)
tail(pretty.table(ratio.skmel29_2.GFP.quantile.99.vs.ratio.saos2.GFP.quantile.99),10)


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
                               tfbs.selection=NULL,
                               subset.name="unknown",
                               binwidth=2

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
    ggsave(paste("histogram_of_shift_for_tf_",tf,"_subset_",subset.name,".png",sep=""),
        plot=p, width=8, height=6, dpi=300)
} 

l <- list("TAa"=sum.cutandrun.tp73.TAa, "DNb"=sum.cutandrun.tp73.DNb, "GFP"=sum.cutandrun.tp73_GFP)

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
