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

require(data.table)

create.lists.for.chromosome <- function(chromosome="22",reportdir="Reports") {

    filename <- paste("TP73_datatable_",chromosome,".bed.gz",sep="")

    # Import of multi-GB large compressed data file
    #m <- fread("TP73_datatable_1.bed.gz")
    m <- fread(filename,fill=FALSE,showProgress=TRUE)

    m.colnames <- colnames(m)
    m.colnames.basename <- basename(m.colnames)
    m.colnames.dirname <- dirname(m.colnames)
    colnames(m) <- m.colnames.basename

    gc()

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
    sum.NumInWindow.tp73.equal.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 == 0]))
    sum.NumInWindow.tp73.filtered.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 >= sum.cutandrun.tp73.quantile.50]))
    sum.NumInWindow.tp73.filtered.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 >= sum.cutandrun.tp73.quantile.75 ]))
    sum.NumInWindow.tp73.filtered.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 >= sum.cutandrun.tp73.quantile.90 ]))
    sum.NumInWindow.tp73.filtered.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 >= sum.cutandrun.tp73.quantile.95 ]))
    sum.NumInWindow.tp73.filtered.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73 >= sum.cutandrun.tp73.quantile.99 ]))
    sum.NumInWindow.pos.equal.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.pos == 0]))
    sum.NumInWindow.pos.filtered.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.pos >= sum.cutandrun.pos.quantile.50]))
    sum.NumInWindow.pos.filtered.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.pos >= sum.cutandrun.pos.quantile.75 ]))
    sum.NumInWindow.pos.filtered.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.pos >= sum.cutandrun.pos.quantile.90 ]))
    sum.NumInWindow.pos.filtered.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.pos >= sum.cutandrun.pos.quantile.95 ]))
    sum.NumInWindow.pos.filtered.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.pos >= sum.cutandrun.pos.quantile.99 ]))
    cat("head fitered@50\n")
    print(head(sort(sum.NumInWindow.tp73.filtered.50),5))
    cat("head fitered@75\n")
    print(head(sort(sum.NumInWindow.tp73.filtered.75),5))
    cat("head fitered@90\n")
    print(head(sort(sum.NumInWindow.tp73.filtered.90),5))
    cat("head fitered@95\n")
    print(head(sort(sum.NumInWindow.tp73.filtered.95),5))
    cat("head fitered@99\n")
    print(head(sort(sum.NumInWindow.tp73.filtered.99),5))

    cat("tail fitered@50\n")
    print(tail(sort(sum.NumInWindow.tp73.filtered.50),5))
    cat("tail fitered@75\n")
    print(tail(sort(sum.NumInWindow.tp73.filtered.75),5))
    cat("tail fitered@90\n")
    print(tail(sort(sum.NumInWindow.tp73.filtered.90),5))
    cat("tail fitered@95\n")
    print(tail(sort(sum.NumInWindow.tp73.filtered.95),5))
    cat("tail fitered@99\n")
    print(tail(sort(sum.NumInWindow.tp73.filtered.99),5))
    gc()

    ratio.tp73.50 <- (sum.NumInWindow.tp73.filtered.50/sum(sum.cutandrun.tp73>=sum.cutandrun.tp73.quantile.50)) / (sum.NumInWindow.tp73.equal.0/sum(sum.cutandrun.tp73==0))
    ratio.tp73.75 <- (sum.NumInWindow.tp73.filtered.75/sum(sum.cutandrun.tp73>=sum.cutandrun.tp73.quantile.75)) / (sum.NumInWindow.tp73.equal.0/sum(sum.cutandrun.tp73==0))
    ratio.tp73.90 <- (sum.NumInWindow.tp73.filtered.90/sum(sum.cutandrun.tp73>=sum.cutandrun.tp73.quantile.90)) / (sum.NumInWindow.tp73.equal.0/sum(sum.cutandrun.tp73==0))
    ratio.tp73.95 <- (sum.NumInWindow.tp73.filtered.95/sum(sum.cutandrun.tp73>=sum.cutandrun.tp73.quantile.95)) / (sum.NumInWindow.tp73.equal.0/sum(sum.cutandrun.tp73==0))

    ratio.tp73.99 <- (
                    sum.NumInWindow.tp73.filtered.99     # Number of predicted occurences of TF in TFBS that are in top 1% of signals of TP73 cut'n'run
                    /
                    sum(sum.cutandrun.tp73>=sum.cutandrun.tp73.quantile.99) # Number of TFBS that constitute the top 1%
                ) / (
                    sum.NumInWindow.tp73.equal.0         # Number of predicted occurences of TF in TFBS that have no coverage in TP73 cut'n'run
                    /
                    sum(sum.cutandrun.tp73==0)      # Number of in TFBS that have no coverage in TP73 cut'n'run
                )

    ratio.pos.50 <- (sum.NumInWindow.pos.filtered.50/sum(sum.cutandrun.pos >= sum.cutandrun.pos.quantile.50)) / (sum.NumInWindow.tp73.equal.0/sum(sum.cutandrun.pos==0))
    ratio.pos.75 <- (sum.NumInWindow.pos.filtered.75/sum(sum.cutandrun.pos >= sum.cutandrun.pos.quantile.75)) / (sum.NumInWindow.tp73.equal.0/sum(sum.cutandrun.pos==0))
    ratio.pos.90 <- (sum.NumInWindow.pos.filtered.90/sum(sum.cutandrun.pos >= sum.cutandrun.pos.quantile.90)) / (sum.NumInWindow.tp73.equal.0/sum(sum.cutandrun.pos==0))
    ratio.pos.95 <- (sum.NumInWindow.pos.filtered.95/sum(sum.cutandrun.pos >= sum.cutandrun.pos.quantile.95)) / (sum.NumInWindow.tp73.equal.0/sum(sum.cutandrun.pos==0))
    ratio.pos.99 <- (sum.NumInWindow.pos.filtered.99/sum(sum.cutandrun.pos >= sum.cutandrun.pos.quantile.99)) / (sum.NumInWindow.tp73.equal.0/sum(sum.cutandrun.pos==0))

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

    pretty.table <- function(x,sep.inner=" (") {
        x.sort <- sort(x)
        d<-data.frame(sapply(strsplit(x=names(x.sort),split="_"),function(X) paste(X[1],sep.inner,X[2],
                                                                                ifelse(sep.inner %in% c("("," ("),")",""),sep="")),
                    round(x.sort,3))
        rownames(d) <- NULL
        colnames(d) <- c("TF","Ratio")
        d
    }

    dir.create(path=reportdir,recursive=TRUE)
    write.table(file=paste(reportdir,paste("report_ratio_tp73_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.tp73.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_tp73_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.tp73.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_tp73_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.tp73.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_tp73_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.tp73.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_tp73_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.tp73.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_pos_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.pos.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_pos_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.pos.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_pos_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.pos.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_pos_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.pos.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_pos_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.pos.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    gc()

    # Enrichment 

    sum.NumInWindow.equal.TAa.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa == 0]))
    sum.NumInWindow.filtered.TAa.quantile.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa >= sum.cutandrun.tp73.TAa.quantile.50]))
    sum.NumInWindow.filtered.TAa.quantile.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa >= sum.cutandrun.tp73.TAa.quantile.75]))
    sum.NumInWindow.filtered.TAa.quantile.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa >= sum.cutandrun.tp73.TAa.quantile.90]))
    sum.NumInWindow.filtered.TAa.quantile.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa >= sum.cutandrun.tp73.TAa.quantile.95]))
    sum.NumInWindow.filtered.TAa.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.TAa >= sum.cutandrun.tp73.TAa.quantile.99]))
    ratio.TAa.quantile.50 <- (sum.NumInWindow.filtered.TAa.quantile.50/sum(sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.50)) / (sum.NumInWindow.equal.TAa.0/sum(sum.cutandrun.tp73.TAa==0))
    ratio.TAa.quantile.75 <- (sum.NumInWindow.filtered.TAa.quantile.75/sum(sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.75)) / (sum.NumInWindow.equal.TAa.0/sum(sum.cutandrun.tp73.TAa==0))
    ratio.TAa.quantile.90 <- (sum.NumInWindow.filtered.TAa.quantile.90/sum(sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.90)) / (sum.NumInWindow.equal.TAa.0/sum(sum.cutandrun.tp73.TAa==0))
    ratio.TAa.quantile.95 <- (sum.NumInWindow.filtered.TAa.quantile.95/sum(sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.95)) / (sum.NumInWindow.equal.TAa.0/sum(sum.cutandrun.tp73.TAa==0))
    ratio.TAa.quantile.99 <- (sum.NumInWindow.filtered.TAa.quantile.99/sum(sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.99)) / (sum.NumInWindow.equal.TAa.0/sum(sum.cutandrun.tp73.TAa==0))
    sum.NumInWindow.equal.DNb.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.DNb == 0]))
    sum.NumInWindow.filtered.DNb.quantile.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.DNb >= sum.cutandrun.tp73.DNb.quantile.50]))
    sum.NumInWindow.filtered.DNb.quantile.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.DNb >= sum.cutandrun.tp73.DNb.quantile.75]))
    sum.NumInWindow.filtered.DNb.quantile.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.DNb >= sum.cutandrun.tp73.DNb.quantile.90]))
    sum.NumInWindow.filtered.DNb.quantile.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.DNb >= sum.cutandrun.tp73.DNb.quantile.95]))
    sum.NumInWindow.filtered.DNb.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.DNb >= sum.cutandrun.tp73.DNb.quantile.99]))
    ratio.DNb.quantile.50 <- (sum.NumInWindow.filtered.DNb.quantile.50/sum(sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.50)) / (sum.NumInWindow.equal.DNb.0/sum(sum.cutandrun.tp73.DNb==0))
    ratio.DNb.quantile.75 <- (sum.NumInWindow.filtered.DNb.quantile.75/sum(sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.75)) / (sum.NumInWindow.equal.DNb.0/sum(sum.cutandrun.tp73.DNb==0))
    ratio.DNb.quantile.90 <- (sum.NumInWindow.filtered.DNb.quantile.99/sum(sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.90)) / (sum.NumInWindow.equal.DNb.0/sum(sum.cutandrun.tp73.DNb==0))
    ratio.DNb.quantile.95 <- (sum.NumInWindow.filtered.DNb.quantile.95/sum(sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.95)) / (sum.NumInWindow.equal.DNb.0/sum(sum.cutandrun.tp73.DNb==0))
    ratio.DNb.quantile.99 <- (sum.NumInWindow.filtered.DNb.quantile.99/sum(sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.99)) / (sum.NumInWindow.equal.DNb.0/sum(sum.cutandrun.tp73.DNb==0))
    sum.NumInWindow.equal.GFP.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.GFP == 0]))
    sum.NumInWindow.filtered.GFP.quantile.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.GFP >= sum.cutandrun.tp73.GFP.quantile.50]))
    sum.NumInWindow.filtered.GFP.quantile.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.GFP >= sum.cutandrun.tp73.GFP.quantile.75]))
    sum.NumInWindow.filtered.GFP.quantile.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.GFP >= sum.cutandrun.tp73.GFP.quantile.90]))
    sum.NumInWindow.filtered.GFP.quantile.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.GFP >= sum.cutandrun.tp73.GFP.quantile.95]))
    sum.NumInWindow.filtered.GFP.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[sum.cutandrun.tp73.GFP >= sum.cutandrun.tp73.GFP.quantile.99]))
    ratio.GFP.quantile.50 <- (sum.NumInWindow.filtered.GFP.quantile.50/sum(sum.cutandrun.tp73.GFP>=sum.cutandrun.tp73.GFP.quantile.50)) / (sum.NumInWindow.equal.GFP.0/sum(sum.cutandrun.tp73.GFP==0))
    ratio.GFP.quantile.75 <- (sum.NumInWindow.filtered.GFP.quantile.75/sum(sum.cutandrun.tp73.GFP>=sum.cutandrun.tp73.GFP.quantile.75)) / (sum.NumInWindow.equal.GFP.0/sum(sum.cutandrun.tp73.GFP==0))
    ratio.GFP.quantile.90 <- (sum.NumInWindow.filtered.GFP.quantile.90/sum(sum.cutandrun.tp73.GFP>=sum.cutandrun.tp73.GFP.quantile.90)) / (sum.NumInWindow.equal.GFP.0/sum(sum.cutandrun.tp73.GFP==0))
    ratio.GFP.quantile.95 <- (sum.NumInWindow.filtered.GFP.quantile.95/sum(sum.cutandrun.tp73.GFP>=sum.cutandrun.tp73.GFP.quantile.95)) / (sum.NumInWindow.equal.GFP.0/sum(sum.cutandrun.tp73.GFP==0))
    ratio.GFP.quantile.99 <- (sum.NumInWindow.filtered.GFP.quantile.99/sum(sum.cutandrun.tp73.GFP>=sum.cutandrun.tp73.GFP.quantile.99)) / (sum.NumInWindow.equal.GFP.0/sum(sum.cutandrun.tp73.GFP==0))

    head(pretty.table(ratio.TAa.quantile.50))
    tail(pretty.table(ratio.TAa.quantile.50))
    head(pretty.table(ratio.DNb.quantile.50))
    tail(pretty.table(ratio.DNb.quantile.50))
    head(pretty.table(ratio.GFP.quantile.50))
    tail(pretty.table(ratio.GFP.quantile.50))
    head(pretty.table(ratio.TAa.quantile.99))
    tail(pretty.table(ratio.TAa.quantile.99))
    head(pretty.table(ratio.DNb.quantile.99))
    tail(pretty.table(ratio.DNb.quantile.99))
    head(pretty.table(ratio.GFP.quantile.99))
    tail(pretty.table(ratio.GFP.quantile.99))

    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_empty_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.TAa.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_empty_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.TAa.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_empty_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.TAa.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_empty_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.TAa.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_empty_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.TAa.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_empty_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.DNb.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_empty_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.DNb.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_empty_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.DNb.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_empty_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.DNb.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_empty_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.DNb.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_empty_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.GFP.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_empty_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.GFP.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_empty_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.GFP.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_empty_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.GFP.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_empty_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.GFP.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)


    ratio.TAa.quantile.50.vs.ratio.DNb.quantile.50 <- ratio.TAa.quantile.50 / ratio.DNb.quantile.50
    head(pretty.table(ratio.TAa.quantile.50.vs.ratio.DNb.quantile.50))
    tail(pretty.table(ratio.TAa.quantile.50.vs.ratio.DNb.quantile.50))
    ratio.TAa.quantile.75.vs.ratio.DNb.quantile.75 <- ratio.TAa.quantile.75 / ratio.DNb.quantile.75
    ratio.TAa.quantile.90.vs.ratio.DNb.quantile.90 <- ratio.TAa.quantile.90 / ratio.DNb.quantile.90
    ratio.TAa.quantile.95.vs.ratio.DNb.quantile.95 <- ratio.TAa.quantile.95 / ratio.DNb.quantile.95
    ratio.TAa.quantile.99.vs.ratio.DNb.quantile.99 <- ratio.TAa.quantile.99 / ratio.DNb.quantile.99
    head(pretty.table(ratio.TAa.quantile.99.vs.ratio.DNb.quantile.99))
    # Auffälligkeit p53
    tail(pretty.table(ratio.TAa.quantile.99.vs.ratio.DNb.quantile.99))


    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs__DNb_vs_empty_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.TAa.quantile.50.vs.ratio.DNb.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs__DNb_vs_empty_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.TAa.quantile.75.vs.ratio.DNb.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs__DNb_vs_empty_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.TAa.quantile.90.vs.ratio.DNb.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs__DNb_vs_empty_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.TAa.quantile.95.vs.ratio.DNb.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs__DNb_vs_empty_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.TAa.quantile.99.vs.ratio.DNb.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)



    sum.NumInWindow.equal.saos2.TAa.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_TA" == 0]))
    sum.NumInWindow.filtered.saos2.TAa.quantile.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_TA" >= sum.cutandrun.tp73.saos2.TAa.quantile.50]))
    sum.NumInWindow.filtered.saos2.TAa.quantile.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_TA" >= sum.cutandrun.tp73.saos2.TAa.quantile.75]))
    sum.NumInWindow.filtered.saos2.TAa.quantile.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_TA" >= sum.cutandrun.tp73.saos2.TAa.quantile.90]))
    sum.NumInWindow.filtered.saos2.TAa.quantile.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_TA" >= sum.cutandrun.tp73.saos2.TAa.quantile.95]))
    sum.NumInWindow.filtered.saos2.TAa.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_TA" >= sum.cutandrun.tp73.saos2.TAa.quantile.99]))
                    # Number of TFBS for particular TF in top 1% of CUT&RUN confirmed sites
                                                        # Number of TFBS in top 1% of confirmed sites
    ratio.saos2.TAa.quantile.50 <- (sum.NumInWindow.filtered.saos2.TAa.quantile.50/sum(m$"tp73_saos2_TA">=sum.cutandrun.tp73.saos2.TAa.quantile.50)) / (sum.NumInWindow.equal.saos2.TAa.0/sum(m$"tp73_saos2_TA"==0))
    ratio.saos2.TAa.quantile.75 <- (sum.NumInWindow.filtered.saos2.TAa.quantile.75/sum(m$"tp73_saos2_TA">=sum.cutandrun.tp73.saos2.TAa.quantile.75)) / (sum.NumInWindow.equal.saos2.TAa.0/sum(m$"tp73_saos2_TA"==0))
    ratio.saos2.TAa.quantile.90 <- (sum.NumInWindow.filtered.saos2.TAa.quantile.90/sum(m$"tp73_saos2_TA">=sum.cutandrun.tp73.saos2.TAa.quantile.90)) / (sum.NumInWindow.equal.saos2.TAa.0/sum(m$"tp73_saos2_TA"==0))
    ratio.saos2.TAa.quantile.95 <- (sum.NumInWindow.filtered.saos2.TAa.quantile.95/sum(m$"tp73_saos2_TA">=sum.cutandrun.tp73.saos2.TAa.quantile.95)) / (sum.NumInWindow.equal.saos2.TAa.0/sum(m$"tp73_saos2_TA"==0))
    ratio.saos2.TAa.quantile.99 <- (sum.NumInWindow.filtered.saos2.TAa.quantile.99/sum(m$"tp73_saos2_TA">=sum.cutandrun.tp73.saos2.TAa.quantile.99)) / (sum.NumInWindow.equal.saos2.TAa.0/sum(m$"tp73_saos2_TA"==0))

    sum.NumInWindow.equal.saos2.DNb.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_DN" == 0]))
    sum.NumInWindow.filtered.saos2.DNb.quantile.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_DN" >= sum.cutandrun.tp73.saos2.DNb.quantile.50]))
    sum.NumInWindow.filtered.saos2.DNb.quantile.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_DN" >= sum.cutandrun.tp73.saos2.DNb.quantile.75]))
    sum.NumInWindow.filtered.saos2.DNb.quantile.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_DN" >= sum.cutandrun.tp73.saos2.DNb.quantile.90]))
    sum.NumInWindow.filtered.saos2.DNb.quantile.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_DN" >= sum.cutandrun.tp73.saos2.DNb.quantile.95]))
    sum.NumInWindow.filtered.saos2.DNb.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_DN" >= sum.cutandrun.tp73.saos2.DNb.quantile.99]))
    ratio.saos2.DNb.quantile.50 <- (sum.NumInWindow.filtered.saos2.DNb.quantile.50/sum(m$"tp73_saos2_DN">=sum.cutandrun.tp73.saos2.DNb.quantile.50)) / (sum.NumInWindow.equal.saos2.DNb.0/sum(m$"tp73_saos2_DN"==0))
    ratio.saos2.DNb.quantile.75 <- (sum.NumInWindow.filtered.saos2.DNb.quantile.75/sum(m$"tp73_saos2_DN">=sum.cutandrun.tp73.saos2.DNb.quantile.75)) / (sum.NumInWindow.equal.saos2.DNb.0/sum(m$"tp73_saos2_DN"==0))
    ratio.saos2.DNb.quantile.90 <- (sum.NumInWindow.filtered.saos2.DNb.quantile.90/sum(m$"tp73_saos2_DN">=sum.cutandrun.tp73.saos2.DNb.quantile.90)) / (sum.NumInWindow.equal.saos2.DNb.0/sum(m$"tp73_saos2_DN"==0))
    ratio.saos2.DNb.quantile.95 <- (sum.NumInWindow.filtered.saos2.DNb.quantile.95/sum(m$"tp73_saos2_DN">=sum.cutandrun.tp73.saos2.DNb.quantile.95)) / (sum.NumInWindow.equal.saos2.DNb.0/sum(m$"tp73_saos2_DN"==0))
    ratio.saos2.DNb.quantile.99 <- (sum.NumInWindow.filtered.saos2.DNb.quantile.99/sum(m$"tp73_saos2_DN">=sum.cutandrun.tp73.saos2.DNb.quantile.99)) / (sum.NumInWindow.equal.saos2.DNb.0/sum(m$"tp73_saos2_DN"==0))

    sum.NumInWindow.equal.saos2.GFP.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_GFP" == 0]))
    sum.NumInWindow.filtered.saos2.GFP.quantile.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_GFP" >= sum.cutandrun.tp73.saos2.GFP.quantile.50]))
    sum.NumInWindow.filtered.saos2.GFP.quantile.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_GFP" >= sum.cutandrun.tp73.saos2.GFP.quantile.75]))
    sum.NumInWindow.filtered.saos2.GFP.quantile.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_GFP" >= sum.cutandrun.tp73.saos2.GFP.quantile.90]))
    sum.NumInWindow.filtered.saos2.GFP.quantile.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_GFP" >= sum.cutandrun.tp73.saos2.GFP.quantile.95]))
    sum.NumInWindow.filtered.saos2.GFP.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_saos2_GFP" >= sum.cutandrun.tp73.saos2.GFP.quantile.99]))
    ratio.saos2.GFP.quantile.50 <- (sum.NumInWindow.filtered.saos2.GFP.quantile.50/sum(m$"tp73_saos2_GFP">=sum.cutandrun.tp73.saos2.GFP.quantile.50)) / (sum.NumInWindow.equal.saos2.GFP.0/sum(m$"tp73_saos2_GFP"==0))
    ratio.saos2.GFP.quantile.75 <- (sum.NumInWindow.filtered.saos2.GFP.quantile.75/sum(m$"tp73_saos2_GFP">=sum.cutandrun.tp73.saos2.GFP.quantile.75)) / (sum.NumInWindow.equal.saos2.GFP.0/sum(m$"tp73_saos2_GFP"==0))
    ratio.saos2.GFP.quantile.90 <- (sum.NumInWindow.filtered.saos2.GFP.quantile.90/sum(m$"tp73_saos2_GFP">=sum.cutandrun.tp73.saos2.GFP.quantile.90)) / (sum.NumInWindow.equal.saos2.GFP.0/sum(m$"tp73_saos2_GFP"==0))
    ratio.saos2.GFP.quantile.95 <- (sum.NumInWindow.filtered.saos2.GFP.quantile.95/sum(m$"tp73_saos2_GFP">=sum.cutandrun.tp73.saos2.GFP.quantile.95)) / (sum.NumInWindow.equal.saos2.GFP.0/sum(m$"tp73_saos2_GFP"==0))
    ratio.saos2.GFP.quantile.99 <- (sum.NumInWindow.filtered.saos2.GFP.quantile.99/sum(m$"tp73_saos2_GFP">=sum.cutandrun.tp73.saos2.GFP.quantile.99)) / (sum.NumInWindow.equal.saos2.GFP.0/sum(m$"tp73_saos2_GFP"==0))

    head(pretty.table(ratio.saos2.TAa.quantile.99))
    tail(pretty.table(ratio.saos2.TAa.quantile.99))
    head(pretty.table(ratio.saos2.DNb.quantile.99))
    tail(pretty.table(ratio.saos2.DNb.quantile.99))
    head(pretty.table(ratio.saos2.GFP.quantile.99))
    tail(pretty.table(ratio.saos2.GFP.quantile.99))

    ratio.saos2.TAa.99.vs.ratio.saos2.DNb.99 <- ratio.saos2.TAa.quantile.99 / ratio.saos2.DNb.quantile.99
    pretty.table(ratio.saos2.TAa.99.vs.ratio.saos2.DNb.99)
    ratio.saos2.TAa.99.vs.ratio.saos2.GFP.99 <- ratio.saos2.TAa.quantile.99 / ratio.saos2.GFP.quantile.99
    pretty.table(ratio.saos2.TAa.99.vs.ratio.saos2.GFP.99)
    ratio.saos2.DNb.99.vs.ratio.saos2.GFP.99 <- ratio.saos2.DNb.quantile.99 / ratio.saos2.GFP.quantile.99
    pretty.table(ratio.saos2.DNb.99.vs.ratio.saos2.GFP.99)

    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_noTAa_in_Saos2_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.TAa.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_noTAa_in_Saos2_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.TAa.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_noTAa_in_Saos2_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.TAa.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_noTAa_in_Saos2_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.TAa.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_noTAa_in_Saos2_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.TAa.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_noDNb_in_Saos2_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.DNb.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_noDNb_in_Saos2_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.DNb.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_noDNb_in_Saos2_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.DNb.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_noDNb_in_Saos2_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.DNb.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_noDNb_in_Saos2_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.DNb.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_noGFP_in_Saos2_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.GFP.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_noGFP_in_Saos2_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.GFP.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_noGFP_in_Saos2_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.GFP.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_noGFP_in_Saos2_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.GFP.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_noGFP_in_Saos2_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.GFP.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)


    sum.NumInWindow.equal.skmel29_2.TAa.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_TA" == 0]))
    sum.NumInWindow.filtered.skmel29_2.TAa.quantile.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_TA" >= sum.cutandrun.tp73.skmel29_2.TAa.quantile.50]))
    sum.NumInWindow.filtered.skmel29_2.TAa.quantile.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_TA" >= sum.cutandrun.tp73.skmel29_2.TAa.quantile.75]))
    sum.NumInWindow.filtered.skmel29_2.TAa.quantile.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_TA" >= sum.cutandrun.tp73.skmel29_2.TAa.quantile.90]))
    sum.NumInWindow.filtered.skmel29_2.TAa.quantile.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_TA" >= sum.cutandrun.tp73.skmel29_2.TAa.quantile.95]))
    sum.NumInWindow.filtered.skmel29_2.TAa.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_TA" >= sum.cutandrun.tp73.skmel29_2.TAa.quantile.99]))
                    # Number of TFBS for particular TF in top 1% of CUT&RUN confirmed sites
                                                        # Number of TFBS in top 1% of confirmed sites
    ratio.skmel29_2.TAa.quantile.50 <- (sum.NumInWindow.filtered.skmel29_2.TAa.quantile.50/sum(m$"tp73_skmel29_2_TA">=sum.cutandrun.tp73.skmel29_2.TAa.quantile.50)) / (sum.NumInWindow.equal.skmel29_2.TAa.0/sum(m$"tp73_skmel29_2_TA"==0))
    ratio.skmel29_2.TAa.quantile.75 <- (sum.NumInWindow.filtered.skmel29_2.TAa.quantile.75/sum(m$"tp73_skmel29_2_TA">=sum.cutandrun.tp73.skmel29_2.TAa.quantile.75)) / (sum.NumInWindow.equal.skmel29_2.TAa.0/sum(m$"tp73_skmel29_2_TA"==0))
    ratio.skmel29_2.TAa.quantile.90 <- (sum.NumInWindow.filtered.skmel29_2.TAa.quantile.90/sum(m$"tp73_skmel29_2_TA">=sum.cutandrun.tp73.skmel29_2.TAa.quantile.90)) / (sum.NumInWindow.equal.skmel29_2.TAa.0/sum(m$"tp73_skmel29_2_TA"==0))
    ratio.skmel29_2.TAa.quantile.95 <- (sum.NumInWindow.filtered.skmel29_2.TAa.quantile.95/sum(m$"tp73_skmel29_2_TA">=sum.cutandrun.tp73.skmel29_2.TAa.quantile.95)) / (sum.NumInWindow.equal.skmel29_2.TAa.0/sum(m$"tp73_skmel29_2_TA"==0))
    ratio.skmel29_2.TAa.quantile.99 <- (sum.NumInWindow.filtered.skmel29_2.TAa.quantile.99/sum(m$"tp73_skmel29_2_TA">=sum.cutandrun.tp73.skmel29_2.TAa.quantile.99)) / (sum.NumInWindow.equal.skmel29_2.TAa.0/sum(m$"tp73_skmel29_2_TA"==0))

    sum.NumInWindow.equal.skmel29_2.DNb.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_DN" == 0]))
    sum.NumInWindow.filtered.skmel29_2.DNb.quantile.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_DN" >= sum.cutandrun.tp73.skmel29_2.DNb.quantile.50]))
    sum.NumInWindow.filtered.skmel29_2.DNb.quantile.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_DN" >= sum.cutandrun.tp73.skmel29_2.DNb.quantile.75]))
    sum.NumInWindow.filtered.skmel29_2.DNb.quantile.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_DN" >= sum.cutandrun.tp73.skmel29_2.DNb.quantile.90]))
    sum.NumInWindow.filtered.skmel29_2.DNb.quantile.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_DN" >= sum.cutandrun.tp73.skmel29_2.DNb.quantile.95]))
    sum.NumInWindow.filtered.skmel29_2.DNb.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_DN" >= sum.cutandrun.tp73.skmel29_2.DNb.quantile.99]))
    ratio.skmel29_2.DNb.quantile.50 <- (sum.NumInWindow.filtered.skmel29_2.DNb.quantile.50/sum(m$"tp73_skmel29_2_DN">=sum.cutandrun.tp73.skmel29_2.DNb.quantile.50)) / (sum.NumInWindow.equal.skmel29_2.DNb.0/sum(m$"tp73_skmel29_2_DN"==0))
    ratio.skmel29_2.DNb.quantile.75 <- (sum.NumInWindow.filtered.skmel29_2.DNb.quantile.75/sum(m$"tp73_skmel29_2_DN">=sum.cutandrun.tp73.skmel29_2.DNb.quantile.75)) / (sum.NumInWindow.equal.skmel29_2.DNb.0/sum(m$"tp73_skmel29_2_DN"==0))
    ratio.skmel29_2.DNb.quantile.90 <- (sum.NumInWindow.filtered.skmel29_2.DNb.quantile.90/sum(m$"tp73_skmel29_2_DN">=sum.cutandrun.tp73.skmel29_2.DNb.quantile.90)) / (sum.NumInWindow.equal.skmel29_2.DNb.0/sum(m$"tp73_skmel29_2_DN"==0))
    ratio.skmel29_2.DNb.quantile.95 <- (sum.NumInWindow.filtered.skmel29_2.DNb.quantile.95/sum(m$"tp73_skmel29_2_DN">=sum.cutandrun.tp73.skmel29_2.DNb.quantile.95)) / (sum.NumInWindow.equal.skmel29_2.DNb.0/sum(m$"tp73_skmel29_2_DN"==0))
    ratio.skmel29_2.DNb.quantile.99 <- (sum.NumInWindow.filtered.skmel29_2.DNb.quantile.99/sum(m$"tp73_skmel29_2_DN">=sum.cutandrun.tp73.skmel29_2.DNb.quantile.99)) / (sum.NumInWindow.equal.skmel29_2.DNb.0/sum(m$"tp73_skmel29_2_DN"==0))

    sum.NumInWindow.equal.skmel29_2.GFP.0 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_GFP" == 0]))
    sum.NumInWindow.filtered.skmel29_2.GFP.quantile.50 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_GFP" >= sum.cutandrun.tp73.skmel29_2.GFP.quantile.50]))
    sum.NumInWindow.filtered.skmel29_2.GFP.quantile.75 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_GFP" >= sum.cutandrun.tp73.skmel29_2.GFP.quantile.75]))
    sum.NumInWindow.filtered.skmel29_2.GFP.quantile.90 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_GFP" >= sum.cutandrun.tp73.skmel29_2.GFP.quantile.90]))
    sum.NumInWindow.filtered.skmel29_2.GFP.quantile.95 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_GFP" >= sum.cutandrun.tp73.skmel29_2.GFP.quantile.95]))
    sum.NumInWindow.filtered.skmel29_2.GFP.quantile.99 <- sapply(m[, ..cols.NumInWindow], function(X) sum(X[m$"tp73_skmel29_2_GFP" >= sum.cutandrun.tp73.skmel29_2.GFP.quantile.99]))
    ratio.skmel29_2.GFP.quantile.50 <- (sum.NumInWindow.filtered.skmel29_2.GFP.quantile.50/sum(m$"tp73_skmel29_2_GFP">=sum.cutandrun.tp73.skmel29_2.GFP.quantile.50)) / (sum.NumInWindow.equal.skmel29_2.GFP.0/sum(m$"tp73_skmel29_2_GFP"==0))
    ratio.skmel29_2.GFP.quantile.75 <- (sum.NumInWindow.filtered.skmel29_2.GFP.quantile.75/sum(m$"tp73_skmel29_2_GFP">=sum.cutandrun.tp73.skmel29_2.GFP.quantile.75)) / (sum.NumInWindow.equal.skmel29_2.GFP.0/sum(m$"tp73_skmel29_2_GFP"==0))
    ratio.skmel29_2.GFP.quantile.90 <- (sum.NumInWindow.filtered.skmel29_2.GFP.quantile.90/sum(m$"tp73_skmel29_2_GFP">=sum.cutandrun.tp73.skmel29_2.GFP.quantile.90)) / (sum.NumInWindow.equal.skmel29_2.GFP.0/sum(m$"tp73_skmel29_2_GFP"==0))
    ratio.skmel29_2.GFP.quantile.95 <- (sum.NumInWindow.filtered.skmel29_2.GFP.quantile.95/sum(m$"tp73_skmel29_2_GFP">=sum.cutandrun.tp73.skmel29_2.GFP.quantile.95)) / (sum.NumInWindow.equal.skmel29_2.GFP.0/sum(m$"tp73_skmel29_2_GFP"==0))
    ratio.skmel29_2.GFP.quantile.99 <- (sum.NumInWindow.filtered.skmel29_2.GFP.quantile.99/sum(m$"tp73_skmel29_2_GFP">=sum.cutandrun.tp73.skmel29_2.GFP.quantile.99)) / (sum.NumInWindow.equal.skmel29_2.GFP.0/sum(m$"tp73_skmel29_2_GFP"==0))

    head(pretty.table(ratio.skmel29_2.TAa.quantile.99),10)
    tail(pretty.table(ratio.skmel29_2.TAa.quantile.99),10)
    head(pretty.table(ratio.skmel29_2.DNb.quantile.99),10)
    tail(pretty.table(ratio.skmel29_2.DNb.quantile.99),10)
    head(pretty.table(ratio.skmel29_2.GFP.quantile.99),10)
    tail(pretty.table(ratio.skmel29_2.GFP.quantile.99),10)

    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_noTAa_in_SKMel29_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_noTAa_in_SKMel29_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_noTAa_in_SKMel29_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_noTAa_in_SKMel29_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_TAa_vs_noTAa_in_SKMel29_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_noDNb_in_SKMel29_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_noDNb_in_SKMel29_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_noDNb_in_SKMel29_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_noDNb_in_SKMel29_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_DNb_vs_noDNb_in_SKMel29_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_noGFP_in_SKMel29_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.GFP.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_noGFP_in_SKMel29_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.GFP.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_noGFP_in_SKMel29_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.GFP.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_noGFP_in_SKMel29_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.GFP.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_GFP_vs_noGFP_in_SKMel29_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.GFP.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    #
    # tissue specific ratios of TFBS
    #

    ratio.skmel29_2.TAa.quantile.50.vs.ratio.skmel29_2.DNb.quantile.50 <- ratio.skmel29_2.TAa.quantile.50 / ratio.skmel29_2.DNb.quantile.50
    ratio.skmel29_2.TAa.quantile.75.vs.ratio.skmel29_2.DNb.quantile.75 <- ratio.skmel29_2.TAa.quantile.75 / ratio.skmel29_2.DNb.quantile.75
    ratio.skmel29_2.TAa.quantile.90.vs.ratio.skmel29_2.DNb.quantile.90 <- ratio.skmel29_2.TAa.quantile.90 / ratio.skmel29_2.DNb.quantile.90
    ratio.skmel29_2.TAa.quantile.95.vs.ratio.skmel29_2.DNb.quantile.95 <- ratio.skmel29_2.TAa.quantile.95 / ratio.skmel29_2.DNb.quantile.95
    ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.DNb.quantile.99 <- ratio.skmel29_2.TAa.quantile.99 / ratio.skmel29_2.DNb.quantile.99

    head(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.DNb.quantile.99),10)
    tail(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.DNb.quantile.99),10)

    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs__DNb_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.50.vs.ratio.skmel29_2.DNb.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs__DNb_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.75.vs.ratio.skmel29_2.DNb.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs__DNb_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.90.vs.ratio.skmel29_2.DNb.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs__DNb_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.95.vs.ratio.skmel29_2.DNb.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs__DNb_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.DNb.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    ratio.skmel29_2.TAa.quantile.50.vs.ratio.skmel29_2.GFP.quantile.50 <- ratio.skmel29_2.TAa.quantile.50 / ratio.skmel29_2.GFP.quantile.50
    ratio.skmel29_2.TAa.quantile.75.vs.ratio.skmel29_2.GFP.quantile.75 <- ratio.skmel29_2.TAa.quantile.75 / ratio.skmel29_2.GFP.quantile.75
    ratio.skmel29_2.TAa.quantile.90.vs.ratio.skmel29_2.GFP.quantile.90 <- ratio.skmel29_2.TAa.quantile.90 / ratio.skmel29_2.GFP.quantile.90
    ratio.skmel29_2.TAa.quantile.95.vs.ratio.skmel29_2.GFP.quantile.95 <- ratio.skmel29_2.TAa.quantile.95 / ratio.skmel29_2.GFP.quantile.95
    ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99 <- ratio.skmel29_2.TAa.quantile.99 / ratio.skmel29_2.GFP.quantile.99
    head(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99),10)
    tail(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99),10)

    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs_GFP_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.50.vs.ratio.skmel29_2.GFP.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs_GFP_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.75.vs.ratio.skmel29_2.GFP.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs_GFP_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.90.vs.ratio.skmel29_2.GFP.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs_GFP_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.95.vs.ratio.skmel29_2.GFP.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_vs_empty__vs_GFP_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    ratio.skmel29_2.DNb.quantile.50.vs.ratio.skmel29_2.GFP.quantile.50 <- ratio.skmel29_2.DNb.quantile.50 / ratio.skmel29_2.GFP.quantile.50
    ratio.skmel29_2.DNb.quantile.75.vs.ratio.skmel29_2.GFP.quantile.75 <- ratio.skmel29_2.DNb.quantile.75 / ratio.skmel29_2.GFP.quantile.75
    ratio.skmel29_2.DNb.quantile.90.vs.ratio.skmel29_2.GFP.quantile.90 <- ratio.skmel29_2.DNb.quantile.90 / ratio.skmel29_2.GFP.quantile.90
    ratio.skmel29_2.DNb.quantile.95.vs.ratio.skmel29_2.GFP.quantile.95 <- ratio.skmel29_2.DNb.quantile.95 / ratio.skmel29_2.GFP.quantile.95
    ratio.skmel29_2.DNb.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99 <- ratio.skmel29_2.DNb.quantile.99 / ratio.skmel29_2.GFP.quantile.99

    head(pretty.table(ratio.skmel29_2.DNb.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99),10)
    tail(pretty.table(ratio.skmel29_2.DNb.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99),10)

    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_DNb_vs_empty__vs_GFP_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.50.vs.ratio.skmel29_2.GFP.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_DNb_vs_empty__vs_GFP_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.75.vs.ratio.skmel29_2.GFP.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_DNb_vs_empty__vs_GFP_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.90.vs.ratio.skmel29_2.GFP.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_DNb_vs_empty__vs_GFP_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.95.vs.ratio.skmel29_2.GFP.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_DNb_vs_empty__vs_GFP_vs_empty_in_SKMel29_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.99.vs.ratio.skmel29_2.GFP.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)


    #
    # cross-tissue comparisons
    #

    ratio.skmel29_2.TAa.quantile.50.vs.ratio.saos2.TAa.quantile.50 <- ratio.skmel29_2.TAa.quantile.50 / ratio.saos2.TAa.quantile.50
    ratio.skmel29_2.TAa.quantile.75.vs.ratio.saos2.TAa.quantile.75 <- ratio.skmel29_2.TAa.quantile.99 / ratio.saos2.TAa.quantile.75
    ratio.skmel29_2.TAa.quantile.90.vs.ratio.saos2.TAa.quantile.90 <- ratio.skmel29_2.TAa.quantile.90 / ratio.saos2.TAa.quantile.90
    ratio.skmel29_2.TAa.quantile.95.vs.ratio.saos2.TAa.quantile.95 <- ratio.skmel29_2.TAa.quantile.95 / ratio.saos2.TAa.quantile.95
    ratio.skmel29_2.TAa.quantile.99.vs.ratio.saos2.TAa.quantile.99 <- ratio.skmel29_2.TAa.quantile.99 / ratio.saos2.TAa.quantile.99

    head(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.saos2.TAa.quantile.99),10)
    tail(pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.saos2.TAa.quantile.99),10)

    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_in_SKMel29_vs_empty__vs__TAa_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.50.vs.ratio.saos2.TAa.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_in_SKMel29_vs_empty__vs__TAa_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.75.vs.ratio.saos2.TAa.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_in_SKMel29_vs_empty__vs__TAa_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.90.vs.ratio.saos2.TAa.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_in_SKMel29_vs_empty__vs__TAa_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.95.vs.ratio.saos2.TAa.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_TAa_in_SKMel29_vs_empty__vs__TAa_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.99.vs.ratio.saos2.TAa.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    ratio.skmel29_2.DNb.quantile.50.vs.ratio.saos2.DNb.quantile.50 <- ratio.skmel29_2.DNb.quantile.50 / ratio.saos2.DNb.quantile.50
    ratio.skmel29_2.DNb.quantile.75.vs.ratio.saos2.DNb.quantile.75 <- ratio.skmel29_2.DNb.quantile.75 / ratio.saos2.DNb.quantile.75
    ratio.skmel29_2.DNb.quantile.90.vs.ratio.saos2.DNb.quantile.90 <- ratio.skmel29_2.DNb.quantile.90 / ratio.saos2.DNb.quantile.90
    ratio.skmel29_2.DNb.quantile.95.vs.ratio.saos2.DNb.quantile.95 <- ratio.skmel29_2.DNb.quantile.95 / ratio.saos2.DNb.quantile.95
    ratio.skmel29_2.DNb.quantile.99.vs.ratio.saos2.DNb.quantile.99 <- ratio.skmel29_2.DNb.quantile.99 / ratio.saos2.DNb.quantile.99

    head(pretty.table(ratio.skmel29_2.DNb.quantile.99.vs.ratio.saos2.DNb.quantile.99),10)
    tail(pretty.table(ratio.skmel29_2.DNb.quantile.99.vs.ratio.saos2.DNb.quantile.99),10)

    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_DNb_in_SKMel29_vs_empty__vs__DNb_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.50.vs.ratio.saos2.DNb.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_DNb_in_SKMel29_vs_empty__vs__DNb_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.75.vs.ratio.saos2.DNb.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_DNb_in_SKMel29_vs_empty__vs__DNb_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.90.vs.ratio.saos2.DNb.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_DNb_in_SKMel29_vs_empty__vs__DNb_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.95.vs.ratio.saos2.DNb.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_DNb_in_SKMel29_vs_empty__vs__DNb_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.DNb.quantile.99.vs.ratio.saos2.DNb.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    ratio.skmel29_2.GFP.quantile.50.vs.ratio.saos2.GFP.quantile.50 <- ratio.skmel29_2.GFP.quantile.50 / ratio.saos2.GFP.quantile.50
    ratio.skmel29_2.GFP.quantile.75.vs.ratio.saos2.GFP.quantile.75 <- ratio.skmel29_2.GFP.quantile.75 / ratio.saos2.GFP.quantile.75
    ratio.skmel29_2.GFP.quantile.90.vs.ratio.saos2.GFP.quantile.90 <- ratio.skmel29_2.GFP.quantile.90 / ratio.saos2.GFP.quantile.90
    ratio.skmel29_2.GFP.quantile.95.vs.ratio.saos2.GFP.quantile.95 <- ratio.skmel29_2.GFP.quantile.95 / ratio.saos2.GFP.quantile.95
    ratio.skmel29_2.GFP.quantile.99.vs.ratio.saos2.GFP.quantile.99 <- ratio.skmel29_2.GFP.quantile.99 / ratio.saos2.GFP.quantile.99

    head(pretty.table(ratio.skmel29_2.GFP.quantile.99.vs.ratio.saos2.GFP.quantile.99),10)
    tail(pretty.table(ratio.skmel29_2.GFP.quantile.99.vs.ratio.saos2.GFP.quantile.99),10)

    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_GFP_in_SKMel29_vs_empty__vs__GFP_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",50,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.GFP.quantile.50.vs.ratio.saos2.GFP.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_GFP_in_SKMel29_vs_empty__vs__GFP_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",75,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.GFP.quantile.75.vs.ratio.saos2.GFP.quantile.75),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_GFP_in_SKMel29_vs_empty__vs__GFP_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",90,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.GFP.quantile.90.vs.ratio.saos2.GFP.quantile.90),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_GFP_in_SKMel29_vs_empty__vs__GFP_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",95,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.GFP.quantile.95.vs.ratio.saos2.GFP.quantile.95),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)
    write.table(file=paste(reportdir,paste("report_ratio_of_ratios_GFP_in_SKMel29_vs_empty__vs__GFP_in_Saos2_vs_empty_chr_",chromosome,"_quantile_",99,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.GFP.quantile.99.vs.ratio.saos2.GFP.quantile.99),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)


    # Singling out TAa and DNb

    mean.NumInWindow.filtered.saos2.TAa.quantile.90.wo.DNb.quantile.50 <- 
        sapply(m[, ..cols.NumInWindow], function(X) mean(X[
            m$"tp73_saos2_TA" >= sum.cutandrun.tp73.saos2.TAa.quantile.90 & m$"tp73_saos2_DN" <= sum.cutandrun.tp73.saos2.DNb.quantile.50
            ]))

    mean.NumInWindow.filtered.saos2.DNb.quantile.90.wo.TAa.quantile.50 <- 
        sapply(m[, ..cols.NumInWindow], function(X) mean(X[
            m$"tp73_saos2_TA" >= sum.cutandrun.tp73.saos2.DNb.quantile.90 & m$"tp73_saos2_DN" <= sum.cutandrun.tp73.saos2.TAa.quantile.50
            ]))

    ratio.saos2.TAa.quantile.90.wo.DNb.quantile.50.vs.saos2.DNb.quantile.90.wo.TAa.quantile.50 <-
        mean.NumInWindow.filtered.saos2.TAa.quantile.90.wo.DNb.quantile.50 / mean.NumInWindow.filtered.saos2.DNb.quantile.90.wo.TAa.quantile.50

    head(pretty.table(ratio.saos2.TAa.quantile.90.wo.DNb.quantile.50.vs.saos2.DNb.quantile.90.wo.TAa.quantile.50),10)
    tail(pretty.table(ratio.saos2.TAa.quantile.90.wo.DNb.quantile.50.vs.saos2.DNb.quantile.90.wo.TAa.quantile.50),10)

    write.table(file=paste(reportdir,paste("report_ratio_TAa_quantile_90_without_DNb_quantile_50__vs__DNb_quantile_90_without_TAa_quantile_50_in_Saos2_chr_",chromosome,".tsv",sep=""),sep="/"),x=pretty.table(ratio.saos2.TAa.quantile.90.wo.DNb.quantile.50.vs.saos2.DNb.quantile.90.wo.TAa.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    mean.NumInWindow.filtered.skmel29_2.TAa.quantile.90.wo.DNb.quantile.50 <- 
        sapply(m[, ..cols.NumInWindow], function(X) mean(X[
            m$"tp73_skmel29_2_TA" >= sum.cutandrun.tp73.skmel29_2.TAa.quantile.90 & m$"tp73_skmel29_2_DN" <= sum.cutandrun.tp73.skmel29_2.DNb.quantile.50
            ]))

    mean.NumInWindow.filtered.skmel29_2.DNb.quantile.90.wo.TAa.quantile.50 <- 
        sapply(m[, ..cols.NumInWindow], function(X) mean(X[
            m$"tp73_skmel29_2_TA" >= sum.cutandrun.tp73.skmel29_2.DNb.quantile.90 & m$"tp73_skmel29_2_DN" <= sum.cutandrun.tp73.skmel29_2.TAa.quantile.50
            ]))

    ratio.skmel29_2.TAa.quantile.90.wo.DNb.quantile.50.vs.skmel29_2.DNb.quantile.90.wo.TAa.quantile.50 <- 
        mean.NumInWindow.filtered.skmel29_2.TAa.quantile.90.wo.DNb.quantile.50 / mean.NumInWindow.filtered.skmel29_2.DNb.quantile.90.wo.TAa.quantile.50

    head(pretty.table(ratio.skmel29_2.TAa.quantile.90.wo.DNb.quantile.50.vs.skmel29_2.DNb.quantile.90.wo.TAa.quantile.50),10)
    tail(pretty.table(ratio.skmel29_2.TAa.quantile.90.wo.DNb.quantile.50.vs.skmel29_2.DNb.quantile.90.wo.TAa.quantile.50),10)

    write.table(file=paste(reportdir,paste("report_ratio_TAa_quantile_90_without_DNb_quantile_50__vs__DNb_quantile_90_without_TAa_quantile_50_in_SKMel29_chr_",chromosome,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.90.wo.DNb.quantile.50.vs.skmel29_2.DNb.quantile.90.wo.TAa.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    ratio.ratio.skmel29_2.vs.saos2.of.TAa.wo.DNb.vs.DNb.wo.TAa <- 
        ratio.skmel29_2.TAa.quantile.90.wo.DNb.quantile.50.vs.skmel29_2.DNb.quantile.90.wo.TAa.quantile.50 /
        ratio.saos2.TAa.quantile.90.wo.DNb.quantile.50.vs.saos2.DNb.quantile.90.wo.TAa.quantile.50

    head(pretty.table(ratio.ratio.skmel29_2.vs.saos2.of.TAa.wo.DNb.vs.DNb.wo.TAa),20)
    tail(pretty.table(ratio.ratio.skmel29_2.vs.saos2.of.TAa.wo.DNb.vs.DNb.wo.TAa),20)    

    write.table(file=paste(reportdir,paste("report_ratio_SKMel_vs_Saos2_of_TAa_without_DNb_vs_DNb_without_TAa_chr_",chromosome,".tsv",sep=""),sep="/"),x=pretty.table(ratio.skmel29_2.TAa.quantile.90.wo.DNb.quantile.50.vs.skmel29_2.DNb.quantile.90.wo.TAa.quantile.50),sep="\t",col.names=TRUE,row.names=FALSE,append=FALSE,quote=FALSE)

    invisible(m)

}

