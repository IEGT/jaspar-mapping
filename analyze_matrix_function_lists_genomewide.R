
m.contexts.all.OneTo18 <- do.call(rbind, lapply(m.contexts, function(x) x[, 1:18, drop = FALSE]))
# Calculate the fraction of entries with value > 0 for "tp73_skmel29_2_DN" and "tp73_skmel29_2_TA" based on "Score"
bin.size <- 3
score_bins <- seq(0, max(m.contexts.all.OneTo18$Score, na.rm = TRUE), by = bin.size)
fractions <- data.frame(
    ScoreBin = score_bins[-length(score_bins)],
    Fraction_DN = sapply(1:(length(score_bins) - 1), function(i) {
        bin <- m.contexts.all.OneTo18$Score >= score_bins[i] & m.contexts.all.OneTo18$Score < score_bins[i + 1]
        sum(m.contexts.all.OneTo18$"tp73_skmel29_2_DN"[bin] > 0, na.rm = TRUE) / sum(bin, na.rm = TRUE)
    }),
    Fraction_TA = sapply(1:(length(score_bins) - 1), function(i) {
        bin <- m.contexts.all.OneTo18$Score >= score_bins[i] & m.contexts.all.OneTo18$Score < score_bins[i + 1]
        sum(m.contexts.all.OneTo18$"tp73_skmel29_2_TA"[bin] > 0, na.rm = TRUE) / sum(bin, na.rm = TRUE)
    }),
      # Add bars to indicate the number of TFBS within each bin
    TFBS_Count = sapply(1:(length(score_bins) - 1), function(i) {
        bin <- m.contexts.all.OneTo18$Score >= score_bins[i] & m.contexts.all.OneTo18$Score < score_bins[i + 1]
        sum(bin, na.rm = TRUE)
    })
)

genome.wide.quantiles.50 <- array(
    unlist(lapply(isoforms, function(isoform) {
        lapply(cell.lines, function(cell.line) {
            lapply(targets, function(target) {
                col.name <- paste0(target, "_", cell.line, "_", isoform)
                if (!col.name %in% names(m.contexts.all.OneTo18)) {
                    stop(paste("Column", col.name, "not found in m.contexts.all.OneTo18"))
                }
                max(quantile(m.contexts.all.OneTo18[[col.name]], probs=0.50, na.rm=TRUE), 0.5)
            })
        })
    })),
    dim = c(length(targets), length(cell.lines), length(isoforms)),
    dimnames = list(targets=targets, cell.lines=cell.lines, isoforms=isoforms)
)

genome.wide.quantiles.90 <- array(
    unlist(lapply(isoforms, function(isoform) {
        lapply(cell.lines, function(cell.line) {
            lapply(targets, function(target) {
                col.name <- paste0(target, "_", cell.line, "_", isoform)
                if (!col.name %in% names(m.contexts.all.OneTo18)) {
                    stop(paste("Column", col.name, "not found in m.contexts.all.OneTo18"))
                }
                max(quantile(m.contexts.all.OneTo18[[col.name]], probs=0.90, na.rm=TRUE), 0.5)
            })
        })
    })),
    dim = c(length(targets), length(cell.lines), length(isoforms)),
    dimnames = list(targets=targets, cell.lines=cell.lines, isoforms=isoforms)
)

genome.wide.quantiles.95 <- array(
    unlist(lapply(isoforms, function(isoform) {
        lapply(cell.lines, function(cell.line) {
            lapply(targets, function(target) {
                col.name <- paste0(target, "_", cell.line, "_", isoform)
                if (!col.name %in% names(m.contexts.all.OneTo18)) {
                    stop(paste("Column", col.name, "not found in m.contexts.all.OneTo18"))
                }
                max(quantile(m.contexts.all.OneTo18[[col.name]], probs=0.95, na.rm=TRUE), 0.5)
            })
        })
    })),
    dim = c(length(targets), length(cell.lines), length(isoforms)),
    dimnames = list(targets=targets, cell.lines=cell.lines, isoforms=isoforms)
)

if (file.exists("ratios.methylated.RData")) {

    load("ratios.methylated.RData")

} else {
 
    compute.ratios.methylated.within.promoters <- function(chromosome, genome.wide.quantiles=genome.wide.quantiles.50,isoforms=c("TA","DN","GFP","any"), cell.lines=c("saos2","skmel29_2")) {
        results <- list()
        cols.NumInWindow <- grepl("_NumInWindow", colnames(m.contexts[[chromosome]]))
        for(cell.line in cell.lines) {
            for(isoform in isoforms) {
                cat("I: Computing ratios for cell line '", cell.line, "' and isoform '", isoform, "' on chromosome '",chromosome,"'.\n", sep="")
                
                if ("any" == isoform) {

                    # If isoform is "any", we need to handle it differently

                    tp73.col.names.GFP <- paste0("tp73_", cell.line, "_", "GFP")
                    tp73.col.names.DN <- paste0("tp73_", cell.line, "_", "DN")
                    tp73.col.names.TA <- paste0("tp73_", cell.line, "_", "TA")

                    pos.col.name.GFP <- paste0("pos_", cell.line, "_", "GFP")
                    pos.col.name.DN <- paste0("pos_", cell.line, "_", "DN")
                    pos.col.name.TA <- paste0("pos_", cell.line, "_", "TA")

                    # Validate column existence
                    if (!tp73.col.names.GFP %in% names(m.contexts[[chromosome]])) {
                        stop(paste("Column", tp73.col.names.GFP, "not found in m.contexts[[",chromosome,"]]"))
                    }
                    if (!tp73.col.names.DN %in% names(m.contexts[[chromosome]])) {
                        stop(paste("Column", tp73.col.names.DN, "not found in m.contexts[[",chromosome,"]]"))
                    }
                    if (!tp73.col.names.TA %in% names(m.contexts[[chromosome]])) {
                        stop(paste("Column", tp73.col.names.TA, "not found in m.contexts[[",chromosome,"]]"))
                    }
                    if (!pos.col.name.GFP %in% names(m.contexts[[chromosome]])) {
                        stop(paste("Column", pos.col.name.GFP, "not found in m.contexts[[",chromosome,"]]"))
                    }
                    if (!pos.col.name.DN %in% names(m.contexts[[chromosome]])) {
                        stop(paste("Column", pos.col.name.DN, "not found in m.contexts[[",chromosome,"]]"))
                    }
                    if (!pos.col.name.TA %in% names(m.contexts[[chromosome]])) {
                        stop(paste("Column", pos.col.name.TA, "not found in m.contexts[[",chromosome,"]]"))
                    }
                
                    methylated.within.promoter <- sapply(m.contexts[[chromosome]][, ..cols.NumInWindow], function(X) {
                        max(sum(X[
                            (
                                m.contexts[[chromosome]][[tp73.col.names.GFP]] >= genome.wide.quantiles["tp73", cell.line, "GFP"] |
                                m.contexts[[chromosome]][[tp73.col.names.DN]] >= genome.wide.quantiles["tp73", cell.line, "DN"] |
                                m.contexts[[chromosome]][[tp73.col.names.TA]] >= genome.wide.quantiles["tp73", cell.line, "TA"]
                            ) & (
                                m.contexts[[chromosome]][[pos.col.name.GFP]] >= genome.wide.quantiles["pos", cell.line, "GFP"] |
                                m.contexts[[chromosome]][[pos.col.name.DN]] >= genome.wide.quantiles["pos", cell.line, "DN"] |
                                m.contexts[[chromosome]][[pos.col.name.TA]] >= genome.wide.quantiles["pos", cell.line, "TA"] 
                            ) &
                            m.contexts[[chromosome]]$InPromoter == TRUE
                        ]>0),0.1)
                    }) / sum(m.contexts[[chromosome]]$InPromoter == TRUE)

                    # Compute not methylated outside promoter ratio
                    not.methylated.outside.promoter <- sapply(m.contexts[[chromosome]][, ..cols.NumInWindow], function(X) {
                        max(sum(X[
                            (
                                m.contexts[[chromosome]][[tp73.col.names.GFP]] < genome.wide.quantiles["tp73", cell.line, "GFP"] &
                                m.contexts[[chromosome]][[tp73.col.names.DN]] < genome.wide.quantiles["tp73", cell.line, "DN"] &
                                m.contexts[[chromosome]][[tp73.col.names.TA]] < genome.wide.quantiles["tp73", cell.line, "TA"]
                            ) & (
                                m.contexts[[chromosome]][[pos.col.name.GFP]] < genome.wide.quantiles["pos", cell.line, "GFP"] &
                                m.contexts[[chromosome]][[pos.col.name.DN]] < genome.wide.quantiles["pos", cell.line, "DN"] &
                                m.contexts[[chromosome]][[pos.col.name.TA]] < genome.wide.quantiles["pos", cell.line, "TA"] 
                            ) &
                            m.contexts[[chromosome]]$InPromoter == FALSE
                        ]>0),0.1)
                    }) / sum(m.contexts[[chromosome]]$InPromoter == FALSE)
                    
                } else {
                    tp73.col.name <- paste0("tp73_", cell.line, "_", isoform)
                    pos.col.name <- paste0("pos_", cell.line, "_", isoform)
                    # Validate column existence
                    if (!tp73.col.name %in% names(m.contexts[[chromosome]])) {
                        stop(paste("Column", tp73.col.name, "not found in m.contexts[[",i,"]]"))
                    }
                    if (!pos.col.name %in% names(m.contexts[[chromosome]])) {
                        stop(paste("Column", pos.col.name, "not found in m.contexts[[",i,"]]"))
                    }

                    # Compute methylated within promoter ratio
                    methylated.within.promoter <- sapply(m.contexts[[chromosome]][, ..cols.NumInWindow], function(X) {
                        max(sum(X[
                            m.contexts[[chromosome]][[tp73.col.name]] >= genome.wide.quantiles["tp73", cell.line, isoform] &
                            m.contexts[[chromosome]][[pos.col.name]] >= genome.wide.quantiles["pos", cell.line, isoform] &
                            m.contexts[[chromosome]]$InPromoter == TRUE
                        ]>0),0.1)
                    }) / sum(m.contexts[[chromosome]]$InPromoter == TRUE)

                    # Compute not methylated outside promoter ratio
                    not.methylated.outside.promoter <- sapply(m.contexts[[chromosome]][, ..cols.NumInWindow], function(X) {
                        max(sum(X[
                            m.contexts[[chromosome]][[tp73.col.name]] < genome.wide.quantiles["tp73", cell.line, isoform] &
                            m.contexts[[chromosome]][[pos.col.name]] < genome.wide.quantiles["pos", cell.line, isoform] &
                            m.contexts[[chromosome]]$InPromoter == FALSE
                        ]>0), 0.1)
                    }) / sum(m.contexts[[chromosome]]$InPromoter == FALSE)

                }
                #methylated.within.promoter[0 == methylated.within.promoter] <- 0.001  # Avoid division by zero, keeping chance to reach value of 1           
                #not.methylated.outside.promoter[0 == not.methylated.outside.promoter] <- 0.001  # Avoid division by zero

                # Store results
                results[[paste0(cell.line, "_", isoform)]] <- methylated.within.promoter / not.methylated.outside.promoter
                print(summary(results[[paste0(cell.line, "_", isoform)]]))
            }
        }
        results
    }

    ratios.methylated <- lapply(names(m.contexts),compute.ratios.methylated.within.promoters,genome.wide.quantiles=genome.wide.quantiles.50)
    names(ratios.methylated) <- names(m.contexts)
    save(ratios.methylated, file="ratios.methylated.RData")
}

# Summarize data across chromosomes for each vector position using R's summary function
summarize_across_chromosomes <- function(ratios.methylated, cell_line_isoform = NULL) {
    # Extract names from the first chromosome's data
    vector_names <- NULL
    if (!is.null(cell_line_isoform)) {
        vector_names <- names(ratios.methylated[[1]][[cell_line_isoform]])
    } else {
        if (is.matrix(ratios.methylated[[1]])) {
            vector_names <- colnames(ratios.methylated[[1]])
        } else {
            vector_names <- names(ratios.methylated[[1]])
        }
    }
    
    # Initialize a list to store summaries for each vector position
    summary_list <- list()
    
    for (name in vector_names) {
        # Collect values across all chromosomes for the current vector position
        values <- NULL
        if (!is.null(cell_line_isoform)) {
            values <- sapply(1:22, function(chromosome) ratios.methylated[[chromosome]][[cell_line_isoform]][[name]])
        } else {
            values <- sapply(1:22, function(chromosome) ratios.methylated[[chromosome]][[name]])
        }
        # Use R's summary function to compute summary statistics
        summary_list[[name]] <- c(summary(values),SD=sd(values))
    }
    
    summary_list
}

# Generate summary for "skmel29_2_any"

nicer.row.names <- function(X) {
    X <- strsplit(X, split="_")
    X <- sapply(X, function(x) {
        if (length(x) >= 2) {
            paste0(x[1]," ","(",x[2],")")
        } else {
            x
        }
    })
}

tf.of.interest <- grep(rownames(summary_data),pattern="^(jun|sp1|rest|yap1|yy1|hey1) ",ignore.case=T,value=T)

summary_data <- summarize_across_chromosomes(ratios.methylated, cell_line_isoform = "skmel29_2_any")
rownames(summary_data) <- nicer.row.names(rownames(summary_data))

# ordering by coefficient of variation (mean/SD)
#summary_data <- summary_data[order(summary_data[,"Mean"]/summary_data[,"SD"]),] 
#summary_data <- summary_data[order(summary_data[,"Mean"]),] 
summary_data <- summary_data[order(summary_data[,"Median"]),] 


head(round(summary_data,3), 30)
cbind(round(which(rownames(summary_data) %in% tf.of.interest)/nrow(summary_data)*100),round(summary_data[tf.of.interest,],3))
tail(round(summary_datam3), 50)

# Print summary
print(summary_data)

# Derive ratios between isoforms for each chromosome
ratio.skmel29_2.TAa.vs.DNb  <- sapply(names(ratios.methylated), function(chromosome) {
    ratios.methylated[[chromosome]][["skmel29_2_TA"]] / ratios.methylated[[chromosome]][["skmel29_2_DN"]]
})

ratio.summary_data <- t(apply(ratio.skmel29_2.TAa.vs.DNb,1,function(X) {
    c(summary(X),SD=sd(X))
}))
rownames(ratio.summary_data) <- nicer.row.names(rownames(ratio.summary_data))
ratio.summary_data <- ratio.summary_data[order(ratio.summary_data[,"Median"]),] 


ratio.skmel29_2.TAa.vs.DNb.log2 <- log2(ratio.skmel29_2.TAa.vs.DNb)
ratio.summary_data.log2 <- t(apply(ratio.skmel29_2.TAa.vs.DNb.log2,1,function(X) {
    c(summary(X),SD=sd(X))
}))
rownames(ratio.summary_data.log2) <- nicer.row.names(rownames(ratio.summary_data.log2))
ratio.summary_data.log2 <- ratio.summary_data.log2[order(ratio.summary_data.log2[,"Median"]),] 
head(round(ratio.summary_data.log2,2), 30)
cbind(round(which(rownames(ratio.summary_data.log2) %in% tf.of.interest)/nrow(ratio.summary_data.log2)*100),round(ratio.summary_data.log2[tf.of.interest,],2))
tail(round(ratio.summary_data.log2,2), 30)


