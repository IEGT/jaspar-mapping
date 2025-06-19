
require(ggplot2)

plot.series.of.ratio.matrices <- function(list.of.sorted.matrices.to.display=NULL,num.from.extrema=rep(10,length(list.of.sorted.matrices.to.display)),tf.patterns.of.interest=NULL,human.only=TRUE) {

    if (is.null(tf.patterns.of.interest)) {
        #tf.patterns.of.interest <- "^(jun|sp1|rest|yap1|yy1|e2f1|e2f7|tp53|tp63|tp73|plag1|rara|odd|irf2|bpc5|irf2|nr2f6|pdr3|pax1|rarg|esr1|odd|elt-3|nfkb1|rela|rora|erf6) "
        tf.patterns.of.interest <- paste0(
            #"^(jun|sp1|rest|yap1|yy1|e2f1|nfkb1|tp53|tp63|tp73|",
            #"BZIP43|IRF2|RXRA--VDR|pan|HAP1|Nr2F6|RORA|E2F3|E2F7|THI2|SPL15|TGIF1|MYB52|ZBED1|NAC083", #  paste(sapply(strsplit(rownames(head(ta_vs_effects_tp73Confirm_posConfirm.sorted,15)),split=" "),function(X) X[1]),collapse="|")
            #"NR3C1|TP63|RREB1|ATHB-9|GLIS1|PLAG1|E2FC|MYB15|DYT1|bZIP911|odd|RARA|BZIP42|BPC5|TCP14", # paste(sapply(strsplit(rownames(tail(ta_vs_effects_tp73Confirm_posConfirm.sorted[!is.na(ta_vs_effects_tp73Confirm_posConfirm.sorted[,"ratio"]),],15)),split=" "),function(X) X[1]),collapse="|")
            "RARA|GLIS1|PLAG1|RREB1|TP63|NR3C1|SRF|NR3C2|EWSR1-FLI1|THRB|TFAP2E|GLIS3|TFAP2C|RXRB|KLF14|",
            "IRF2|RXRA|RORA|E2F7|E2F3|TGIF1|ZBED1|IRF4|RFX3|MEF2D|STAT1|MEF2B|NFIC|IRF8|CREB3L1",
        ") ")
    }

    if (is.null(list.of.sorted.matrices.to.display)) {
        cat("I: No sorted matrix provided, using default sorted matrices.\n")
        list.of.sorted.matrices.to.display <- list(
                effect_nonmethylated=ta_vs_dn_nonMethylated.sorted,
                effect_methylated=ta_vs_dn_w_methylation.sorted
        )
    } else if (!is.list(list.of.sorted.matrices.to.display)) {
        if (is.matrix(list.of.sorted.matrices.to.display) || is.data.frame(sorted.matrix.to.display)) {
            list.of.sorted.matrices.to.display <- list(unnamed=sorted.matrix.to.display)
        } else {
            cat("E: Invalid input for sorted.matrix.to.display, expected a matrix, data frame, or list.\n")
            stop("Exiting due to invalid input.")
        }
        cat("I: Using provided sorted matrix list.\n")
    }

    if (0 == length(list.of.sorted.matrices.to.display)) {
        cat("E: No sorted matrices to display, exiting.\n")
        stop("Exiting due to no sorted matrices to display.")
    } else {
           cat("I: Found ", length(list.of.sorted.matrices.to.display), " sorted matrices to display:\n", sep="")
           cat("  ", paste(names(list.of.sorted.matrices.to.display), collapse=", "), "\n", sep="")
    }

    pdf("vertical_plot.pdf", width=3*length(list.of.sorted.matrices.to.display)+1, height=8)

    for(i in 1:length(list.of.sorted.matrices.to.display)) {

        sorted.matrix.to.display.tag <- names(list.of.sorted.matrices.to.display)[i]
        add.to.existing.graph <- ( i > 1 )

        cat("I: Processing sorted matrix for tag: ", sorted.matrix.to.display.tag, "\n", sep="")

        print_legend <- (i == length(list.of.sorted.matrices.to.display))

        if (add.to.existing.graph) {
            cat("I: Adding to existing graph for tag: ", sorted.matrix.to.display.tag, "\n", sep="")
        } else {
            cat("I: Creating new graph for tag: ", sorted.matrix.to.display.tag, "\n", sep="")
        }

        sorted.matrix.to.display <- list.of.sorted.matrices.to.display[[sorted.matrix.to.display.tag]]
        print(dim(sorted.matrix.to.display))
        sorted.matrix.to.display <- sorted.matrix.to.display[which(is.human.jaspar.id(rownames(sorted.matrix.to.display))),,drop=FALSE]
        print(dim(sorted.matrix.to.display))

        sorted.matrix.to.display <- sorted.matrix.to.display[order(sorted.matrix.to.display[, "ratio"],decreasing=TRUE), , drop=FALSE]

        if (is.null(sorted.matrix.to.display)) {
            cat("E: No data for tag ", sorted.matrix.to.display.tag, " - skipping\n", sep="")
            next
        }

        sorted.matrix.to.display.valid <- rowSums(sorted.matrix.to.display[, c("mean.TA","mean.DN")])>0.0005

        # Get top and bottom
        top5 <- head(sorted.matrix.to.display[sorted.matrix.to.display.valid,], num.from.extrema[i])
        bottom5 <- tail(sorted.matrix.to.display[sorted.matrix.to.display.valid,], num.from.extrema[i])

        # Genes of interest (tf.of.interest) that are not in top5/bottom5
        tf_interest_names <- grep(rownames(ta_vs_dn_nonMethylated.sorted),pattern=tf.patterns.of.interest,ignore.case=T,value=T)
        print(tf_interest_names)

        tf_interest_set <- setdiff(tf_interest_names, c(rownames(top5), rownames(bottom5)))
        tf_interest <- sorted.matrix.to.display[rownames(sorted.matrix.to.display) %in% tf_interest_set, , drop=FALSE] # extra complicated to preserve order

        # Assign y positions: top5 (1:5), gap, tf_interest (by their rank), gap, bottom5 (last 5)
        n_top <- nrow(top5)
        n_bottom <- nrow(bottom5)
        gap_size <- 2

        # Get the ranks of tf_interest in the full sorted table
        tf_interest_ranks <- match(rownames(tf_interest), rownames(sorted.matrix.to.display))

        # Compose plotting data
        plot_data <- rbind(
            data.frame(Name=rownames(top5), mean.TA=top5[,1], mean.DN=top5[,2], ratio=top5[,3], group="Top", y=seq(1,n_top)),
            data.frame(Name=rownames(tf_interest), mean.TA=tf_interest[,1], mean.DN=tf_interest[,2], ratio=tf_interest[,3], group="Interest", y=tf_interest_ranks),
            data.frame(Name=rownames(bottom5), mean.TA=bottom5[,1], mean.DN=bottom5[,2], ratio=bottom5[,3], group="Bottom",
                        y=seq(sum(sorted.matrix.to.display.valid)-n_bottom+1,sum(sorted.matrix.to.display.valid)))
        )

        # Assign y for genes of interest: scale to fit in the gap
        gap_start <- 0
        gap_end <- -gap_size
        if(nrow(tf_interest)>0) {
            # Scale ranks to fit in the gap
            plot_data$y[plot_data$group=="Expression/Literature"] <- gap_start - (tf_interest_ranks - min(tf_interest_ranks)) /
                (max(tf_interest_ranks)-min(tf_interest_ranks)+1) * gap_size
        }

        # Add a column for circle size (sum of means)
        plot_data$circle_size <- plot_data$mean.TA + plot_data$mean.DN

        # Assign y positions with reserved space: 20% for top5, 20% for bottom5, 60% for interest
        n_total <- sum(sorted.matrix.to.display.valid)
        n_top <- nrow(top5)
        n_bottom <- nrow(bottom5)
        n_interest <- nrow(tf_interest)

        # Define y-axis range
        y_min <- 0
        y_max <- 1

        # Allocate space
        top_space <- num.from.extrema[i]/50
        bottom_space <- num.from.extrema[i]/50
        interest_space <- 1-top_space-bottom_space

        # Calculate y positions for each group
        if (n_top > 1) {
            y_top <- seq(y_max, y_max - top_space, length.out = n_top + 1)[-1]
        } else {
            y_top <- y_max - top_space / 2
        }
        if (n_bottom > 1) {
            y_bottom <- seq(y_min + bottom_space, y_min , length.out = n_bottom+1)[-1]
        } else {
            y_bottom <- y_min + bottom_space / 2
        }
        if (n_interest > 0) {
            y_interest <- seq(y_max - top_space, y_min + bottom_space, length.out = n_interest + 2)[-c(1, n_interest + 2)]
        } else {
            y_interest <- numeric(0)
        }

        # Assign y to plot_data
        plot_data$y <- NA
        plot_data$y[plot_data$group == "Top"] <- y_top
        plot_data$y[plot_data$group == "Bottom"] <- y_bottom
        plot_data$y[plot_data$group == "Interest"] <- y_interest

        # Plot
        # Filter out rows with NA y values to avoid warnings and missing labels
        plot_data_valid <- plot_data[!is.na(plot_data$y), ]

        # Instead of using labs(title=...), add the tag as a text annotation at the correct x/y position
        # Place the label at the center top of the current column (x = (i-1)*2, y = max y + offset)
        p <- ggplot(plot_data_valid, aes(x=(i-1)*1.5, y=y)) +
            geom_point(aes(size=circle_size, fill=group), shape=21, color="black", alpha=0.8) +
            geom_text(aes(label=Name), hjust=0, nudge_x=0.12, size=2) +
            geom_text(aes(label=round(log2(ratio),2)), hjust=1, vjust=0.5, size=2, nudge_x=-0.12, fontface="bold") +
            scale_size_continuous(breaks = c(0.05,0.2,1,5), labels = c("0.05","0.2","1","5")) +
            scale_fill_manual(values=c("Top"="#1b9e77", "Interest"="#d95f02", "Bottom"="#7570b3")) +
            theme_minimal() +
            labs(size="Frequency", fill="Group") +
            theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
            panel.grid=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),
            legend.position=if (exists("print_legend") && print_legend) "right" else "none") +
            xlim(c(-0.75, (length(list.of.sorted.matrices.to.display)-1)*1.5+0.75)) +
            ylim(min(plot_data_valid$y, na.rm=TRUE), max(plot_data_valid$y, na.rm=TRUE)+0.2) +
            annotate("text", x=0.5+(i-1)*1.5, y=max(plot_data_valid$y, na.rm=TRUE)+0.1, na.rm=TRUE,
                 label=sorted.matrix.to.display.tag, fontface="bold", size=4, hjust=0.5, vjust=0)

        if (add.to.existing.graph) {
            print(p, newpage=FALSE)
        } else {
            print(p)
        }
        add.to.existing.graph <- TRUE

    }

    dev.off()
    cat("I: Vertical plot saved to 'vertical_plot.pdf'.\n")
}


list.of.matrices <- list(
    "p73 binding"=ta_vs_dn_tp73Confirm.sorted,
    "in promoter"= ta_vs_dn_tp73Confirm_inPromoter.sorted,
    "also methylated"=ta_vs_dn_tp73Confirm_inPromoter_posConfirm.sorted,
    # "effect nonmethylated"=ta_vs_dn_nonMethylated.sorted,
    #"effects expression"=ta_vs_dn_w_methylation.sorted
    "effects expression"=ta_vs_effects_tp73Confirm_posConfirm.sorted,
    #"geneset EMT"=ta_vs_dn_genesetEMT_tp73Confirm_posConfirm.sorted
    "TA/DN-induced EMT"=ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm.sorted
)

plot.series.of.ratio.matrices(
    list.of.sorted.matrices.to.display=list.of.matrices,
    num.from.extrema=c(10,10,10,15,15), # Adjusted to match the number of matrices
)
