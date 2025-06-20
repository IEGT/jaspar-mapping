
#
# Heatmap for gene selection
#

# Retrieves all contexts for binding sites in promoter regions, which is the background for the analysis.
all_inPromoter_tp73ConfirmAny <- retrieve_context_data_by_chromosome(NULL, confirmation=c("tp73","promoter"),TA.or.DN="any")
all_inPromoter_tp73ConfirmAny.colSums <- colSums(all_inPromoter_tp73ConfirmAny$context_data_binary,na.rm=T)

require(ggplot2)

#gene.selection <- promoterBedTables$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.promoter.bed$Gene
#gene.selection <- gene.selection.EMT.TAorDN

#gene.selection <- e.ta.up
#gene.selection <- e.dn.up

relevant.cofactors.ta.names <- names(which(relevant.cofactors.ta))
relevant.cofactors.dn.names <- names(which(relevant.cofactors.dn))

my.heatmap(target.gene.selection.TA=e.ta.up, target.gene.selection.DN=e.dn.up) 

my.heatmap.and.corrplot <- function(base.filename="gene_selection_tp73ConfirmAny_unnamed",
                                    formats=c("pdf","png","svg"),
                                    target.gene.selection.TA, target.gene.selection.DN, target.gene.selection.always=c("TP73"),
                                    relevant.cofactors.selection.ta=relevant.cofactors.ta.names,
                                    relevant.cofactors.selection.dn=relevant.cofactors.dn.names,
                                    relevant.cofactors.selection.always=c("TP73 (MA0861.1)"),
                                    combined.expression.data = combined.expression.data,
                                    universe.prior.to.selection = all_inPromoter_tp73ConfirmAny,
                                    threshold.min.frequency = 0.05,
                                    threshold.min.log2.enrichment = 0.1,
                                    threshold.min.log2.depletion = -threshold.min.log2.enrichment,
                                    #threshold.min.log2.depletion = NA,
                                    filter.for.human.jaspar=TRUE,
                                    filter.for.strong.enrichment.num=NA,
                                    show.cutandrun=TRUE,
                                    show.score=TRUE,
                                    combined.expression.data.selected=combined.expression.data[,"Gene Symbol"] %in% c(target.gene.selection.TA, target.gene.selection.DN),
                                    plot.corrplot=TRUE,
                                    plot.heatmap=TRUE,
                                    verbose=TRUE) {

    stopifnot(all(formats %in% c("pdf","png","svg")))

    target.gene.selection_tp73ConfirmAny <- retrieve_context_data_by_chromosome(combined.expression.data.selected,
                                                                         confirmation=c("tp73"),TA.or.DN="any")
    if (verbose) {
        cat("I: Distribution of number of matches by chromosome:\n")
        print(target.gene.selection_tp73ConfirmAny$num_matches_by_chromosome)
        cat("I: Total number of matches across all chromosomes: ",sum(target.gene.selection_tp73ConfirmAny$num_matches_by_chromosome),"\n",sep="")
    }

    target.gene.selection_tp73ConfirmAny.colSums <- colSums(target.gene.selection_tp73ConfirmAny$context_data_binary,na.rm=T)
    print(target.gene.selection_tp73ConfirmAny.colSums)

    target.gene.vs.universe.prior.to.selection <- cbind(
        mean.enriched=target.gene.selection_tp73ConfirmAny$mean_total_binary,
        mean.background=universe.prior.to.selection$mean_total_binary,
        log2.ratio=log2(target.gene.selection_tp73ConfirmAny$mean_total_binary / universe.prior.to.selection$mean_total_binary)
    )

    if (FALSE) {
        tf.patterns.of.interest <- paste0(
                    #"^(jun|sp1|rest|yap1|yy1|e2f1|nfkb1|tp53|tp63|tp73|",
                    #"BZIP43|IRF2|RXRA--VDR|pan|HAP1|Nr2F6|RORA|E2F3|E2F7|THI2|SPL15|TGIF1|MYB52|ZBED1|NAC083", #  paste(sapply(strsplit(rownames(head(ta_vs_effects_tp73Confirm_posConfirm.sorted,15)),split=" "),function(X) X[1]),collapse="|")
                    #"NR3C1|TP63|RREB1|ATHB-9|GLIS1|PLAG1|E2FC|MYB15|DYT1|bZIP911|odd|RARA|BZIP42|BPC5|TCP14", # paste(sapply(strsplit(rownames(tail(ta_vs_effects_tp73Confirm_posConfirm.sorted[!is.na(ta_vs_effects_tp73Confirm_posConfirm.sorted[,"ratio"]),],15)),split=" "),function(X) X[1]),collapse="|")
                    "RARA|GLIS1|PLAG1|RREB1|TP63|NR3C1|SRF|NR3C2|EWSR1-FLI1|THRB|TFAP2E|GLIS3|TFAP2C|RXRB|KLF14|",
                    "IRF2|RXRA|RORA|E2F7|E2F3|TFIF1|ZBED1|IRF4|RFX3|MEF2D|STAT1|MEF2B|NFIC|IRF8|CREB3L1",
                ") ")
    }

    # Ranking genes for their ratio (== enrichment)
    target.gene.vs.universe.prior.to.selection.sorted <- target.gene.vs.universe.prior.to.selection[order(target.gene.vs.universe.prior.to.selection[, "log2.ratio"],decreasing=T), ]
    if (verbose) cat("I: Number of target genes vs universe prior to selection: ", nrow(target.gene.vs.universe.prior.to.selection.sorted), "\n")
    rownames(target.gene.vs.universe.prior.to.selection.sorted) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection.sorted))

    if (filter.for.human.jaspar) {
        if (verbose) cat("I: Filtering for human JASPAR IDs...\n")
        target.gene.vs.universe.prior.to.selection.sorted <- target.gene.vs.universe.prior.to.selection.sorted[is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection.sorted)),]
        if (verbose) cat("I: Number of relevant cofactors remaining: ", nrow(target.gene.vs.universe.prior.to.selection.sorted), "\n")
    }

    if (!is.na(threshold.min.log2.enrichment) && threshold.min.log2.enrichment > 0) {
        if (verbose) cat("I: Filtering for log2 ratio > ",threshold.min.log2.enrichment," if above 0\n",sep="")
        target.gene.vs.universe.prior.to.selection.sorted <- target.gene.vs.universe.prior.to.selection.sorted[
            target.gene.vs.universe.prior.to.selection.sorted[,"log2.ratio"] < 0 |
            target.gene.vs.universe.prior.to.selection.sorted[,"log2.ratio"] > threshold.min.log2.enrichment,]
        if (verbose) cat("I: Number of relevant cofactors remaining: ", nrow(target.gene.vs.universe.prior.to.selection.sorted), "\n")
    }    

    if (!is.na(threshold.min.log2.depletion) && threshold.min.log2.depletion < 0) {
        if (verbose) cat("I: Filtering for log2 ratio < ",threshold.min.log2.depletion," if below 0\n",sep="")
        target.gene.vs.universe.prior.to.selection.sorted <- target.gene.vs.universe.prior.to.selection.sorted[
            target.gene.vs.universe.prior.to.selection.sorted[,"log2.ratio"] > 0 |
            target.gene.vs.universe.prior.to.selection.sorted[,"log2.ratio"] < threshold.min.log2.depletion,]
        if (verbose) cat("I: Number of relevant cofactors remaining: ", nrow(target.gene.vs.universe.prior.to.selection.sorted), "\n")
    }

    if(!is.na(threshold.min.frequency) && threshold.min.frequency > 0) {
        if (verbose) cat("I: Filtering for mean background > ",threshold.min.frequency,".\n",sep="")
        target.gene.vs.universe.prior.to.selection.sorted <- target.gene.vs.universe.prior.to.selection.sorted[
            target.gene.vs.universe.prior.to.selection.sorted[,"mean.enriched"] > threshold.min.frequency,]
        if (verbose) cat("I: Number of relevant cofactors remaining: ", nrow(target.gene.vs.universe.prior.to.selection.sorted), "\n")
    }

    if (!is.na(filter.for.strong.enrichment.num) && filter.for.strong.enrichment.num > 0) {
        if (verbose) cat("I: Filtering for strong enrichment > ",filter.for.strong.enrichment.num,".\n",sep="")
        target.gene.vs.universe.prior.to.selection.sorted <- target.gene.vs.universe.prior.to.selection.sorted[1:filter.for.strong.enrichment.num,]
        if (verbose) cat("I: Number of relevant cofactors remaining: ", nrow(target.gene.vs.universe.prior.to.selection.sorted), "\n")
    }

    cofactor.names <- unique(c(relevant.cofactors.selection.dn,relevant.cofactors.selection.always,
                               rownames(target.gene.vs.universe.prior.to.selection.sorted),
                               relevant.cofactors.selection.ta))

    #?
    #combined.expression.data.reduced <- combined.expression.data[combined.expression.data.selected,]
    #nrow(target.gene.selection_tp73ConfirmAny$context_data)

# Print all possible six combinations of "pos"/"tp73" with "skmel29_2" and "TA"/"DN"/"GFP"
    combinations <- sort(as.vector(outer(c("pos", "tp73"), c("TA", "DN", "GFP"), function(x, y) paste(x, "skmel29_2", y, sep = "_"))))
    #print(combinations)

    target.gene.selection_tp73ConfirmAny_context <- target.gene.selection_tp73ConfirmAny$context_data_binary
    colnames(target.gene.selection_tp73ConfirmAny_context) <- prettyIdentifierJaspar(colnames(target.gene.selection_tp73ConfirmAny_context))
    target.gene.selection_tp73ConfirmAny_context_interest <- target.gene.selection_tp73ConfirmAny_context[, cofactor.names, drop = FALSE]
    target.gene.selection_tp73ConfirmAny_context_interest.colsums <- colSums(target.gene.selection_tp73ConfirmAny_context_interest)
    target.gene.selection_tp73ConfirmAny_context_interest.informativeColumNames <- colnames(target.gene.selection_tp73ConfirmAny_context_interest)[target.gene.selection_tp73ConfirmAny_context_interest.colsums > 0]
    target.gene.selection_tp73ConfirmAny_context_interest.informative <- target.gene.selection_tp73ConfirmAny_context_interest[,target.gene.selection_tp73ConfirmAny_context_interest.informativeColumNames, drop = FALSE]


    if (plot.corrplot) {    
        # Plot correlation matrix of all transcription factors (columns) against each other

        require(corrplot)

        m <- target.gene.selection_tp73ConfirmAny_context_interest.informative
        m.sd <- apply(m,2,sd)
        m <- m[,m.sd>0]

        # Use the context data for all promoters with TP73 confirmation as the matrix
        cor_matrix <- cor(m, use="pairwise.complete.obs", method="pearson")

        pdf("correlation_matrix_TFBS.pdf", width=2+9*ncol(m)/60, height=1 + 9*ncol(m)/60)
        #op <- par(mar=c(0,0,2,0)) # Set margins for the plot
        corrplot(cor_matrix, is.corr=TRUE, method="color", type="upper", order="hclust",
                tl.col= 1 + colnames(m) %in% relevant.cofactors.selection.dn + 
                    2*colnames(m) %in% relevant.cofactors.selection.ta + 
                    4*colnames(m) %in% relevant.cofactors.selection.always,
                tl.cex=0.7, number.cex=0.35, number.digits=1,insig="n",
                addCoef.col=NULL,
                title="Correlation Matrix of TFBS (all promoters, TP73 confirmed)")
        legend("bottomleft", legend=c("enriched DN","enriched TA","enriched DN+TA","newly found in data"), fill=c(1+1,1+2,1+2+4,1), bty="n", cex=0.9, inset=0.02, xpd=TRUE)
        #par(op) # Restore previous margins
        dev.off()
        cat("I: Correlation matrix plot saved to 'correlation_matrix_TFBS.pdf'.\n")
    }
    debugme <- data.frame(colnames(m),col=1 + colnames(m) %in% relevant.cofactors.selection.dn +  2*colnames(m) %in% relevant.cofactors.selection.ta +  4*colnames(m) %in% relevant.cofactors.selection.always)
    target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input <- cbind(
        Gene=target.gene.selection_tp73ConfirmAny$downstream_genes,
        target.gene.selection_tp73ConfirmAny_context_interest.informative
    )
    if (show.score) {
        target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input <- cbind(
            target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input,
            Score=gene.selection_tp73ConfirmAny$coordinates[,c(5)]
        )
    }

    if (show.cutandrun) {
        target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input <- cbind(
            target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input,
            target.gene.selection_tp73ConfirmAny$cutandrun_data[,..combinations] 
        )
    }

    # FIXME: Is this avoidable
    target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant <- target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input[!duplicated(target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input[,"Gene"]),]
    rownames(target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant) <- target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant$Gene

    require(gplots)

    n <- target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant$Gene
    n.show <- as.matrix(target.gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant[,-c(1)])

    if (FALSE) {
        n.show <- apply(n.show, 2, function(x) {
            x[is.na(x)] <- 0
            x = x/max(x, na.rm=TRUE)  # Scale to [0, 1]
            return(x)
        })
    }
    if (FALSE) {
        n.show <- apply(n.show, c(1,2), function(x) {
            if (x>0.33) x=1
            else x=0
            x
        })
    }

    for (f in formats) {
        file <- paste0(base.filename,".", f)
        if (f == "png") {
            png(file)
        } else if (f == "pdf") {
            pdf(file, width=round(18*ncol(n.show)/50,0), height=round(8*nrow(n.show)/40,0))
        } else if (f == "svg") {
            svg(file)
        } else {
            stop(paste("Unsupported format:", f))
        }

        cat("I: Plotting heatmap to '",getwd(),"/",file,"'.\n",sep="")

        heatmap.2(n.show,
                trace="none", key=FALSE,
                col = colorRampPalette(c("white", "gray50", "black"))(100),
                scale="none", margins=c(10,10),
                dendrogram="none",
                Rowv=TRUE, Colv=TRUE,
                main="Heatmap of Gene Selection for EMT in TP73 Confirmed Contexts",
                labRow=n,
                colRow=1 + n %in% target.gene.selection.DN + 2 * n %in% target.gene.selection.TA + 4 * n %in% target.gene.selection.always,
                labCol=colnames(n.show),
                colCol=1 + colnames(n.show) %in% relevant.cofactors.selection.dn + 
                    2*colnames(n.show) %in% relevant.cofactors.selection.ta + 
                    4*colnames(n.show) %in% relevant.cofactors.selection.always,
                cexRow=0.75, cexCol=0.75
        )
        # Add legend for red label color
        legend("topright", legend=c("DN","TA","DN+TA"), text.col=c(2,3,4), bty="n", cex=0.9)
        dev.off()
    }
}