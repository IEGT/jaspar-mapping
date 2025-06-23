
#
# Heatmap for gene selection
#

# Retrieves all contexts for binding sites in promoter regions, which is the background for the analysis.
all_inPromoter_tp73ConfirmAny <- retrieve_context_data_by_chromosome(NULL, confirmation=c("tp73","promoter"),TA.or.DN="any")
all_inPromoter_tp73ConfirmAny.colSums <- colSums(all_inPromoter_tp73ConfirmAny$context_data_binary,na.rm=T)


require(corrplot)
# Use ggrepel to avoid overlapping labels and dots
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(ggrepel))

my.corrplot <- function(m,file="correlation_matrix_TFBS",formats="none",
                        cofactors.selection.and.order=NULL, order="hclust", type="upper",
                        title="Correlation Matrix of TFBS (all promoters, TP73 confirmed)",
                        relevant.cofactors.selection.ta=c(), relevant.cofactors.selection.dn=c(), external.relevant.cofactors.selection.always=c()) {
        #m <- target.gene.selection_tp73ConfirmAny_context_interest.informative
        m.sd <- apply(m,2,sd)
        m <- m[,m.sd>0]

        stopifnot(ncol(m) > 0)
        stopifnot(nrow(m) > 0)

        cor_matrix <- NULL
 
        if (is.null(cofactors.selection.and.order)) {
            cor_matrix <- cor(x=m, use="pairwise.complete.obs", method="pearson")
        } else {
            cofactors.selection.and.order <- cofactors.selection.and.order[cofactors.selection.and.order %in% colnames(m)]
            if (length(cofactors.selection.and.order) == 0) {
                stop("E: No cofactors selected for correlation matrix.")
            } else {
                cat("I: Selected for ",length(cofactors.selection.and.order)," cofactors for correlation matrix.\n")
            }
            cor_matrix <- cor(x=m[,cofactors.selection.and.order], use="pairwise.complete.obs", method="pearson")
            cor_matrix <- cor_matrix[cofactors.selection.and.order, cofactors.selection.and.order]
        }

        # Use the context data for all promoters with TP73 confirmation as the matrix

        for(f in formats) {
            
            file <- paste0(file,".",f)
            if (f == "png") {
                png(file, width=2+9*ncol(m)/60, height=1 + 9*ncol(m)/60, res=300)
            } else if (f == "svg") {
                svg(file)
            } else if (f == "pdf") {
                pdf(file=file, width=2+9*ncol(m)/60, height=1 + 9*ncol(m)/60)
            } else if (f == "jpeg") {
                jpeg(file, width=2+9*ncol(m)/60, height=1 + 9*ncol(m)/60, res=300)
            } else if (f == "tiff") {
                tiff(file, width=2+9*ncol(m)/60, height=1 + 9*ncol(m)/60, res=300)
            } else if (f == "none") {
            } else {
                stop(paste("Unsupported format:", f))
            }
        
            # Increase the bottom margin to make the title fully visible
            op <- par(mar = c(0, 0, 6, 0)) # Increase bottom margin (third value)
            corrplot(cor_matrix, is.corr=TRUE, order=order, method="color", type=type,
                tl.col= 1 + colnames(cor_matrix) %in% relevant.cofactors.selection.dn + 
                    2*colnames(cor_matrix) %in% relevant.cofactors.selection.ta + 
                    4*colnames(cor_matrix) %in% external.relevant.cofactors.selection.always,
                tl.cex=0.7, number.cex=0.35, number.digits=1,insig="n",
                addCoef.col=NULL,
                title=title)
            par(op) # Restore previous margins
            #legend("bottomleft", legend=c("enriched DN","enriched TA","enriched DN+TA","newly found in data"), fill=c(1+1,1+2,1+2+4,1), bty="n", cex=0.9, inset=0.02, xpd=TRUE)

            if ("none" != f) {
                dev.off()
                cat("I: Correlation matrix plot saved to '",file,"'.\n")
            }            
        }
        invisible(cor_matrix)
}

require(ggplot2)

#gene.selection <- promoterBedTables$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.promoter.bed$Gene
#gene.selection <- gene.selection.EMT.TAorDN

target.gene.selection.TA <- e.ta.up
target.gene.selection.DN <- e.dn.up

relevant.cofactors.ta.names <- names(which(relevant.cofactors.ta))
relevant.cofactors.dn.names <- names(which(relevant.cofactors.dn))


my.heatmap.and.corrplot <- function(base.filename="gene_selection_tp73ConfirmAny_unnamed",
                                    formats=c("pdf","png","svg"),
                                    target.gene.selection.TA,
                                    target.gene.selection.DN,
                                    target.gene.selection.always=c("TP73"),
                                    combined.expression.data = combined.expression.data,
                                    external.relevant.cofactors.selection.ta=NULL,
                                    external.relevant.cofactors.selection.dn=NULL,
                                    external.relevant.cofactors.selection.always=c("TP73 (MA0861.1)"),
                                    universe.prior.to.selection.TA = all_inPromoter_tp73ConfirmAny,
                                    universe.prior.to.selection.DN = all_inPromoter_tp73ConfirmAny,
                                    threshold.min.frequency = 0.05,
                                    threshold.min.log2.enrichment = 0.1,
                                    threshold.min.log2.depletion = -threshold.min.log2.enrichment,
                                    #threshold.min.log2.depletion = NA,
                                    filter.for.human.jaspar=TRUE,
                                    filter.for.strong.enrichment.num=NA,
                                    show.cutandrun=TRUE,
                                    show.score=TRUE,
                                    plot.corrplot=TRUE,
                                    plot.heatmap=TRUE,
                                    verbose=TRUE) {

    stopifnot(all(formats %in% c("pdf","png","svg","none")))
    stopifnot("Gene Symbol" %in% colnames(combined.expression.data))
    combined.expression.data.selected.TA = combined.expression.data[,"Gene Symbol"] %in% target.gene.selection.TA
    combined.expression.data.selected.DN = combined.expression.data[,"Gene Symbol"] %in% target.gene.selection.DN

    cat("I: Retrieving context data for target genes in TP73 confirmed contexts...takes a while.\n")
    target.gene.selection_tp73ConfirmTA <- retrieve_context_data_by_chromosome(combined.expression.data.selected.TA, confirmation=c("tp73"),TA.or.DN="TA", verbose=verbose)
    target.gene.selection_tp73ConfirmDN <- retrieve_context_data_by_chromosome(combined.expression.data.selected.DN, confirmation=c("tp73"),TA.or.DN="DN", verbose=verbose)

    if (verbose) {
        cat("I: Distribution of number of matches by chromosome:\n")
        cat("TA: "); print(target.gene.selection_tp73ConfirmTA$num_matches_by_chromosome)
        cat("DN: "); print(target.gene.selection_tp73ConfirmDN$num_matches_by_chromosome)
        cat("I: Total number of matches across all chromosomes:\n",
            "   TA:",sum(target.gene.selection_tp73ConfirmTA$num_matches_by_chromosome),"\n",
            "   DN:",sum(target.gene.selection_tp73ConfirmDN$num_matches_by_chromosome),"\n",sep="")
    }

    target.gene.selection_tp73ConfirmTA.colSums <- colSums(target.gene.selection_tp73ConfirmTA$context_data_binary,na.rm=T)
    if (verbose) print(target.gene.selection_tp73ConfirmTA.colSums)[1:10]
    target.gene.selection_tp73ConfirmDN.colSums <- colSums(target.gene.selection_tp73ConfirmDN$context_data_binary,na.rm=T)
    if (verbose) print(target.gene.selection_tp73ConfirmDN.colSums)[1:10]

    target.gene.vs.universe.prior.to.selection_TA <- cbind(
        mean.enriched=target.gene.selection_tp73ConfirmTA$mean_total_binary,
        mean.background=universe.prior.to.selection$mean_total_binary,
        log2.ratio=log2(target.gene.selection_tp73ConfirmTA$mean_total_binary / universe.prior.to.selection$mean_total_binary)
    )

    target.gene.vs.universe.prior.to.selection_DN <- cbind(
        mean.enriched=target.gene.selection_tp73ConfirmDN$mean_total_binary,
        mean.background=universe.prior.to.selection$mean_total_binary,
        log2.ratio=log2(target.gene.selection_tp73ConfirmDN$mean_total_binary / universe.prior.to.selection$mean_total_binary)
    )

    # Ranking genes for their ratio (== enrichment)
    target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection[order(target.gene.vs.universe.prior.to.selection_TA[, "log2.ratio"],decreasing=T), ]
    if (verbose) cat("I: Number of target genes vs universe prior to selection: ", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted), "\n")
    rownames(target.gene.vs.universe.prior.to.selection_TA.sorted) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_TA.sorted))

    target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection[order(target.gene.vs.universe.prior.to.selection_DN[, "log2.ratio"],decreasing=T), ]
    if (verbose) cat("I: Number of target genes vs universe prior to selection: ", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted), "\n")
    rownames(target.gene.vs.universe.prior.to.selection_DN.sorted) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_DN.sorted))

    if (filter.for.human.jaspar) {
        if (verbose) cat("I: Filtering for human JASPAR IDs...\n")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection_TA.sorted)),]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection_DN.sorted)),]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),"\n")
        }
    }

    if (!is.na(threshold.min.log2.enrichment) && threshold.min.log2.enrichment > 0) {
        if (verbose) cat("I: Filtering for log2 ratio > ",threshold.min.log2.enrichment," if above 0\n",sep="")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"log2.ratio"] < 0 |
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"log2.ratio"] > threshold.min.log2.enrichment,]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"log2.ratio"] < 0 |
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"log2.ratio"] > threshold.min.log2.enrichment,]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),"\n")
        }
    }    

    if (!is.na(threshold.min.log2.depletion) && threshold.min.log2.depletion < 0) {
        if (verbose) cat("I: Filtering for log2 ratio < ",threshold.min.log2.depletion," if below 0\n",sep="")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"log2.ratio"] > 0 |
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"log2.ratio"] < threshold.min.log2.depletion,]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"log2.ratio"] > 0 |
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"log2.ratio"] < threshold.min.log2.depletion,]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),"\n")
        }
    }

    if(!is.na(threshold.min.frequency) && threshold.min.frequency > 0) {
        if (verbose) cat("I: Filtering for mean background > ",threshold.min.frequency,".\n",sep="")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"mean.enriched"] > threshold.min.frequency,]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"mean.enriched"] > threshold.min.frequency,]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),"\n")
        }
    }

    if (!is.na(filter.for.strong.enrichment.num) && filter.for.strong.enrichment.num > 0) {
        if (verbose) cat("I: Filtering for strong enrichment > ",filter.for.strong.enrichment.num,".\n",sep="")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[1:filter.for.strong.enrichment.num,]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[1:filter.for.strong.enrichment.num,]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),"\n")
        }
    }

    cofactor.names.TA <- unique(c(external.relevant.cofactors.selection.ta,external.relevant.cofactors.selection.always,
                               rownames(target.gene.vs.universe.prior.to.selection_TA.sorted)))
    cofactor.names.DN <- unique(c(external.relevant.cofactors.selection.dn,external.relevant.cofactors.selection.always,
                               rownames(target.gene.vs.universe.prior.to.selection_DN.sorted)))
    cofactor.names <- unique(c(cofactor.names.TA, cofactor.names.DN))


    # Print all possible six combinations of "pos"/"tp73" with "skmel29_2" and "TA"/"DN"/"GFP"
    combinations <- sort(as.vector(outer(c("pos", "tp73"), c("TA", "DN", "GFP"), function(x, y) paste(x, "skmel29_2", y, sep = "_"))))
    #print(combinations)

    target.gene.selection_tp73ConfirmTA_context <- target.gene.selection_tp73ConfirmTA$context_data_binary
    colnames(target.gene.selection_tp73ConfirmTA_context) <- prettyIdentifierJaspar(colnames(target.gene.selection_tp73ConfirmTA_context))
    target.gene.selection_tp73ConfirmTA_context_interest <- target.gene.selection_tp73ConfirmTA_context[, cofactor.names, drop = FALSE]
    target.gene.selection_tp73ConfirmTA_context_interest.colsums <- colSums(target.gene.selection_tp73ConfirmTA_context_interest)
    target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames <- colnames(target.gene.selection_tp73ConfirmTA_context_interest)[target.gene.selection_tp73ConfirmTA_context_interest.colsums > 0]
    target.gene.selection_tp73ConfirmTA_context_interest.informative <- target.gene.selection_tp73ConfirmTA_context_interest[,target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames, drop = FALSE]

    target.gene.selection_tp73ConfirmDN_context <- target.gene.selection_tp73ConfirmDN$context_data_binary
    colnames(target.gene.selection_tp73ConfirmDN_context) <- prettyIdentifierJaspar(colnames(target.gene.selection_tp73ConfirmDN_context))
    target.gene.selection_tp73ConfirmDN_context_interest <- target.gene.selection_tp73ConfirmDN_context[, cofactor.names, drop = FALSE]
    target.gene.selection_tp73ConfirmDN_context_interest.colsums <- colSums(target.gene.selection_tp73ConfirmDN_context_interest)
    target.gene.selection_tp73ConfirmDN_context_interest.informativeColumNames <- colnames(target.gene.selection_tp73ConfirmDN_context_interest)[target.gene.selection_tp73ConfirmDN_context_interest.colsums > 0]
    target.gene.selection_tp73ConfirmDN_context_interest.informative <- target.gene.selection_tp73ConfirmDN_context_interest[,target.gene.selection_tp73ConfirmDN_context_interest.informativeColumNames, drop = FALSE]

    m.cor.ta <- m.cor.dn <- NULL

    if (plot.corrplot) {
        pdf("correlation_matrix_TA_and_DN_separate.pdf", width=21, height=22)

        # Plot correlation matrix of all transcription factors (columns) against each other
        m.cor.ta <- my.corrplot(m=target.gene.selection_tp73ConfirmTA_context_interest.informative, relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                    title="Correlation Matrix of TA-confirmed p73 BS cofactors in genes upregulated up TA overexpression")
        m.cor.dn <- my.corrplot(m=target.gene.selection_tp73ConfirmDN_context_interest.informative, relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                    title="Correlation Matrix of DN-confirmed p73 BS cofactors in gene upregulated up DN overexpression")
        cofactors.joined <- intersect(colnames(m.cor.ta), colnames(m.cor.dn))
        m.cor.dn.ordered <- my.corrplot(m=target.gene.selection_tp73ConfirmDN_context_interest.informative[,cofactors.joined,drop=F],order="original",
                                cofactors.selection.and.order=cofactors.joined,
                                relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                                title="Correlation Matrix of DN-confirmed p73 BS cofactors in gene upregulated up DN overexpression - order adjusted to match TA")
        m.cor.ta.minus.dn <- m.cor.ta[cofactors.joined,cofactors.joined]-m.cor.dn.ordered[cofactors.joined,cofactors.joined]
        #quantile(m.cor.ta.minus.dn, na.rm=TRUE)
        corrplot(corr=m.cor.ta.minus.dn/max(abs(m.cor.ta.minus.dn)), is.corr=TRUE, order="original",addCoef.col=NULL,type="lower",tl.cex=0.7, number.cex=0.35, number.digits=1,insig="n",
                title=paste0("Difference of correlation matrix of TA-confirmed p73 BS cofactors in genes upregulated up TA overexpression - DN / scaled by ",max(abs(m.cor.ta.minus.dn))))

        x_vals <- colSums(target.gene.selection_tp73ConfirmTA_context_interest.informative[,cofactors.joined]) / nrow(target.gene.selection_tp73ConfirmTA_context_interest.informative)
        y_vals <- colSums(target.gene.selection_tp73ConfirmDN_context_interest.informative[,cofactors.joined]) / nrow(target.gene.selection_tp73ConfirmDN_context_interest.informative)
    
        df <- data.frame(
            x = x_vals,
            y = y_vals,
            label = cofactors.joined
        )
        p <- ggplot(df, aes(x = x, y = y)) +
            geom_point(color = "blue", size = 2) +
            geom_abline(intercept = 0, slope = 1, color = "black") +
            geom_text_repel(aes(label = label), size = 3, color = "darkred", max.overlaps = Inf) +
            labs(
            x = "Relative Abundance in TA",
            y = "Relative Abundance in DN",
            title = "TA and DN-enriched cofactors in contexts of confirmed p73 binding sites upstream of genes upregulated upon isoform overexpression"
            ) +
            theme_bw()
        print(p)

        dev.off()
    }
    #debugme <- data.frame(colnames(m),col=1 + colnames(m) %in% relevant.cofactors.selection.dn +  2*colnames(m) %in% relevant.cofactors.selection.ta +  4*colnames(m) %in% relevant.cofactors.selection.always)

    if (plot.heatmap) {

        target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input <- cbind(
            Gene=target.gene.selection_tp73ConfirmTA$downstream_genes,
            target.gene.selection_tp73ConfirmTA_context_interest.informative
        )
        target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input <- cbind(
            Gene=target.gene.selection_tp73ConfirmDN$downstream_genes,
            target.gene.selection_tp73ConfirmDN_context_interest.informative
        )

        if (show.score) {
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input,
                Score=gene.selection_tp73ConfirmTA$coordinates[,c(5)]
            )
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input,
                Score=gene.selection_tp73ConfirmDN$coordinates[,c(5)]
            )
        }

        if (show.cutandrun) {
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input,
                target.gene.selection_tp73ConfirmTA$cutandrun_data[,..combinations] 
            )
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input,
                target.gene.selection_tp73ConfirmDN$cutandrun_data[,..combinations] 
            )
        }

        # FIXME: Is this avoidable
        target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input_nonRedundant <-
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input[!duplicated(target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input[,"Gene"]),]
        rownames(target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input_nonRedundant) <- target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input_nonRedundant$Gene
        target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input_nonRedundant <-
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input[!duplicated(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input[,"Gene"]),]
        rownames(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input_nonRedundant) <- target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input_nonRedundant$Gene



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
    invisible(
        list(
            target.gene.selection.TA=target.gene.selection_tp73ConfirmTA_context_interest.informative,
            target.gene.selection.DN=target.gene.selection_tp73ConfirmDN_context_interest.informative,
            target.gene.selection.always=target.gene.selection_tp73ConfirmAny_context_interest.informative,
            target.gene.vs.universe.prior.to.selection.TA.sorted=target.gene.vs.universe.prior.to.selection_TA.sorted,
            target.gene.vs.universe.prior.to.selection.DN.sorted=target.gene.vs.universe.prior.to.selection_DN.sorted,
            m.cor.ta=m.cor.ta,
            m.cor.dn=m.cor.dn
        )
    )
}


pdf("correlation_upregulated_TA_DN.pdf", width=12, height=10)
r <- my.heatmap.and.corrplot(target.gene.selection.TA=e.ta.up, target.gene.selection.DN=e.dn.up,formats="none",combined.expression.data = combined.expression.data)
dev.off()
