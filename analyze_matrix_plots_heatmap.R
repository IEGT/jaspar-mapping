
#
# Heatmap for gene selection
#

suppressPackageStartupMessages(require(corrplot))

# Use ggrepel to avoid overlapping labels and dots
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(ggrepel))

#' Visualizes a correlation matrix using a customizable plot.
#'
#' The `my.corrplot` function creates a graphical display of a correlation matrix,
#' allowing users to easily interpret the strength and direction of relationships
#' between variables. The function supports various customization options for
#' colors, labels, and layout to enhance interpretability.
#'
#' @param m A correlation matrix or a matrix to be visualized.
#' @param formats Character string indicating the file formats for the correlation plot
#' @param type Character string specifying the type of plot ("full", "upper", or "lower").
#' @param cofactors.selection.and.order
#' @param order Character string specifying the ordering method for the correlation matrix.
#' @param title Title for the plot.
#' @param relevant.cofactors.selection.ta A vector of relevant transcription factors (TFs) for the variant TA (transcriptional activation).
#' @param relevant.cofactors.selection.dn A vector of relevant transcription factors (TFs) for the variant missing the TA (delta N).
#' @param external.relevant.cofactors.selection.always A vector of external relevant transcription factors (TFs) to always include.
#'
#' @return A correlation plot is displayed. The function is called for its side effect.
#'
#' @seealso \code{\link[corrplot]{corrplot}}
#' @export
my.corrplot <- function(m,file="correlation_matrix_TFBS",formats="none",
                        mar = c(0, 0, 6, 0),   # Increase top margin (third value)
                        is.corr=FALSE,
                        cofactors.selection.and.order=NULL, order="hclust", type="upper",
                        title="Correlation Matrix of TFBS (all promoters, TP73 confirmed)",
                        relevant.cofactors.selection.ta=c(), relevant.cofactors.selection.dn=c(), external.relevant.cofactors.selection.always=c()) {
        m.sd <- apply(m,2,sd)
        # Remove constant columns
        m <- m[,m.sd>0]

        stopifnot(ncol(m) > 0)
        stopifnot(nrow(m) > 0)

        cor_matrix <- NULL
        if (is.corr) {
            cor_matrix <- m
        } else {
            if (is.null(cofactors.selection.and.order)) {
                cor_matrix <- cor(x=m, use="pairwise.complete.obs", method="pearson")
            } else {
                # cofactors.selection.and.order is a vector of TFs that need to be described by the correlation matrix m
                cofactors.selection.and.order <- cofactors.selection.and.order[cofactors.selection.and.order %in% colnames(m)]
                if (length(cofactors.selection.and.order) == 0) {
                    stop("E: No cofactors selected for correlation matrix.")
                } else {
                    cat("I: Selected for ",length(cofactors.selection.and.order)," cofactors for correlation matrix.\n")
                }
                cor_matrix <- cor(x=m[,cofactors.selection.and.order], use="pairwise.complete.obs", method="pearson")
                cor_matrix <- cor_matrix[cofactors.selection.and.order, cofactors.selection.and.order]
            }
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
                stop(paste0("E: Unsupported format: '", f, "'"))
            }
        
            tl.col <- 1 # black
            if ("original" == order) {
                # Use original order of columns
                tl.col <- ifelse("original" == order, 1 + colnames(cor_matrix) %in% relevant.cofactors.selection.dn + 
                    2*colnames(cor_matrix) %in% relevant.cofactors.selection.ta + 
                    4*colnames(cor_matrix) %in% external.relevant.cofactors.selection.always, 1)
            }

            # Increase the bottom margin to make the title fully visible
            ret <- corrplot(cor_matrix, is.corr=TRUE, order=order, method="color", type=type, mar=mar,
                tl.col= tl.col,
                tl.cex=ifelse(nrow(cor_matrix)>180,0.25,0.45), number.cex=0.25, number.digits=1,insig="n",
                addCoef.col=NULL,
                title=title)
            legend("bottomleft", legend=c("enriched DN","enriched TA","enriched DN+TA"), fill=c(1+1,1+2,1+1+2), bty="n", cex=0.9, inset=0.02, xpd=TRUE)

            if ("none" != f) {
                dev.off()
                cat("I: Correlation matrix plot saved to '",file,"'.\n")
            }            
        }
        invisible(ret)
}


if (FALSE) {

    # For manual selection of target genes while debugging

    #target.gene.selection <- promoterBedTables$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.promoter.bed$Gene
    #target.gene.selection <- gene.selection.EMT.TAorDN

    target.gene.selection.TA <- e.ta.up
    target.gene.selection.DN <- e.dn.up
    # or
    target.gene.selection.TA <- gene.selection.EMT.TA
    target.gene.selection.DN <- gene.selection.EMT.DN
    # copying internal data from prior run on upregulated genes for EMT subset
    external.relevant.cofactors.selection.ta=rownames(target.gene.vs.universe.prior.to.selection_TA.sorted)
    external.relevant.cofactors.selection.dn=rownames(target.gene.vs.universe.prior.to.selection_DN.sorted)

    #relevant.cofactors.ta.names <- names(which(relevant.cofactors.ta))
    #relevant.cofactors.dn.names <- names(which(relevant.cofactors.dn))
}


if (0 == length(ls(pattern="all_inPromoter_tp73ConfirmAny"))) {
    cat("I: No data found for all_inPromoter_tp73ConfirmAny, retrieving context data...\n")
    source("analyze_matrix_function_retrieve_context.R")
}

my.heatmap.and.corrplot <- function(base.filename="gene_selection_tp73ConfirmAny_unnamed",
                                    formats=c("pdf","png","svg"),
                                    target.gene.selection.TA,
                                    target.gene.selection.DN,
                                    target.gene.selection.TAandDN=c("TP73"),
                                    expression.data=combined.expression.data,
                                    external.relevant.cofactors.selection.ta=NULL,
                                    external.relevant.cofactors.selection.dn=NULL,
                                    external.relevant.cofactors.selection.always=c("TP73 (MA0861.1)"),
                                    universe.prior.to.selection.TA = all_inPromoter_tp73ConfirmTA,
                                    universe.prior.to.selection.DN = all_inPromoter_tp73ConfirmDN,
                                    universe.prior.to.selection.TAandDN = all_inPromoter_tp73ConfirmAll,
                                    universe.prior.to.selection.TAorDN = all_inPromoter_tp73ConfirmAny,
                                    threshold.min.frequency = 0.05,
                                    threshold.min.log2.enrichment = 0.1,
                                    threshold.min.log2.depletion = -threshold.min.log2.enrichment,
                                    filter.for.human.jaspar=TRUE,
                                    filter.for.strong.enrichment.num=NA,
                                    show.cutandrun=TRUE,
                                    show.score=TRUE,
                                    plot.corrplot=TRUE,
                                    plot.heatmap=TRUE,
                                    title.heatmap=NULL,
                                    verbose=TRUE) {

    stopifnot(all(formats %in% c("pdf","png","svg","none")))
    stopifnot(!is.null(target.gene.selection.TA) && is.vector(target.gene.selection.TA))
    stopifnot(!is.null(target.gene.selection.DN) && is.vector(target.gene.selection.DN))
    # The intersection of TA and DN may be empty
    stopifnot(is.null(target.gene.selection.TAandDN) || is.vector(target.gene.selection.TAandDN))

    stopifnot(!is.null(expression.data))
    stopifnot("Gene Symbol" %in% colnames(expression.data))

    stopifnot(!is.null(universe.prior.to.selection.TA))
    stopifnot(!is.null(universe.prior.to.selection.DN))
    stopifnot(!is.null(universe.prior.to.selection.TAandDN))
    stopifnot(!is.null(universe.prior.to.selection.TAorDN))

    # Subset combined expression data for target genes - array of >77500 logicals, indicating rows of expression data
    expression.data.selected.TA = expression.data[,"Gene Symbol"] %in% target.gene.selection.TA
    expression.data.selected.DN = expression.data[,"Gene Symbol"] %in% target.gene.selection.DN
    # for joined representation in Heatmap
    expression.data.selected.TAandDN = expression.data[,"Gene Symbol"] %in% target.gene.selection.TAandDN |
                                         (expression.data[,"Gene Symbol"] %in% target.gene.selection.TA & 
                                          expression.data[,"Gene Symbol"] %in% target.gene.selection.DN)

    expression.data.selected.TAorDN = expression.data[,"Gene Symbol"] %in% c(target.gene.selection.TA,target.gene.selection.DN,target.gene.selection.TAandDN)


    # Context data retrieval for genes  - presented as list of properties
    # > names(target.gene.selection_tp73ConfirmTA)
    # [1] "context_data"                  "context_data_binary"           "context_shifts"
    # [4] "context_matches"               "downstream_genes"              "cutandrun_data"
    # [7] "coordinates"                   "col_sums_by_chromosome"        "col_sums_total"
    # [10] "col_sums_by_chromosome_binary" "col_sums_binary_total"         "num_matches_by_chromosome"
    # [13] "num_matches_total"             "mean_by_chromosome"            "mean_total"
    # [16] "mean_total_binary"
    cat("I: Retrieving context data for target genes in TP73 confirmed contexts...takes a while.\n")
    target.gene.selection_tp73ConfirmTA <- retrieve_context_data_by_chromosome(expression.data.selected.TA, confirmation=c("tp73"),TA.or.DN="TA", verbose=verbose)
    target.gene.selection_tp73ConfirmDN <- retrieve_context_data_by_chromosome(expression.data.selected.DN, confirmation=c("tp73"),TA.or.DN="DN", verbose=verbose)
    target.gene.selection_tp73ConfirmTAandDN <- retrieve_context_data_by_chromosome(expression.data.selected.TAandDN, confirmation=c("tp73"),TA.or.DN="all", verbose=verbose)
    target.gene.selection_tp73ConfirmTAorDN <- retrieve_context_data_by_chromosome(expression.data.selected.TAorDN, confirmation=c("tp73"),TA.or.DN="any", verbose=verbose)

    if (verbose) {
        cat("I: Distribution of number of matches by chromosome:\n")
        cat("TA:\n"); print(target.gene.selection_tp73ConfirmTA$num_matches_by_chromosome)
        cat("DN:\n"); print(target.gene.selection_tp73ConfirmDN$num_matches_by_chromosome)
        cat("TA&DN:\n"); print(target.gene.selection_tp73ConfirmTAandDN$num_matches_by_chromosome)
        cat("TA|DN:\n"); print(target.gene.selection_tp73ConfirmTAorDN$num_matches_by_chromosome)

        cat("I: Total number of matches across all chromosomes:\n",
            "   TA:",sum(target.gene.selection_tp73ConfirmTA$num_matches_by_chromosome),"\n",
            "   DN:",sum(target.gene.selection_tp73ConfirmDN$num_matches_by_chromosome),"\n",
            "   TA&DN:",sum(target.gene.selection_tp73ConfirmTAandDN$num_matches_by_chromosome),"\n",
            "   TA|DN:",sum(target.gene.selection_tp73ConfirmTAorDN$num_matches_by_chromosome),"\n",
            sep="")
    }

    # Check if the context data is available
    cat("I: Overview on the number of observations made for first 10 TFs in neighbourhood of confirmed p73 binding sites.\n")
    target.gene.selection_tp73ConfirmTA.colSums <- colSums(target.gene.selection_tp73ConfirmTA$context_data_binary,na.rm=T)
    if (verbose) {cat("TA:\n"); print(target.gene.selection_tp73ConfirmTA.colSums[1:10])}
    target.gene.selection_tp73ConfirmDN.colSums <- colSums(target.gene.selection_tp73ConfirmDN$context_data_binary,na.rm=T)
    if (verbose) {cat("DN:\n"); print(target.gene.selection_tp73ConfirmDN.colSums[1:10])}
    target.gene.selection_tp73ConfirmTAandDN.colSums <- colSums(target.gene.selection_tp73ConfirmTAandDN$context_data_binary,na.rm=T)
    if (verbose) {cat("TA&DN:\n"); print(target.gene.selection_tp73ConfirmTAandDN.colSums[1:10])}
    target.gene.selection_tp73ConfirmTAorDN.colSums <- colSums(target.gene.selection_tp73ConfirmTAorDN$context_data_binary,na.rm=T)
    if (verbose) {cat("TA|DN:\n"); print(target.gene.selection_tp73ConfirmTAorDN.colSums[1:10])}


    # In analogy to Volcano plots - determine enrichment in comparison to background of arbitrary p73 binding sites in promoter regions
    target.gene.vs.universe.prior.to.selection_TA <- cbind(
        mean.enriched=target.gene.selection_tp73ConfirmTA$mean_total_binary,
        mean.background=universe.prior.to.selection.TA$mean_total_binary,
        log2.ratio=log2(target.gene.selection_tp73ConfirmTA$mean_total_binary / universe.prior.to.selection.TA$mean_total_binary)
    )

    target.gene.vs.universe.prior.to.selection_DN <- cbind(
        mean.enriched=target.gene.selection_tp73ConfirmDN$mean_total_binary,
        mean.background=universe.prior.to.selection.DN$mean_total_binary,
        log2.ratio=log2(target.gene.selection_tp73ConfirmDN$mean_total_binary / universe.prior.to.selection.DN$mean_total_binary)
    )
    
     target.gene.vs.universe.prior.to.selection_TAandDN <- cbind(
        mean.enriched=target.gene.selection_tp73ConfirmTAandDN$mean_total_binary,
        mean.background=universe.prior.to.selection.TAandDN$mean_total_binary,
        log2.ratio=log2(target.gene.selection_tp73ConfirmTAandDN$mean_total_binary / universe.prior.to.selection.TAandDN$mean_total_binary)
    )

    target.gene.vs.universe.prior.to.selection_TAorDN <- cbind(
        mean.enriched=target.gene.selection_tp73ConfirmTAorDN$mean_total_binary,
        mean.background=universe.prior.to.selection.TAorDN$mean_total_binary,
        log2.ratio=log2(target.gene.selection_tp73ConfirmTAorDN$mean_total_binary / universe.prior.to.selection.TAorDN$mean_total_binary)
    )

    # Ranking cofactors for their ratio (== enrichment)
    target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA[order(target.gene.vs.universe.prior.to.selection_TA[, "log2.ratio"],decreasing=T), ]
    head(target.gene.vs.universe.prior.to.selection_TA.sorted,5)
    if (verbose) cat("I: TA Number of target genes vs universe prior to selection: ", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted), "\n")
    rownames(target.gene.vs.universe.prior.to.selection_TA.sorted) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_TA.sorted))

    target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN[order(target.gene.vs.universe.prior.to.selection_DN[, "log2.ratio"],decreasing=T), ]
    head(target.gene.vs.universe.prior.to.selection_DN.sorted,5)
    if (verbose) cat("I: DN Number of target genes vs universe prior to selection: ", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted), "\n")
    rownames(target.gene.vs.universe.prior.to.selection_DN.sorted) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_DN.sorted))

    target.gene.vs.universe.prior.to.selection_TAandDN.sorted <- target.gene.vs.universe.prior.to.selection_TAandDN[order(target.gene.vs.universe.prior.to.selection_TAandDN[, "log2.ratio"],decreasing=T), ]
    head(target.gene.vs.universe.prior.to.selection_TAandDN.sorted,5)
    if (verbose) cat("I: TA&DN Number of target genes vs universe prior to selection: ", nrow(target.gene.vs.universe.prior.to.selection_TAandDN.sorted), "\n")
    rownames(target.gene.vs.universe.prior.to.selection_TAandDN.sorted) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_TAandDN.sorted))

    target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN[order(target.gene.vs.universe.prior.to.selection_TAorDN[, "log2.ratio"],decreasing=T), ]
    head(target.gene.vs.universe.prior.to.selection_TAorDN.sorted,5)
    if (verbose) cat("I: TA|DN Number of target genes vs universe prior to selection: ", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
    rownames(target.gene.vs.universe.prior.to.selection_TAorDN.sorted) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_TAorDN.sorted))

    # Focus on human JASPAR IDs
    if (filter.for.human.jaspar) {
        if (verbose) cat("I: Filtering for human JASPAR IDs...\n")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection_TA.sorted)),]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection_DN.sorted)),]
        target.gene.vs.universe.prior.to.selection_TAandDN.sorted <- target.gene.vs.universe.prior.to.selection_TAandDN.sorted[is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection_TAandDN.sorted)),]
        target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN.sorted[is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection_TAorDN.sorted)),]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),
                ", TA&DN:", nrow(target.gene.vs.universe.prior.to.selection_TAandDN.sorted),
                ", TA|DN:", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
        }
    }

    # In analogy to Volcano plots - filter out those with a too low enrichment
    if (!is.na(threshold.min.log2.enrichment) && threshold.min.log2.enrichment > 0) {
        if (verbose) cat("I: Filtering for log2 ratio > ",threshold.min.log2.enrichment," if above 0\n",sep="")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"log2.ratio"] < 0 |
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"log2.ratio"] > threshold.min.log2.enrichment,]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"log2.ratio"] < 0 |
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"log2.ratio"] > threshold.min.log2.enrichment,]
        target.gene.vs.universe.prior.to.selection_TAandDN.sorted <- target.gene.vs.universe.prior.to.selection_TAandDN.sorted[
            target.gene.vs.universe.prior.to.selection_TAandDN.sorted[,"log2.ratio"] < 0 |
            target.gene.vs.universe.prior.to.selection_TAandDN.sorted[,"log2.ratio"] > threshold.min.log2.enrichment,]
        target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN.sorted[
            target.gene.vs.universe.prior.to.selection_TAorDN.sorted[,"log2.ratio"] < 0 |
            target.gene.vs.universe.prior.to.selection_TAorDN.sorted[,"log2.ratio"] > threshold.min.log2.enrichment,]

        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),
                ", TA&DN:", nrow(target.gene.vs.universe.prior.to.selection_TAandDN.sorted),
                ", TA|DN:", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
        }
    }    

    # In analogy to Volcano plots - filter out those with a too low depletion
    if (!is.na(threshold.min.log2.depletion) && threshold.min.log2.depletion < 0) {
        if (verbose) cat("I: Filtering for log2 ratio < ",threshold.min.log2.depletion," if below 0\n",sep="")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"log2.ratio"] > 0 |
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"log2.ratio"] < threshold.min.log2.depletion,]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"log2.ratio"] > 0 |
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"log2.ratio"] < threshold.min.log2.depletion,]
        target.gene.vs.universe.prior.to.selection_TAandDN.sorted <- target.gene.vs.universe.prior.to.selection_TAandDN.sorted[
            target.gene.vs.universe.prior.to.selection_TAandDN.sorted[,"log2.ratio"] > 0 |
            target.gene.vs.universe.prior.to.selection_TAandDN.sorted[,"log2.ratio"] < threshold.min.log2.depletion,]
        target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN.sorted[
            target.gene.vs.universe.prior.to.selection_TAorDN.sorted[,"log2.ratio"] > 0 |
            target.gene.vs.universe.prior.to.selection_TAorDN.sorted[,"log2.ratio"] < threshold.min.log2.depletion,]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),
                ", TA&DN:", nrow(target.gene.vs.universe.prior.to.selection_TAandDN.sorted),
                ", TA|DN:", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
        }
    }

    if(!is.na(threshold.min.frequency) && threshold.min.frequency > 0) {
        if (verbose) cat("I: Filtering for mean background > ",threshold.min.frequency,".\n",sep="")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"mean.enriched"] > threshold.min.frequency,]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"mean.enriched"] > threshold.min.frequency,]
        target.gene.vs.universe.prior.to.selection_TAandDN.sorted <- target.gene.vs.universe.prior.to.selection_TAandDN.sorted[
            target.gene.vs.universe.prior.to.selection_TAandDN.sorted[,"mean.enriched"] > threshold.min.frequency,]
        target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN.sorted[
            target.gene.vs.universe.prior.to.selection_TAorDN.sorted[,"mean.enriched"] > threshold.min.frequency,]

        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),
                ", TA&DN:", nrow(target.gene.vs.universe.prior.to.selection_TAandDN.sorted),
                ", TA|DN:", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
        }
    }

    #grep(x=rownames(target.gene.vs.universe.prior.to.selection_TA.sorted),pattern="PAT")
    if (!is.na(filter.for.strong.enrichment.num) && filter.for.strong.enrichment.num > 0) {
        if (verbose) cat("I: Filtering for strong enrichment > ",filter.for.strong.enrichment.num,".\n",sep="")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[1:filter.for.strong.enrichment.num,]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[1:filter.for.strong.enrichment.num,]
        target.gene.vs.universe.prior.to.selection_TAandDN.sorted <- target.gene.vs.universe.prior.to.selection_TAandDN.sorted[1:filter.for.strong.enrichment.num,]
        target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN.sorted[1:filter.for.strong.enrichment.num,]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),
                ", TA&DN:", nrow(target.gene.vs.universe.prior.to.selection_TAandDN.sorted),
                ", TA|DN:", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
        }
    }

    cofactor.names.TA <- unique(c(external.relevant.cofactors.selection.ta,external.relevant.cofactors.selection.always,
                               rownames(target.gene.vs.universe.prior.to.selection_TA.sorted)))
    if (any(is.na(cofactor.names.TA))) {
        cat("W: NA values found in cofactor.names.TA - removing those\n")
        cofactor.names.TA <- cofactor.names.TA[!is.na(cofactor.names.TA)]
    }
    cofactor.names.DN <- unique(c(external.relevant.cofactors.selection.dn,external.relevant.cofactors.selection.always,
                               rownames(target.gene.vs.universe.prior.to.selection_DN.sorted)))
    if (any(is.na(cofactor.names.DN))) {
        cat("W: NA values found in cofactor.names.DN - removing those\n")
        cofactor.names.DN <- cofactor.names.DN[!is.na(cofactor.names.DN)]
    }
    cofactor.names.TAandDN <- unique(c(external.relevant.cofactors.selection.always,
                               rownames(target.gene.vs.universe.prior.to.selection_TAandDN.sorted)))
    if (any(is.na(cofactor.names.TAandDN))) {
        cat("W: NA values found in cofactor.names.TAandDN - removing those\n")
        cofactor.names.TAandDN <- cofactor.names.TAandDN[!is.na(cofactor.names.TAandDN)]
    }
    cofactor.names.TAorDN <- unique(c(external.relevant.cofactors.selection.always,
                               rownames(target.gene.vs.universe.prior.to.selection_TAorDN.sorted)))
    if (any(is.na(cofactor.names.TAorDN))) {
        cat("W: NA values found in cofactor.names.TAorDN - removing those\n")
        cofactor.names.TAorDN <- cofactor.names.TAorDN[!is.na(cofactor.names.TAorDN)]
    }
    # The order matters here - for the looks only
    cofactor.names.union <- unique(c(cofactor.names.TA, cofactor.names.TAandDN, cofactor.names.TAorDN, cofactor.names.DN))
    cofactor.names.union <- unique(c(cofactor.names.TA, cofactor.names.DN))
    
    if (any(is.na(cofactor.names.union))) {
        cat("E: NA values found in cofactor.names.union - removing those, but meant to already have done so in inputs.\n")
        cofactor.names.union <- cofactor.names.union[!is.na(cofactor.names.union)]
    }

    cat("I: # TF of TA in DN: "); print(sum(cofactor.names.TA %in% cofactor.names.DN))
    cat("I: # TF of DN in TA: "); print(sum(cofactor.names.DN %in% cofactor.names.TA))
    cat("I: # TF of Union in TA: "); print(sum(cofactor.names.union %in% cofactor.names.TA))
    cat("I: # TF of Union in DN: "); print(sum(cofactor.names.union %in% cofactor.names.DN))
    cat("I: # TF of Union in TA&DN: "); print(sum(cofactor.names.union %in% cofactor.names.TAandDN))
    cat("I: # TF of Union in TA|DN: "); print(sum(cofactor.names.union %in% cofactor.names.TAorDN))

     # binary representation of context in p73 binding sites
    target.gene.selection_tp73ConfirmTA_context <- target.gene.selection_tp73ConfirmTA$context_data_binary
    # human-readable column names
    colnames(target.gene.selection_tp73ConfirmTA_context) <- prettyIdentifierJaspar(colnames(target.gene.selection_tp73ConfirmTA_context))
    # all the transcription factors in either of the sets should be present after name conversion
    stopifnot(all(cofactor.names.union %in% colnames(target.gene.selection_tp73ConfirmTA_context)))
    # choose the same set of TFs for TA and DN confirmed sites, so we can compare
    target.gene.selection_tp73ConfirmTA_context_interest <- target.gene.selection_tp73ConfirmTA_context[, cofactor.names.union, drop = FALSE]
    # number of observations of TFs within vicinity of confirmed binding sites
    target.gene.selection_tp73ConfirmTA_context_interest.colsums <- colSums(target.gene.selection_tp73ConfirmTA_context_interest)
    # Reduce to those with > 0 observations (important for correlation plot, does not look good in heatmap)
    target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames <- colnames(target.gene.selection_tp73ConfirmTA_context_interest)[target.gene.selection_tp73ConfirmTA_context_interest.colsums > 0]
    # Perform reduction
    target.gene.selection_tp73ConfirmTA_context_interest.informative <- target.gene.selection_tp73ConfirmTA_context_interest[,target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames, drop = FALSE]

    # analogously for DN
    target.gene.selection_tp73ConfirmDN_context <- target.gene.selection_tp73ConfirmDN$context_data_binary
    colnames(target.gene.selection_tp73ConfirmDN_context) <- prettyIdentifierJaspar(colnames(target.gene.selection_tp73ConfirmDN_context))
    # all the transcription factors in either of the sets should be present after name conversion
    stopifnot(all(cofactor.names.union %in% colnames(target.gene.selection_tp73ConfirmDN_context)))
    target.gene.selection_tp73ConfirmDN_context_interest <- target.gene.selection_tp73ConfirmDN_context[, cofactor.names.union, drop = FALSE]
    target.gene.selection_tp73ConfirmDN_context_interest.colsums <- colSums(target.gene.selection_tp73ConfirmDN_context_interest)
    target.gene.selection_tp73ConfirmDN_context_interest.informativeColumNames <- colnames(target.gene.selection_tp73ConfirmDN_context_interest)[target.gene.selection_tp73ConfirmDN_context_interest.colsums > 0]
    target.gene.selection_tp73ConfirmDN_context_interest.informative <- target.gene.selection_tp73ConfirmDN_context_interest[,target.gene.selection_tp73ConfirmDN_context_interest.informativeColumNames, drop = FALSE]

    # analogously for TAandDN
    target.gene.selection_tp73ConfirmTAandDN_context <- target.gene.selection_tp73ConfirmTAandDN$context_data_binary
    colnames(target.gene.selection_tp73ConfirmTAandDN_context) <- prettyIdentifierJaspar(colnames(target.gene.selection_tp73ConfirmTAandDN_context))
    # all the transcription factors in either of the sets should be present after name conversion
    stopifnot(all(cofactor.names.union %in% colnames(target.gene.selection_tp73ConfirmTAandDN_context)))
    target.gene.selection_tp73ConfirmTAandDN_context_interest <- target.gene.selection_tp73ConfirmTAandDN_context[, cofactor.names.union, drop = FALSE]
    target.gene.selection_tp73ConfirmTAandDN_context_interest.colsums <- colSums(target.gene.selection_tp73ConfirmTAandDN_context_interest)
    target.gene.selection_tp73ConfirmTAandDN_context_interest.informativeColumNames <- colnames(target.gene.selection_tp73ConfirmTAandDN_context_interest)[target.gene.selection_tp73ConfirmTAandDN_context_interest.colsums > 0]
    target.gene.selection_tp73ConfirmTAandDN_context_interest.informative <- target.gene.selection_tp73ConfirmTAandDN_context_interest[,target.gene.selection_tp73ConfirmTAandDN_context_interest.informativeColumNames, drop = FALSE]

    # Preparation for heatmap and correlation matrix
    stopifnot(all(colnames(target.gene.selection_tp73ConfirmTA_context_interest) == colnames(target.gene.selection_tp73ConfirmDN_context_interest)))

    # Deriving informative column names for TAorDN - which is not explicitly passed as an argument to the function
    target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames.TAorDN <- colnames(target.gene.selection_tp73ConfirmTA_context_interest)[target.gene.selection_tp73ConfirmTA_context_interest.colsums > 0 | target.gene.selection_tp73ConfirmDN_context_interest.colsums > 0 | target.gene.selection_tp73ConfirmTAandDN_context_interest.colsums > 0]
    target.gene.selection_tp73ConfirmDN_context_interest.informativeColumNames.TAorDN <- colnames(target.gene.selection_tp73ConfirmDN_context_interest)[target.gene.selection_tp73ConfirmTA_context_interest.colsums > 0 | target.gene.selection_tp73ConfirmDN_context_interest.colsums > 0 | target.gene.selection_tp73ConfirmTAandDN_context_interest.colsums > 0]
    stopifnot(all(target.gene.selection_tp73ConfirmDN_context_interest.informativeColumNames.TAorDN ==target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames.TAorDN ))

    # Start analogously to the TA and DN sets
    target.gene.selection_tp73ConfirmTAorDN_context <- target.gene.selection_tp73ConfirmTAorDN$context_data_binary
    colnames(target.gene.selection_tp73ConfirmTAorDN_context) <- prettyIdentifierJaspar(colnames(target.gene.selection_tp73ConfirmTAorDN_context))
    # all the transcription factors in either of the sets should be present after name conversion
    stopifnot(all(cofactor.names.union %in% colnames(target.gene.selection_tp73ConfirmTAorDN_context)))
    target.gene.selection_tp73ConfirmTAorDN_context_interest <- target.gene.selection_tp73ConfirmTAorDN_context[, cofactor.names.union, drop = FALSE]
    # Do not compute colums, just take the same TFs that have been determined for TA and DN individually
    target.gene.selection_tp73ConfirmTAorDN_context_interest.informative <- target.gene.selection_tp73ConfirmTAorDN_context_interest[,target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames.TAorDN, drop = FALSE]       

    if (plot.corrplot) {

        m.cor.ta <- m.cor.dn <- NULL

        # Plot correlation matrix of all transcription factors (columns) against each other
        ret.corrplot.ta <- my.corrplot(m=target.gene.selection_tp73ConfirmTA_context_interest.informative,
                    relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                    title="Correlation Matrix of TA-confirmed p73 BS cofactors in genes upregulated up TA overexpression")
        m.cor.ta <- ret.corrplot.ta$corr
        ret.corrplot.dn <- my.corrplot(m=target.gene.selection_tp73ConfirmDN_context_interest.informative,
                    relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                    title="Correlation Matrix of DN-confirmed p73 BS cofactors in gene upregulated up DN overexpression")
        m.cor.dn <- ret.corrplot.dn$corr

        cofactors.union <- unique(c(cofactor.names.TA, cofactor.names.DN, rownames(target.gene.vs.universe.prior.to.selection_TAorDN.sorted)))

        cofactors.union.order.ta <- rownames(ret.corrplot.ta$corr)[rownames(ret.corrplot.ta$corr) %in% cofactors.union]

        # control - should look like the first TA one
        ret.corrplot.ta.ordered.ta <- my.corrplot(m=target.gene.selection_tp73ConfirmTA_context_interest.informative[,cofactors.union.order.ta,drop=F],
                        order="original",
                        cofactors.selection.and.order=cofactors.union.order.ta,
                        relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                        title="Control plot - should look like first - Correlation Matrix of TA")
 
        # now the DN one, but ordered to match the TA one
        ret.corrplot.dn.ordered.ta <- my.corrplot(m=target.gene.selection_tp73ConfirmDN_context_interest.informative[,cofactors.union.order.ta,drop=F],
                        order="original",
                        cofactors.selection.and.order=cofactors.union.order.ta,
                        relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                        title="Correlation Matrix of DN-confirmed p73 BS cofactors in gene upregulated up DN overexpression - order adjusted to match TA")
        # Now we have two correlation matrices, one for TA and one for DN, both ordered by the same set of cofactors, in the same order - explicitly
        m.cor.ta.minus.dn <- ret.corrplot.ta.ordered.ta$corr[cofactors.union.order.ta,cofactors.union.order.ta]-ret.corrplot.dn.ordered.ta$corr[cofactors.union.order.ta,cofactors.union.order.ta]

        # Transform m.cor.ta.minus.dn into a data frame of pairs and values
        m.cor.ta.minus.dn.df <- as.data.frame(as.table(m.cor.ta.minus.dn))
        m.cor.ta.df <- as.data.frame(as.table(ret.corrplot.ta.ordered.ta$corr[cofactors.union.order.ta,cofactors.union.order.ta]))
        m.cor.dn.df <- as.data.frame(as.table(ret.corrplot.dn.ordered.ta$corr[cofactors.union.order.ta,cofactors.union.order.ta]))
        colnames(m.cor.ta.minus.dn.df) <- c("TF1", "TF2", "Difference")
        colnames(m.cor.ta.df) <- c("TF1", "TF2", "TA_Correlation")
        colnames(m.cor.dn.df) <- c("TF1", "TF2", "DN_Correlation")
        stopifnot(all(m.cor.dn.df[,1] == m.cor.ta.df[,1])) # should be TRUE
        stopifnot(all(m.cor.dn.df[,2] == m.cor.ta.df[,2])) # should be TRUE
        stopifnot(all(m.cor.ta.minus.dn.df[,1] == m.cor.ta.df[,1])) # should be TRUE
        stopifnot(all(m.cor.ta.minus.dn.df[,2] == m.cor.ta.df[,2])) # should be TRUE
        # attach the correlation values to the difference data frame
        m.cor.ta.minus.dn.df <- cbind(m.cor.ta.minus.dn.df,correlation.TA=m.cor.ta.df[,3], correlation.DN=m.cor.dn.df[,3])
        # Optionally, filter to upper triangle (excluding diagonal) to avoid duplicates
        m.cor.ta.minus.dn.df <- m.cor.ta.minus.dn.df[as.numeric(m.cor.ta.minus.dn.df$TF1) < as.numeric(m.cor.ta.minus.dn.df$TF2), ]
        print(head(m.cor.ta.minus.dn.df))
        write.table(m.cor.ta.minus.dn.df, file="correlation_difference_TA_minus_DN.tsv", sep="\t", row.names=FALSE, quote=FALSE, dec=",",col.names=TRUE)

        #quantile(m.cor.ta.minus.dn, na.rm=TRUE)
        ret <- my.corrplot(m=m.cor.ta.minus.dn/max(abs(m.cor.ta.minus.dn)), is.corr=TRUE, order="original",
                #addCoef.col=NULL,
                type="lower",
                #tl.cex=ifelse(nrow(m.cor.ta.minus.dn)>180,0.25,0.5), number.cex=0.25, number.digits=1,insig="n",
                #cofactors.selection.and.order=cofactors.union.order.ta,
                relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                title=paste0("Difference of correlation matrix of TA-confirmed p73 BS cofactors in genes upregulated up TA overexpression - DN / scaled by ",max(abs(m.cor.ta.minus.dn)))
            )

        x_vals <- colSums(target.gene.selection_tp73ConfirmTA_context_interest.informative[,cofactors.union]) / nrow(target.gene.selection_tp73ConfirmTA_context_interest.informative)
        y_vals <- colSums(target.gene.selection_tp73ConfirmDN_context_interest.informative[,cofactors.union]) / nrow(target.gene.selection_tp73ConfirmDN_context_interest.informative)
    
        df <- data.frame(
            x = x_vals,
            y = y_vals,
            label = cofactors.union
        )
        p <- ggplot(df, aes(x = x, y = y)) +
            geom_point(color = "blue", size = 2) +
            geom_abline(intercept = 0, slope = 1, color = "black") +
            geom_text_repel(aes(label = label), size = 2.5, color = "darkred", max.overlaps = Inf) +
            labs(
            x = "Relative Abundance in TA",
            y = "Relative Abundance in DN",
            title = "TA and DN-enriched cofactors in contexts of confirmed p73 binding sites upstream of genes upregulated upon isoform overexpression"
            ) +
            theme_bw()
        print(p)

    }
    #debugme <- data.frame(colnames(m),col=1 + colnames(m) %in% relevant.cofactors.selection.dn +  2*colnames(m) %in% relevant.cofactors.selection.ta +  4*colnames(m) %in% relevant.cofactors.selection.always)

    if (plot.heatmap) {

        cat("I: Preparing heatmap...\n")

        # Genes found within 500 bp downstream of TA-confirmed p73 binding sites - list has duplicatess
        target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input.genes <- target.gene.selection_tp73ConfirmTA$downstream_genes
        target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input <- target.gene.selection_tp73ConfirmTA_context_interest.informative * 1 # transform to numeric
        #stopifnot(nrow(target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input) == length(target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input.genes)  )

        target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input.genes <- target.gene.selection_tp73ConfirmDN$downstream_genes
        target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input <- target.gene.selection_tp73ConfirmDN_context_interest.informative * 1 # transform to numeric
        #stopifnot(nrow(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input) == length(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input.genes)  )

        target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input.genes <- target.gene.selection_tp73ConfirmTAandDN$downstream_genes
        target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input <- target.gene.selection_tp73ConfirmTAandDN_context_interest.informative * 1 # transform to numeric
        #stopifnot(nrow(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input) == length(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input.genes)  )

        target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input.genes <- target.gene.selection_tp73ConfirmTAorDN$downstream_genes
        target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input <- target.gene.selection_tp73ConfirmTAorDN_context_interest.informative * 1 # transform to numeric
        #stopifnot(nrow(target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input) == length(target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input.genes)  )

        cat("D: a\n")

        if (show.score) {
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input,
                Score=target.gene.selection_tp73ConfirmTA$coordinates[,c(5)]
            )
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input,
                Score=target.gene.selection_tp73ConfirmDN$coordinates[,c(5)]
            )
            target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input,
                Score=target.gene.selection_tp73ConfirmTAandDN$coordinates[,c(5)]
            )
            target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input,
                Score=target.gene.selection_tp73ConfirmTAorDN$coordinates[,c(5)]
            )
        }


        if (show.cutandrun) {
            
            # Print all possible six combinations of "pos"/"tp73" with "skmel29_2" and "TA"/"DN"/"GFP"
            combinations <- sort(as.vector(outer(c("pos", "tp73"), c("TA", "DN", "GFP"), function(x, y) paste(x, "skmel29_2", y, sep = "_"))))
            cat("D: Adding CUT&RUN age data for p73 binding site\n") ; print(combinations)
            
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input,
                target.gene.selection_tp73ConfirmTA$cutandrun_data[,..combinations] 
            )
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input,
                target.gene.selection_tp73ConfirmDN$cutandrun_data[,..combinations] 
            )
            target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input,
                target.gene.selection_tp73ConfirmTAandDN$cutandrun_data[,..combinations] 
            )
            target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input,
                target.gene.selection_tp73ConfirmTAorDN$cutandrun_data[,..combinations] 
            )
        }

        
        cat("D: c\n")

        # Somewhat unfortunate - removing duplicates so we can uniquely name rows of matrix
        target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input_nonRedundant <-
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input[!duplicated(target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input.genes),]
        rownames(target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input_nonRedundant) <-
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input.genes[!duplicated(target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input.genes)]
        target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input_nonRedundant <-
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input[!duplicated(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input.genes),]
        rownames(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input_nonRedundant) <-
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input.genes[!duplicated(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input.genes)]
        target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input_nonRedundant <-
            target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input[!duplicated(target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input.genes),]
        rownames(target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input_nonRedundant) <-
            target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input.genes[!duplicated(target.gene.selection_tp73ConfirmTAandDN_context_interest_heatmap_input.genes)]
        target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input_nonRedundant <-
            target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input[!duplicated(target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input.genes),]
        rownames(target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input_nonRedundant) <-
            target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input.genes[!duplicated(target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input.genes)]

        cat("D: d\n")

        require(gplots)
        n <- rownames(target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input_nonRedundant)
        n.show <- target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input_nonRedundant
        stopifnot(length(n) == nrow(n.show))

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

        
        cat("D: e\n")

        for (f in formats) {
            file <- paste0(base.filename,".", f)
            if (f == "png") {
                png(file)
            } else if (f == "pdf") {
                if (nrow(n.show) < ncol(n.show)) {
                    # Transposing since more columns than rows
                    # If more rows than columns, adjust size accordingly
                    pdf(file, width=round(2*nrow(n.show)/25+10,0), height=round(3*ncol(n.show)/25,0))
                } else {
                    # If more columns than rows, adjust size accordingly
                    pdf(file, width=round(2*ncol(n.show)/25+10,0), height=round(3*nrow(n.show)/25,0))
                }
            } else if (f == "svg") {
                svg(file)
            } else {
                stop(paste("Unsupported format:", f))
            }

            cat("I: Plotting heatmap to '",getwd(),"/",file,"'.\n",sep="")

            colors <- colorRampPalette(c("white", "gray20"))(100)
            margins <- c(10,10)
            if (nrow(n.show) < ncol(n.show)) {

                cat("D: Transpose if more rows than columns\n")
                heatmap.2(t(n.show),
                    trace="none", key=FALSE,
                    col = colors,
                    scale="none", margins=margins,
                    dendrogram="none",
                    Rowv=TRUE, Colv=TRUE,
                    main=title.heatmap,
                    labRow=colnames(n.show),
                    colRow= 1 + colnames(n.show) %in% cofactor.names.TA + 
                        2*colnames(n.show) %in% cofactor.names.DN,
                    labCol=rownames(n.show),
                    colCol=1 + rownames(n.show) %in% target.gene.selection.TA + 2 * rownames(n.show) %in% target.gene.selection.DN + 4 * rownames(n.show) %in% target.gene.selection.TAandDN,
                    cexRow=ifelse(ncol(n.show)>50,0.55,0.7),
                    cexCol=ifelse(nrow(n.show)>180,ifelse(nrow(n.show)>250,0.35,0.5),0.75)
                )
             
            } else {

                heatmap.2(n.show,
                    trace="none", key=FALSE,
                    col = colors,
                    scale="none", margins=margins,
                    dendrogram="none",
                    Rowv=TRUE, Colv=TRUE,
                    main=title.heatmap,
                    labRow=n,
                    colRow=1 + rownames(n.show) %in% target.gene.selection.TA + 2 * rownames(n.show) %in% target.gene.selection.DN + 4 * rownames(n.show) %in% target.gene.selection.TAandDN,
                    labCol=colnames(n.show),
                    colCol=1 + colnames(n.show) %in% cofactor.names.TA + 
                        2*colnames(n.show) %in% cofactor.names.DN,
                    cexRow=ifelse(nrow(n.show)>50,0.55,0.7),
                    cexCol=ifelse(ncol(n.show)>180,ifelse(ncol(n.show)>250,0.35,0.5),0.75)
                )
            }

            # Add legend for red label color
            legend("topright", legend=c("TA","DN","DN+TA"), text.col=c(1+1+0,1+0+2,1+1+2), bty="n", cex=0.9)
            dev.off()
        }
        cat("I: Finished plotting heatmap for ", title.heatmap, "\n", sep="")
        write.table(rowSums(n.show), file=paste0(base.filename,"_gene_sums_of_cofactors_present.tsv"), sep="\t", quote=FALSE)
        write.table(colSums(n.show), file=paste0(base.filename,"_cofactor_sums_of_genes_targeted.tsv"), sep="\t", quote=FALSE)

    }
    invisible(
        list(
            target.gene.selection.TA=target.gene.selection_tp73ConfirmTA_context_interest.informative,
            target.gene.selection.DN=target.gene.selection_tp73ConfirmDN_context_interest.informative,
            #target.gene.selection.TAandDN=target.gene.selection_tp73ConfirmAll_context_interest.informative,
            #target.gene.selection.TAorDN=target.gene.selection_tp73ConfirmAny_context_interest.informative,
            target.gene.vs.universe.prior.to.selection.TA.sorted=target.gene.vs.universe.prior.to.selection_TA.sorted,
            target.gene.vs.universe.prior.to.selection.DN.sorted=target.gene.vs.universe.prior.to.selection_DN.sorted,
            m.cor.ta=m.cor.ta,
            m.cor.dn=m.cor.dn,
            cofactors.TA=cofactor.names.TA,
            cofactors.DN=cofactor.names.DN,
            cofactors.intersection=intersect(cofactor.names.TA,cofactor.names.DN),
            cofactors.union=cofactors.union
        )
    )
}