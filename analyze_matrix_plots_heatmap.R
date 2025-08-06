
#
# Heatmap for gene selection
#

## Figure 4

# Retrieves all contexts for binding sites in promoter regions, which is the background for the analysis.
all_inPromoter_tp73ConfirmAny <- retrieve_context_data_by_chromosome(NULL, confirmation=c("tp73","promoter"),TA.or.DN="any")
all_inPromoter_tp73ConfirmAny.colSums <- colSums(all_inPromoter_tp73ConfirmAny$context_data_binary,na.rm=T)

all_inPromoter_tp73ConfirmTA <- retrieve_context_data_by_chromosome(NULL, confirmation=c("tp73","promoter"),TA.or.DN="TA")
all_inPromoter_tp73ConfirmTA.colSums <- colSums(all_inPromoter_tp73ConfirmTA$context_data_binary,na.rm=T)

all_inPromoter_tp73ConfirmDN <- retrieve_context_data_by_chromosome(NULL, confirmation=c("tp73","promoter"),TA.or.DN="DN")
all_inPromoter_tp73ConfirmDN.colSums <- colSums(all_inPromoter_tp73ConfirmDN$context_data_binary,na.rm=T)

require(ggplot2)
require(ggrepel)

target.gene.selection <- list()
gene.set.names <- c("HALLMARK_ANGIOGENESIS",
"HALLMARK_APOPTOSIS",
"HALLMARK_E2F_TARGETS",
"HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
"HALLMARK_HYPOXIA",
"HALLMARK_IL6_JAK_STAT3_SIGNALING",
"HALLMARK_INFLAMMATORY_RESPONSE",
"HALLMARK_INTERFERON_ALPHA_RESPONSE",
"HALLMARK_INTERFERON_GAMMA_RESPONSE",
"HALLMARK_P53_PATHWAY",
"HALLMARK_PI3K_AKT_MTOR_SIGNALING",
"HALLMARK_TGF_BETA_SIGNALING",
"HALLMARK_WNT_BETA_CATENIN_SIGNALING",
"HP_MELANOMA",
"JIANG_MELANOMA_TRM_HIGH_SURVIVAL_22_GENE_SIGNATURE",
"KAUFFMANN_MELANOMA_RELAPSE_DN",
"KAUFFMANN_MELANOMA_RELAPSE_UP",
"KEGG_APOPTOSIS",
"KEGG_MELANOMA",
"LIN_MELANOMA_COPY_NUMBER_DN",
"LIN_MELANOMA_COPY_NUMBER_UP",
"NABA_MATRISOME_HIGHLY_METASTATIC_MELANOMA",
"NABA_MATRISOME_HIGHLY_METASTATIC_MELANOMA_TUMOR_CELL_DERIVED",
"NABA_MATRISOME_POORLY_METASTATIC_MELANOMA",
"NABA_MATRISOME_POORLY_METASTATIC_MELANOMA_TUMOR_CELL_DERIVED",
"ONKEN_UVEAL_MELANOMA_DN",
"ONKEN_UVEAL_MELANOMA_UP",
"REACTOME_APOPTOSIS_INDUCED_DNA_FRAGMENTATION",
"REACTOME_APOPTOSIS",
"REACTOME_DEFECTIVE_INTRINSIC_PATHWAY_FOR_APOPTOSIS",
"REACTOME_EXTRINSIC_PATHWAY_FOR_APOPTOSIS",
"REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS",
"REACTOME_NEURONAL_SYSTEM",
"REACTOME_REGULATION_OF_APOPTOSIS",
"REACTOME_REGULATION_OF_MITF_M_DEPENDENT_GENES_INVOLVED_IN_APOPTOSIS",
"REACTOME_SUPPRESSION_OF_APOPTOSIS",
"REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_CYCLE_GENES",
"REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES",
"REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_DEATH_RECEPTORS_AND_LIGANDS",
"REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_GENES_INVOLVED_IN_G1_CELL_CYCLE_ARREST",
"REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_GENES_INVOLVED_IN_G2_CELL_CYCLE_ARREST",
"REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_SEVERAL_ADDITIONAL_CELL_DEATH_GENES_WHOSE_SPECIFIC_ROLES_IN_P53_DEPENDENT_APOPTOSIS_REMAIN_UNCERTAIN",
"REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53",
"TANG_SENESCENCE_TP53_TARGETS_DN",
"TANG_SENESCENCE_TP53_TARGETS_UP",
"WINNEPENNINCKX_MELANOMA_METASTASIS_DN",
"WINNEPENNINCKX_MELANOMA_METASTASIS_UP",
"WP_APOPTOSIS_MODULATION_AND_SIGNALING",
"WP_APOPTOSIS_MODULATION_BY_HSP70",
"WP_APOPTOSIS",
"WP_HALLMARK_OF_CANCER_METASTASIS_AND_EPITHELIALTOMESENCHYMAL_TRANSITION",
"WP_HALLMARK_OF_CANCER_NONMUTATIONAL_EPIGENETIC_REPROGRAMMING",
"WP_HALLMARK_OF_CANCER_SUSTAINING_PROLIFERATIVE_SIGNALING",
"WP_INHIBITION_OF_GNAQREGULATED_SIGNALING_IN_UVEAL_MELANOMA",
"WP_MELANOMA",
"WP_MIRNA_REGULATION_OF_P53_PATHWAY_IN_PROSTATE_CANCER",
"WP_P53_TRANSCRIPTIONAL_GENE_NETWORK",
"WP_QUERCETIN_AND_NFKB_AP1_INDUCED_APOPTOSIS",
"WP_TP53_NETWORK")

for (hallmarkID in gene.set.names) {
    cat("I: Retrieving genes for ",hallmarkID,"...\n",sep="")
    hallmarkID.selection <- paste(hallmarkID,".promoter.bed",sep="")
    if (!hallmarkID.selection %in% names(promoterBedTables)) {
        cat("E: ",hallmarkID.selection," not found in promoterBedTables. Skipping...\n",sep="")
        break
    }
    target.gene.selection[[hallmarkID]] <- promoterBedTables[[hallmarkID.selection]]$Gene
    cat("D: Added '",hallmarkID,"' with ",length(target.gene.selection[[hallmarkID]])," genes to target.gene.selection.\n",sep="")
}
print(names(target.gene.selection))

pdf(paste0("enrichment_TAvsDN_for_hallmarks.pdf"), width=25, height=25)

for (hallmarkID in names(target.gene.selection)) {
    if (length(target.gene.selection[[hallmarkID]]) == 0) {
        cat("E: No genes found for ",hallmarkID,". Skipping...\n",sep="")
        break
    }
    cat("I: ",hallmarkID,": ",length(target.gene.selection[[hallmarkID]])," genes selected.\n",sep="")
    combined.expression.data.hallmark <- combined.expression.data[,"Gene Symbol"] %in% target.gene.selection[[hallmarkID]]
    if (sum(combined.expression.data.hallmark) == 0) {
        cat("E: No genes found in combined expression data for ",hallmarkID,". Skipping...\n",sep="")
        break
    }
    cat("D: Found ",sum(combined.expression.data.hallmark)," genes in combined expression data of ",
            length(target.gene.selection[[hallmarkID]]), " for ",hallmarkID,".\n",sep="")

    target.gene.selection.hallmark.tp73ConfirmTA <- retrieve_context_data_by_chromosome(combined.expression.data.hallmark, confirmation=c("tp73"),TA.or.DN="TA", verbose=verbose)
    target.gene.selection.hallmark.tp73ConfirmDN <- retrieve_context_data_by_chromosome(combined.expression.data.hallmark, confirmation=c("tp73"),TA.or.DN="DN", verbose=verbose)

    target.gene.selection.hallmark.tp73ConfirmTA.colSums <- colSums(target.gene.selection.hallmark.tp73ConfirmTA$context_data_binary,na.rm=T)
    target.gene.selection.hallmark.tp73ConfirmDN.colSums <- colSums(target.gene.selection.hallmark.tp73ConfirmDN$context_data_binary,na.rm=T)

    if (verbose) {
        print(target.gene.selection.hallmark.tp73ConfirmTA.colSums[1:10])
        print(target.gene.selection.hallmark.tp73ConfirmDN.colSums[1:10])
    }

    target.gene.vs.universe.prior.to.selection_TA <- cbind(
        mean.enriched=target.gene.selection.hallmark.tp73ConfirmTA$mean_total_binary,
        mean.background=all_inPromoter_tp73ConfirmTA$mean_total_binary,
        log2.ratio=log2(target.gene.selection.hallmark.tp73ConfirmTA$mean_total_binary / all_inPromoter_tp73ConfirmTA$mean_total_binary)
    )
    rownames(target.gene.vs.universe.prior.to.selection_TA) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_TA))

    target.gene.vs.universe.prior.to.selection_DN <- cbind(
        mean.enriched=target.gene.selection.hallmark.tp73ConfirmDN$mean_total_binary,
        mean.background=all_inPromoter_tp73ConfirmDN$mean_total_binary,
        log2.ratio=log2(target.gene.selection.hallmark.tp73ConfirmDN$mean_total_binary / all_inPromoter_tp73ConfirmDN$mean_total_binary)
    )
    rownames(target.gene.vs.universe.prior.to.selection_DN) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_DN))


    target.gene.vs.universe.prior.to.selection.hallmark.eligible <- target.gene.vs.universe.prior.to.selection_TA[,"mean.enriched"]>0 | target.gene.vs.universe.prior.to.selection_DN[,"mean.enriched"]>0
    target.gene.vs.universe.prior.to.selection.hallmark.eligible <- target.gene.vs.universe.prior.to.selection.hallmark.eligible & (
        (is.finite(target.gene.vs.universe.prior.to.selection_TA[,"log2.ratio"]) & abs(target.gene.vs.universe.prior.to.selection_TA[,"log2.ratio"]) > 0.1 ) &
        (is.finite(target.gene.vs.universe.prior.to.selection_DN[,"log2.ratio"]) & abs(target.gene.vs.universe.prior.to.selection_DN[,"log2.ratio"]) > 0.1 )
    )
    target.gene.vs.universe.prior.to.selection.hallmark.eligible <- target.gene.vs.universe.prior.to.selection.hallmark.eligible & (
        is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection_TA))
    )

    write.table(cbind(TA=target.gene.vs.universe.prior.to.selection_TA,DN=target.gene.vs.universe.prior.to.selection_DN),
        file=paste0("enrichment_TAvsDN_for_hallmark_",hallmarkID,".tsv"),
        col.names=NA, row.names=TRUE, sep="\t", quote=FALSE, dec=",")


    # Plot the enrichment for TA vs DN for Hallmark gene selection 

    df <- data.frame(
        log2_enrichment_TA = target.gene.vs.universe.prior.to.selection_TA[,"log2.ratio"][target.gene.vs.universe.prior.to.selection.hallmark.eligible],
        log2_enrichment_DN = target.gene.vs.universe.prior.to.selection_DN[,"log2.ratio"][target.gene.vs.universe.prior.to.selection.hallmark.eligible],
        max_mean_enriched = apply(cbind(target.gene.vs.universe.prior.to.selection_TA[,"mean.enriched"][target.gene.vs.universe.prior.to.selection.hallmark.eligible],
                                target.gene.vs.universe.prior.to.selection_DN[,"mean.enriched"][target.gene.vs.universe.prior.to.selection.hallmark.eligible]),1,max),
        TF = prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_TA)[target.gene.vs.universe.prior.to.selection.hallmark.eligible])
    )

    print(quantile(df$max_mean_enriched, probs = 0:20/20))

    # Color and size by max_mean_enriched (continuous)
    p <- ggplot(df, aes(x = log2_enrichment_TA, y = log2_enrichment_DN, label = TF, 
                color = max_mean_enriched, size = max_mean_enriched)) +
        geom_point() +
        geom_text_repel(size = 3, max.overlaps = Inf) +
        scale_color_gradient(low = "red", high = "blue", name = "Max frequency (%)") +
        scale_size_continuous(range = c(2, 8), name = "Max frequency (%)") +
        labs(
            x = "log2 enrichment TA",
            y = "log2 enrichment DN",
            title = paste0("Enrichment of TFBS in ",hallmarkID," genes (TA vs DN)")
        ) +
        theme_bw()

    print(p)
    ggsave(filename = paste0("enrichment_TAvsDN_for_hallmark_",hallmarkID,".png"), plot = p, width = 25, height = 25, units = "in", dpi = 300)
    ggsave(filename = paste0("enrichment_TAvsDN_for_hallmark_",hallmarkID,".svg"), plot = p, width = 25, height = 25, units = "in")

}
dev.off()




require(corrplot)
# Use ggrepel to avoid overlapping labels and dots
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(ggrepel))

my.corrplot <- function(m,file="correlation_matrix_TFBS",formats="none",
                        mar=c(0, 2, 0, 0), # Increase top margin for title
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
            ret <- corrplot(cor_matrix, is.corr=TRUE, order=order, method="color", type=type, mar=mar,
                tl.col= 1 + colnames(cor_matrix) %in% relevant.cofactors.selection.dn + 
                    2*colnames(cor_matrix) %in% relevant.cofactors.selection.ta + 
                    4*colnames(cor_matrix) %in% external.relevant.cofactors.selection.always,
                tl.cex=ifelse(nrow(cor_matrix)>180,0.5,0.7), number.cex=0.35, number.digits=1,insig="n",
                addCoef.col=NULL,
                title=title)
            par(op) # Restore previous margins
            legend("bottomleft", legend=c("enriched DN","enriched TA","enriched DN+TA"), fill=c(1+1,1+2,1+1+2), bty="n", cex=0.9, inset=0.02, xpd=TRUE)

            if ("none" != f) {
                dev.off()
                cat("I: Correlation matrix plot saved to '",file,"'.\n")
            }            
        }
        invisible(ret)
}

require(ggplot2)

#target.gene.selection <- promoterBedTables$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.promoter.bed$Gene
#target.gene.selection <- gene.selection.EMT.TAorDN

#target.gene.selection.TA <- e.ta.up
#target.gene.selection.DN <- e.dn.up


target.gene.selection.TA <- gene.selection.EMT.TA
target.gene.selection.DN <- gene.selection.EMT.DN
# copying internal data from prior run on upregulated genes for EMT subset
external.relevant.cofactors.selection.ta=rownames(target.gene.vs.universe.prior.to.selection_TA.sorted)
external.relevant.cofactors.selection.dn=rownames(target.gene.vs.universe.prior.to.selection_DN.sorted)

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
    # for joined representation in Heatmap
    combined.expression.data.selected.TAorDN = combined.expression.data[,"Gene Symbol"] %in% c(target.gene.selection.TA,target.gene.selection.DN)


    cat("I: Retrieving context data for target genes in TP73 confirmed contexts...takes a while.\n")
    target.gene.selection_tp73ConfirmTA <- retrieve_context_data_by_chromosome(combined.expression.data.selected.TA, confirmation=c("tp73"),TA.or.DN="TA", verbose=verbose)
    target.gene.selection_tp73ConfirmDN <- retrieve_context_data_by_chromosome(combined.expression.data.selected.DN, confirmation=c("tp73"),TA.or.DN="DN", verbose=verbose)
    target.gene.selection_tp73ConfirmTAorDN <- retrieve_context_data_by_chromosome(combined.expression.data.selected.TAorDN, confirmation=c("tp73"),TA.or.DN="any", verbose=verbose)

    if (verbose) {
        cat("I: Distribution of number of matches by chromosome:\n")
        cat("TA:\n"); print(target.gene.selection_tp73ConfirmTA$num_matches_by_chromosome)
        cat("DN:\n"); print(target.gene.selection_tp73ConfirmDN$num_matches_by_chromosome)
        cat("TA/DN:\n"); print(target.gene.selection_tp73ConfirmTAorDN$num_matches_by_chromosome)

        cat("I: Total number of matches across all chromosomes:\n",
            "   TA:",sum(target.gene.selection_tp73ConfirmTA$num_matches_by_chromosome),"\n",
            "   DN:",sum(target.gene.selection_tp73ConfirmDN$num_matches_by_chromosome),"\n",
            "   TA/DN:",sum(target.gene.selection_tp73ConfirmTAorDN$num_matches_by_chromosome),"\n",
            sep="")
    }

    target.gene.selection_tp73ConfirmTA.colSums <- colSums(target.gene.selection_tp73ConfirmTA$context_data_binary,na.rm=T)
    if (verbose) print(target.gene.selection_tp73ConfirmTA.colSums)[1:10]
    target.gene.selection_tp73ConfirmDN.colSums <- colSums(target.gene.selection_tp73ConfirmDN$context_data_binary,na.rm=T)
    if (verbose) print(target.gene.selection_tp73ConfirmDN.colSums)[1:10]
    target.gene.selection_tp73ConfirmTAorDN.colSums <- colSums(target.gene.selection_tp73ConfirmTAorDN$context_data_binary,na.rm=T)
    if (verbose) print(target.gene.selection_tp73ConfirmTAorDN.colSums)[1:10]


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
    
    target.gene.vs.universe.prior.to.selection_TAorDN <- cbind(
        mean.enriched=target.gene.selection_tp73ConfirmTAorDN$mean_total_binary,
        mean.background=universe.prior.to.selection$mean_total_binary,
        log2.ratio=log2(target.gene.selection_tp73ConfirmTAorDN$mean_total_binary / universe.prior.to.selection$mean_total_binary)
    )

    # Ranking genes for their ratio (== enrichment)
    target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA[order(target.gene.vs.universe.prior.to.selection_TA[, "log2.ratio"],decreasing=T), ]
    head(target.gene.vs.universe.prior.to.selection_TA.sorted,3)
    if (verbose) cat("I: Number of target genes vs universe prior to selection: ", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted), "\n")
    rownames(target.gene.vs.universe.prior.to.selection_TA.sorted) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_TA.sorted))

    target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN[order(target.gene.vs.universe.prior.to.selection_DN[, "log2.ratio"],decreasing=T), ]
    head(target.gene.vs.universe.prior.to.selection_DN.sorted,3)
    if (verbose) cat("I: Number of target genes vs universe prior to selection: ", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted), "\n")
    rownames(target.gene.vs.universe.prior.to.selection_DN.sorted) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_DN.sorted))

    target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN[order(target.gene.vs.universe.prior.to.selection_TAorDN[, "log2.ratio"],decreasing=T), ]
    head(target.gene.vs.universe.prior.to.selection_TAorDN.sorted,3)
    if (verbose) cat("I: Number of target genes vs universe prior to selection: ", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
    rownames(target.gene.vs.universe.prior.to.selection_TAorDN.sorted) <- prettyIdentifierJaspar(rownames(target.gene.vs.universe.prior.to.selection_TAorDN.sorted))

    if (filter.for.human.jaspar) {
        if (verbose) cat("I: Filtering for human JASPAR IDs...\n")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection_TA.sorted)),]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection_DN.sorted)),]
        target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN.sorted[is.human.jaspar.id(rownames(target.gene.vs.universe.prior.to.selection_TAorDN.sorted)),]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),
                ", TA/DN:", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
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
        target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN.sorted[
            target.gene.vs.universe.prior.to.selection_TAorDN.sorted[,"log2.ratio"] < 0 |
            target.gene.vs.universe.prior.to.selection_TAorDN.sorted[,"log2.ratio"] > threshold.min.log2.enrichment,]

        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),
                ", TA/DN:", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
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
        target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN.sorted[
            target.gene.vs.universe.prior.to.selection_TAorDN.sorted[,"log2.ratio"] > 0 |
            target.gene.vs.universe.prior.to.selection_TAorDN.sorted[,"log2.ratio"] < threshold.min.log2.depletion,]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),
                ", TA/DN:", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
        }
    }

    if(!is.na(threshold.min.frequency) && threshold.min.frequency > 0) {
        if (verbose) cat("I: Filtering for mean background > ",threshold.min.frequency,".\n",sep="")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[
            target.gene.vs.universe.prior.to.selection_TA.sorted[,"mean.enriched"] > threshold.min.frequency,]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[
            target.gene.vs.universe.prior.to.selection_DN.sorted[,"mean.enriched"] > threshold.min.frequency,]
        target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN.sorted[
            target.gene.vs.universe.prior.to.selection_TAorDN.sorted[,"mean.enriched"] > threshold.min.frequency,]

        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),
                ", TA/DN:", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
        }
    }

    if (!is.na(filter.for.strong.enrichment.num) && filter.for.strong.enrichment.num > 0) {
        if (verbose) cat("I: Filtering for strong enrichment > ",filter.for.strong.enrichment.num,".\n",sep="")
        target.gene.vs.universe.prior.to.selection_TA.sorted <- target.gene.vs.universe.prior.to.selection_TA.sorted[1:filter.for.strong.enrichment.num,]
        target.gene.vs.universe.prior.to.selection_DN.sorted <- target.gene.vs.universe.prior.to.selection_DN.sorted[1:filter.for.strong.enrichment.num,]
        target.gene.vs.universe.prior.to.selection_TAorDN.sorted <- target.gene.vs.universe.prior.to.selection_TAorDN.sorted[1:filter.for.strong.enrichment.num,]
        if (verbose) {
            cat("I: Number of relevant cofactors remaining:",
                "  TA:", nrow(target.gene.vs.universe.prior.to.selection_TA.sorted),
                ", DN:", nrow(target.gene.vs.universe.prior.to.selection_DN.sorted),
                ", TA/DN:", nrow(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), "\n")
        }
    }

    cofactor.names.TA <- unique(c(external.relevant.cofactors.selection.ta,external.relevant.cofactors.selection.always,
                               rownames(target.gene.vs.universe.prior.to.selection_TA.sorted)))
    cofactor.names.DN <- unique(c(external.relevant.cofactors.selection.dn,external.relevant.cofactors.selection.always,
                               rownames(target.gene.vs.universe.prior.to.selection_DN.sorted)))
    cofactor.names.union <- unique(c(cofactor.names.TA, rownames(target.gene.vs.universe.prior.to.selection_TAorDN.sorted), cofactor.names.DN))

    cat("I: All of TA in DN: "); print(all(cofactor.names.TA %in% cofactor.names.DN))
    cat("I: All of DN in TA: "); print(all(cofactor.names.DN %in% cofactor.names.TA))

    # Print all possible six combinations of "pos"/"tp73" with "skmel29_2" and "TA"/"DN"/"GFP"
    combinations <- sort(as.vector(outer(c("pos", "tp73"), c("TA", "DN", "GFP"), function(x, y) paste(x, "skmel29_2", y, sep = "_"))))
    #print(combinations)

    target.gene.selection_tp73ConfirmTA_context <- target.gene.selection_tp73ConfirmTA$context_data_binary
    colnames(target.gene.selection_tp73ConfirmTA_context) <- prettyIdentifierJaspar(colnames(target.gene.selection_tp73ConfirmTA_context))
    target.gene.selection_tp73ConfirmTA_context_interest <- target.gene.selection_tp73ConfirmTA_context[, cofactor.names.union, drop = FALSE]
    target.gene.selection_tp73ConfirmTA_context_interest.colsums <- colSums(target.gene.selection_tp73ConfirmTA_context_interest)
    target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames <- colnames(target.gene.selection_tp73ConfirmTA_context_interest)[target.gene.selection_tp73ConfirmTA_context_interest.colsums > 0]
    target.gene.selection_tp73ConfirmTA_context_interest.informative <- target.gene.selection_tp73ConfirmTA_context_interest[,target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames, drop = FALSE]


    target.gene.selection_tp73ConfirmDN_context <- target.gene.selection_tp73ConfirmDN$context_data_binary
    colnames(target.gene.selection_tp73ConfirmDN_context) <- prettyIdentifierJaspar(colnames(target.gene.selection_tp73ConfirmDN_context))
    target.gene.selection_tp73ConfirmDN_context_interest <- target.gene.selection_tp73ConfirmDN_context[, cofactor.names.union, drop = FALSE]
    target.gene.selection_tp73ConfirmDN_context_interest.colsums <- colSums(target.gene.selection_tp73ConfirmDN_context_interest)
    target.gene.selection_tp73ConfirmDN_context_interest.informativeColumNames <- colnames(target.gene.selection_tp73ConfirmDN_context_interest)[target.gene.selection_tp73ConfirmDN_context_interest.colsums > 0]
    target.gene.selection_tp73ConfirmDN_context_interest.informative <- target.gene.selection_tp73ConfirmDN_context_interest[,target.gene.selection_tp73ConfirmDN_context_interest.informativeColumNames, drop = FALSE]

    # Preparation for heatmap and correlation matrix
    stopifnot(all(colnames(target.gene.selection_tp73ConfirmTA_context_interest) == colnames(target.gene.selection_tp73ConfirmDN_context_interest)))
    target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames.TAorDN <- colnames(target.gene.selection_tp73ConfirmTA_context_interest)[target.gene.selection_tp73ConfirmTA_context_interest.colsums > 0 | target.gene.selection_tp73ConfirmDN_context_interest.colsums > 0]
    target.gene.selection_tp73ConfirmDN_context_interest.informativeColumNames.TAorDN <- colnames(target.gene.selection_tp73ConfirmDN_context_interest)[target.gene.selection_tp73ConfirmTA_context_interest.colsums > 0 | target.gene.selection_tp73ConfirmDN_context_interest.colsums > 0]
    stopifnot(all(target.gene.selection_tp73ConfirmDN_context_interest.informativeColumNames.TAorDN ==target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames.TAorDN ))
    target.gene.selection_tp73ConfirmTAorDN_context_interest.informativeColumNames.TAorDN <- target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames.TAorDN

    target.gene.selection_tp73ConfirmTAorDN_context <- target.gene.selection_tp73ConfirmTAorDN$context_data_binary
    colnames(target.gene.selection_tp73ConfirmTAorDN_context) <- prettyIdentifierJaspar(colnames(target.gene.selection_tp73ConfirmTAorDN_context))
    target.gene.selection_tp73ConfirmTAorDN_context_interest <- target.gene.selection_tp73ConfirmTAorDN_context[, cofactor.names.union, drop = FALSE]
    target.gene.selection_tp73ConfirmTAorDN_context_interest.informative <- target.gene.selection_tp73ConfirmTAorDN_context_interest[,target.gene.selection_tp73ConfirmTA_context_interest.informativeColumNames.TAorDN, drop = FALSE]       

    m.cor.ta <- m.cor.dn <- NULL

    if (plot.corrplot) {
        pdf("correlation_matrix_TA_and_DN_separate.pdf", width=21, height=26)

        # Plot correlation matrix of all transcription factors (columns) against each other
        ret.corrplot.ta <- my.corrplot(m=target.gene.selection_tp73ConfirmTA_context_interest.informative,
                    relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                    title="Correlation Matrix of TA-confirmed p73 BS cofactors in genes upregulated up TA overexpression")
        m.cor.ta <- ret.corrplot.ta$corr
        ret.corrplot.dn <- my.corrplot(m=target.gene.selection_tp73ConfirmDN_context_interest.informative,
                    relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                    title="Correlation Matrix of DN-confirmed p73 BS cofactors in gene upregulated up DN overexpression")
        m.cor.dn <- ret.corrplot.dn$corr

        cofactors.union.order.ta <- rownames(ret.corrplot.ta$corr)[rownames(ret.corrplot.ta$corr) %in% cofactors.union]

        # control - should look like the first TA one
        ret.corrplot.ta.ordered.ta <- my.corrplot(m=target.gene.selection_tp73ConfirmTA_context_interest.informative[,cofactors.union.order.ta,drop=F],order="original",
                        cofactors.selection.and.order=cofactors.union.order.ta,
                        relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                        title="Control experiment - Correlation Matrix of TA -confirmed p73 BS cofactors in gene upregulated up TA overexpression - order adjusted to match TA - should look like the first TA one")
 
        # now the DN one, but ordered to match the TA one
        ret.corrplot.dn.ordered.ta <- my.corrplot(m=target.gene.selection_tp73ConfirmDN_context_interest.informative[,cofactors.union.order.ta,drop=F],order="original",
                                cofactors.selection.and.order=cofactors.union.order.ta,
                                relevant.cofactors.selection.ta=cofactor.names.TA,  relevant.cofactors.selection.dn=cofactor.names.DN, formats="none",
                                title="Correlation Matrix of DN-confirmed p73 BS cofactors in gene upregulated up DN overexpression - order adjusted to match TA")
        # Now we have two correlation matrices, one for TA and one for DN, both ordered by the same set of cofactors, in the same order - explicitly
        m.cor.ta.minus.dn <- ret.corrplot.ta.ordered.ta$corr[cofactors.union,cofactors.union]-ret.corrplot.dn.ordered.ta$corr[cofactors.union,cofactors.union]

        # Transform m.cor.ta.minus.dn into a data frame of pairs and values
        m.cor.ta.minus.dn.df <- as.data.frame(as.table(m.cor.ta.minus.dn))
        m.cor.ta.df <- as.data.frame(as.table(ret.corrplot.ta.ordered.ta$corr[cofactors.union, cofactors.union]))
        m.cor.dn.df <- as.data.frame(as.table(ret.corrplot.dn.ordered.ta$corr[cofactors.union,cofactors.union]))
        colnames(m.cor.ta.minus.dn.df) <- c("TF1", "TF2", "Difference")
        colnames(m.cor.ta.df) <- c("TF1", "TF2", "TA_Correlation")
        colnames(m.cor.dn.df) <- c("TF1", "TF2", "DN_Correlation")
        stopifnot(all(m.cor.dn.df[,1] == m.cor.ta.df[,1])) # should be TRUE
        stopifnot(all(m.cor.dn.df[,2] == m.cor.ta.df[,2])) # should be TRUE
        stopifnot(all(m.cor.ta.minus.dn.df[,1] == m.cor.ta.df[,1])) # should be TRUE
        stopifnot(all(m.cor.ta.minus.dn.df[,2] == m.cor.ta.df[,2]))= # should be TRUE
        # attach the correlation values to the difference data frame
        m.cor.ta.minus.dn.df <- cbind(m.cor.ta.minus.dn.df,correlation.TA=m.cor.ta.df[,3], correlation.DN=m.cor.dn.df[,3])
        # Optionally, filter to upper triangle (excluding diagonal) to avoid duplicates
        m.cor.ta.minus.dn.df <- m.cor.ta.minus.dn.df[as.numeric(m.cor.ta.minus.dn.df$TF1) < as.numeric(m.cor.ta.minus.dn.df$TF2), ]
        print(head(m.cor.ta.minus.dn.df))
        write.table(m.cor.ta.minus.dn.df, file="correlation_difference_TA_minus_DN.tsv", sep="\t", row.names=FALSE, quote=FALSE, dec=",",col.names=TRUE)

        #quantile(m.cor.ta.minus.dn, na.rm=TRUE)
        ret <- corrplot(corr=m.cor.ta.minus.dn/max(abs(m.cor.ta.minus.dn)), is.corr=TRUE, order="original",addCoef.col=NULL,type="lower",
                tl.cex=ifelse(nrow(m.cor.ta.minus.dn)>180,0.5,0.7), number.cex=0.35, number.digits=1,insig="n",
                title=paste0("Difference of correlation matrix of TA-confirmed p73 BS cofactors in genes upregulated up TA overexpression - DN / scaled by ",max(abs(m.cor.ta.minus.dn))))

    

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
            geom_text_repel(aes(label = label), size = 3, color = "darkred", max.overlaps = Inf) +
            labs(
            x = "Relative Abundance in TA",
            y = "Relative Abundance in DN",
            title = "TA and DN-enriched cofactors in contexts of confirmed p73 binding sites upstream of genes upregulated upon isoform overexpression"
            ) +
            theme_bw()
        print(p)

        dev.off()
        cat("I: Correlation matrix plots saved to 'correlation_matrix_TA_and_DN_separate.pdf'.\n")
    }
    #debugme <- data.frame(colnames(m),col=1 + colnames(m) %in% relevant.cofactors.selection.dn +  2*colnames(m) %in% relevant.cofactors.selection.ta +  4*colnames(m) %in% relevant.cofactors.selection.always)

    if (plot.heatmap) {

        if (FALSE) {
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input.genes <- target.gene.selection_tp73ConfirmTA$downstream_genes
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input <- target.gene.selection_tp73ConfirmTA_context_interest.informative * 1 # transform to numeric

            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input.genes <- target.gene.selection_tp73ConfirmDN$downstream_genes
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input <- target.gene.selection_tp73ConfirmDN_context_interest.informative * 1 # transform to numeric
        }

        target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input.genes <- target.gene.selection_tp73ConfirmTAorDN$downstream_genes
        target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input <- target.gene.selection_tp73ConfirmTAorDN_context_interest.informative * 1 # transform to numeric

        if (show.score) {
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input,
                Score=target.gene.selection_tp73ConfirmTA$coordinates[,c(5)]
            )
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input,
                Score=target.gene.selection_tp73ConfirmDN$coordinates[,c(5)]
            )
            target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input,
                Score=target.gene.selection_tp73ConfirmTAorDN$coordinates[,c(5)]
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
            target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input <- cbind(
                target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input,
                target.gene.selection_tp73ConfirmTAorDN$cutandrun_data[,..combinations] 
            )

        }

        # FIXME: Is this avoidable
        target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input_nonRedundant <-
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input[!duplicated(target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input.genes),]
        rownames(target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input_nonRedundant) <-
            target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input.genes[!duplicated(target.gene.selection_tp73ConfirmTA_context_interest_heatmap_input.genes)]
        target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input_nonRedundant <-
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input[!duplicated(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input.genes),]
        rownames(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input_nonRedundant) <-
            target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input.genes[!duplicated(target.gene.selection_tp73ConfirmDN_context_interest_heatmap_input.genes)]
        target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input_nonRedundant <-
            target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input[!duplicated(target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input.genes),]
        rownames(target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input_nonRedundant) <-
            target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input.genes[!duplicated(target.gene.selection_tp73ConfirmTAorDN_context_interest_heatmap_input.genes)]



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

            colors <- colorRampPalette(c("white", "gray50"))(100)

            if (nrow(n.show) < ncol(n.show)) {

                cat("D: Transpose if more rows than columns\n")
                heatmap.2(t(n.show),
                    trace="none", key=FALSE,
                    col = colors,
                    scale="none", margins=c(10,10),
                    dendrogram="none",
                    Rowv=TRUE, Colv=TRUE,
                    main="Heatmap of Gene Selection for EMT in TP73 Confirmed Contexts",
                    labRow=colnames(n.show),
                    colRow=1 + colnames(n.show) %in% relevant.cofactors.selection.ta + 
                        2*colnames(n.show) %in% relevant.cofactors.selection.dn + 
                        4*colnames(n.show) %in% relevant.cofactors.selection.always,
                    labCol=rownames(n.show),
                    colCol=1 + rownames(n.show) %in% target.gene.selection.TA + 2 * rownames(n.show) %in% target.gene.selection.DN + 4 * rownames(n.show) %in% target.gene.selection.always,
                    cexRow=ifelse(ncol(n.show)>50,0.6,0.75),
                    cexCol=ifelse(nrow(n.show)>180,ifelse(nrow(n.show)>250,0.35,0.5),0.75)
                )
             
            } else {

                heatmap.2(n.show,
                    trace="none", key=FALSE,
                    col = colors,
                    scale="none", margins=c(10,10),
                    dendrogram="none",
                    Rowv=TRUE, Colv=TRUE,
                    main="Heatmap of Gene Selection for EMT in TP73 Confirmed Contexts",
                    labRow=n,
                    colRow=1 + rownames(n.show) %in% target.gene.selection.TA + 2 * rownames(n.show) %in% target.gene.selection.DN + 4 * rownames(n.show) %in% target.gene.selection.always,
                    labCol=colnames(n.show),
                    colCol=1 + colnames(n.show) %in% relevant.cofactors.selection.ta + 
                        2*colnames(n.show) %in% relevant.cofactors.selection.dn + 
                        4*colnames(n.show) %in% relevant.cofactors.selection.always,
                    cexRow=ifelse(nrow(n.show)>50,0.6,0.75),
                    cexCol=ifelse(ncol(n.show)>180,ifelse(ncol(n.show)>250,0.35,0.5),0.75)
                )
            }

            # Add legend for red label color
            legend("topright", legend=c("DN","TA","DN+TA"), text.col=c(1+1,1+2,1+1+2), bty="n", cex=0.9)
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
            m.cor.dn=m.cor.dn,
            cofactors.TA=cofactor.names.TA,
            cofactors.DN=cofactor.names.DN,
            cofactors.intersection=intersect(cofactor.names.TA,cofactor.names.DN),
            cofactors.union=cofactors.union
        )
    )
}


pdf("correlation_upregulated_TA_DN.pdf", width=12, height=10)
r <- my.heatmap.and.corrplot(target.gene.selection.TA=e.ta.up, target.gene.selection.DN=e.dn.up,formats="none",combined.expression.data = combined.expression.data)
dev.off()
