
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

#
#  Figure 4 - end
#