
#
# Heatmap for gene selection
#

# Retrieves all contexts for binding sites in promoter regions, which is the background for the analysis.
all_inPromoter_tp73ConfirmAny <- retrieve_context_data_by_chromosome(NULL, confirmation=c("tp73","promoter"),TA.or.DN="any")
all_inPromoter_tp73ConfirmAny.colSums <- colSums(all_inPromoter_tp73ConfirmAny$context_data,na.rm=T)

require(ggplot2)

#gene.selection <- promoterBedTables$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.promoter.bed$Gene
#gene.selection <- gene.selection.EMT.TAorDN

#gene.selection <- e.ta.up
gene.selection <- e.dn.up



combined.expression.data.selected <- combined.expression.data[,"Gene Symbol"] %in% gene.selection
gene.selection_tp73ConfirmAny <- retrieve_context_data_by_chromosome(combined.expression.data.selected, confirmation=c("tp73"),TA.or.DN="any")


cat("I: Distribution of number of matches by chromosome:\n")
print(gene.selection_tp73ConfirmAny$num_matches_by_chromosome)
cat("I: Total number of matches across all chromosomes: ",sum(gene.selection_tp73ConfirmAny$num_matches_by_chromosome),"\n",sep="")

gene.selection_tp73ConfirmAny.colSums <- colSums(gene.selection_tp73ConfirmAny$context_data,na.rm=T)

print(gene.selection_tp73ConfirmAny.colSums)

gene.vs.all_inPromoter_tp73ConfirmAny <- cbind(
    mean.EMT=gene.selection_tp73ConfirmAny$mean_total,
    mean.background=all_inPromoter_tp73ConfirmAny$mean_total,
    ratio=gene.selection_tp73ConfirmAny$mean_total / all_inPromoter_tp73ConfirmAny$mean_total
)

tf.patterns.of.interest <- paste0(
            #"^(jun|sp1|rest|yap1|yy1|e2f1|nfkb1|tp53|tp63|tp73|",
            #"BZIP43|IRF2|RXRA--VDR|pan|HAP1|Nr2F6|RORA|E2F3|E2F7|THI2|SPL15|TGIF1|MYB52|ZBED1|NAC083", #  paste(sapply(strsplit(rownames(head(ta_vs_effects_tp73Confirm_posConfirm.sorted,15)),split=" "),function(X) X[1]),collapse="|")
            #"NR3C1|TP63|RREB1|ATHB-9|GLIS1|PLAG1|E2FC|MYB15|DYT1|bZIP911|odd|RARA|BZIP42|BPC5|TCP14", # paste(sapply(strsplit(rownames(tail(ta_vs_effects_tp73Confirm_posConfirm.sorted[!is.na(ta_vs_effects_tp73Confirm_posConfirm.sorted[,"ratio"]),],15)),split=" "),function(X) X[1]),collapse="|")
            "RARA|GLIS1|PLAG1|RREB1|TP63|NR3C1|SRF|NR3C2|EWSR1-FLI1|THRB|TFAP2E|GLIS3|TFAP2C|RXRB|KLF14|",
            "IRF2|RXRA|RORA|E2F7|E2F3|TFIF1|ZBED1|IRF4|RFX3|MEF2D|STAT1|MEF2B|NFIC|IRF8|CREB3L1",
        ") ")

relevant.cofactors.selection.dn

gene.vs.all_inPromoter_tp73ConfirmAny.sorted <- gene.vs.all_inPromoter_tp73ConfirmAny[order(gene.vs.all_inPromoter_tp73ConfirmAny[, "ratio"],decreasing=T), ]
rownames(gene.vs.all_inPromoter_tp73ConfirmAny.sorted) <- prettyIdentifierJaspar(rownames(gene.vs.all_inPromoter_tp73ConfirmAny.sorted))
if (FALSE) {
    gene.vs.all_inPromoter_tp73ConfirmAny.rownames.interest <- c(rownames(gene.vs.all_inPromoter_tp73ConfirmAny.sorted)[is.human.jaspar.id(rownames(gene.vs.all_inPromoter_tp73ConfirmAny.sorted))][1:50],
        grep(x=prettyIdentifierJaspar(rownames(gene.vs.all_inPromoter_tp73ConfirmAny.sorted)),pattern=tf.patterns.of.interest,ignore.case=T,value=T)
        )
    gene.vs.all_inPromoter_tp73ConfirmAny.rownames.interest <- unique(gene.vs.all_inPromoter_tp73ConfirmAny.rownames.interest)
}
gene.vs.all_inPromoter_tp73ConfirmAny.rownames.interest <- c(relevant.cofactors.selection.dn,relevant.cofactors.selection.ta)

combined.expression.data.reduced <- combined.expression.data[combined.expression.data.selected,]

nrow(gene.selection_tp73ConfirmAny$context_data)

# Print all possible six combinations of "pos"/"tp73" with "skmel29_2" and "TA"/"DN"/"GFP"
combinations <- sort(as.vector(outer(c("pos", "tp73"), c("TA", "DN", "GFP"), function(x, y) paste(x, "skmel29_2", y, sep = "_"))))
print(combinations)

gene.selection_tp73ConfirmAny_context <- gene.selection_tp73ConfirmAny$context_data
colnames(gene.selection_tp73ConfirmAny_context) <- prettyIdentifierJaspar(colnames(gene.selection_tp73ConfirmAny_context))
gene.selection_tp73ConfirmAny_context_interest <- gene.selection_tp73ConfirmAny_context[, ..gene.vs.all_inPromoter_tp73ConfirmAny.rownames.interest, drop = FALSE]
gene.selection_tp73ConfirmAny_context_interest.colsums <- colSums(gene.selection_tp73ConfirmAny_context_interest)
gene.selection_tp73ConfirmAny_context_interest.informativeColumNames <- colnames(gene.selection_tp73ConfirmAny_context_interest)[gene.selection_tp73ConfirmAny_context_interest.colsums > 0]
gene.selection_tp73ConfirmAny_context_interest.informative <- gene.selection_tp73ConfirmAny_context_interest[,..gene.selection_tp73ConfirmAny_context_interest.informativeColumNames, drop = FALSE]

gene.selection_tp73ConfirmAny_context_interest_heatmap_input <- cbind(Gene=gene.selection_tp73ConfirmAny$downstream_genes,
    #ene.selection_tp73ConfirmAny$coordinates[,c(5)],
      #gene.selection_tp73ConfirmAny$cutandrun_data[,..combinations],
      gene.selection_tp73ConfirmAny_context_interest.informative
)
gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant <- gene.selection_tp73ConfirmAny_context_interest_heatmap_input[!duplicated(gene.selection_tp73ConfirmAny_context_interest_heatmap_input[,"Gene"]),]

rownames(gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant) <- gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant$Gene

require(gplots)

pdf("heatmap_gene_selection_tp73ConfirmAny_testing.pdf", width=18, height=8)
n <- gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant$Gene
n.show <- as.matrix(gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant[,-c(1)])
n.show <- apply(n.show, 2, function(x) {
    x[is.na(x)] <- 0
    x = x/max(x, na.rm=TRUE)  # Scale to [0, 1]
    return(x)
})

if (FALSE) {
    n.show <- apply(n.show, c(1,2), function(x) {
        if (x>0.33) x=1
        else x=0
        x
    })
}
heatmap.2(n.show,
          trace="none", key=FALSE,
          col = colorRampPalette(c("white", "gray50", "black"))(100),
          scale="none", margins=c(10,10),
          dendrogram="none",
          Rowv=TRUE, Colv=TRUE,
          main="Heatmap of Gene Selection for EMT in TP73 Confirmed Contexts",
          labRow=n,
          colRow=1 + n %in% gene.selection.EMT.DNonly,
          labCol=colnames(n.show),
          colCol=1 + colnames(n.show) %in% relevant.cofactors.selection.dn,
          cexRow=0.75, cexCol=0.75
)
# Add legend for red label color
legend("topright", legend=c("Red label: selection for DN","Black: selection for TA"), text.col=c("red","black"), bty="n", cex=0.9)
dev.off()
