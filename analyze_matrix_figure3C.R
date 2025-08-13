# Figure 3C - preparation

# TA-/DN-specific gene expresssion change from Venn diagrams

require(xlsx)
#e <- read.xlsx("GeneLists/Venn Diagram (Gene mit hoher p73 Bindung) für Korrelationsgrafik.xlsx",header=F,sheetIndex=2,endRow=8)
e <- read.xlsx("GeneLists/Genlisten_fuer_p73_Bindung_EMT_Gene_GSEA_20250721.xlsx",header=F,sheetIndex=2,startRow=2,endRow=4)
e.ta.up<-unlist(e[1,!is.na(e[1,]),drop=T])[-1]
e.dn.up<-unlist(e[2,!is.na(e[2,]),drop=T])[-1]
e.both.up<-unlist(e[3,!is.na(e[3,]),drop=T])[-1]

target.gene.selection.TA <- e.ta.up
target.gene.selection.DN <- e.dn.up

if (0 == length(ls(pattern="all_inPromoter_tp73ConfirmAny"))) {
    cat("I: No data found for all_inPromoter_tp73ConfirmAny, retrieving context data...\n")
    source("analyze_matrix_function_retrieve_context.R")
}

pdf("correlation_upregulated_TA_DN.pdf", width=12, height=10)
r <- my.heatmap.and.corrplot(target.gene.selection.TA=e.ta.up,
                             target.gene.selection.DN=e.dn.up,formats="none",
                             combined.expression.data = combined.expression.data,
                             plot.corrplot=TRUE,
                             plot.heatmap=FALSE,
                             title.heatmap="Cofactors of genes upregulated in TP73 Confirmed Contexts")
dev.off()
