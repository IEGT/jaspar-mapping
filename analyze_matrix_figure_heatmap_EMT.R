# Manually copied from Nico's Excel table

# EMT-enriched genes, upregulated

TA.only <- c("CRLF1","SERPINE1","MATN2","GPC1","EDIL3","VCAN","MYLK","IGFBP3","BMP1",
             "SDC4","SLIT2","TGFB1","LAMA3","PTHLH","TAGLN","FBLN2","RHOB","ACTA2","TPM1","EFEMP2",
             "HTRA1","GEM","QSOX1","COMP","CALD1","TNFRSF12A","CD44","COPA")
DN.only <- c("SFRP1","DKK1","COL4A1","PMP22","FN1","EMP3","COLGALT1","DAB2")
TA.and.DN <- c("LGALS1","PMEPA1","SDC1","SERPINH1","CADM1","ID2","THBS1","GADD45B",
             "ITGAV","FLNA","PPIB","MCM7","SPARC","SLC6A8","PLOD1")


pdf("heatmap_emt.pdf")


my.heatmap.and.corrplot(
    base.filename="heatmap_emt",
    title.heatmap="Cofactors enriched for EMT-associated genes",
    target.gene.selection.TA=TA.only,
    target.gene.selection.DN=DN.only,
    target.gene.selection.TAandDN=TA.and.DN,
    formats="pdf",
    plot.corrplot=FALSE,
    plot.heatmap=TRUE,
    show.cutandrun=FALSE,
    show.score=FALSE,
    filter.for.strong.enrichment.num=25
)

dev.off()