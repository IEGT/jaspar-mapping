
# 
# Preparation of Vulcano Plot
#

# Also defines
#  relevant.cofactors.selection.ta
#  relevant.cofactors.selection.dn
#  relevant.cofactors.selection

require(ggplot2)

threshold.num.reads <- 3

#rm(m); m <- m.contexts[[22]]
# Derive the number of TFBS per cofactor in each context, not counting multiple occurrences in same TFBS context
m.contexts.num.binary.sum <- sapply(m.contexts, function(m) {
    if (is.null(m)) {
        return(NA)
    }
    # Sum the number of non-NA entries in columns with "_Score" in their name
    #cols.Score.name <- grep("_Score", colnames(m), value = TRUE)
    cols.Score.index <- grep("_NumInWindow$", colnames(m), value = FALSE)
    colSums(m[,..cols.Score.index,drop = FALSE] > 0, na.rm = TRUE) # number of TFBS featuring that cofactor
})
print(dim(m.contexts.num.binary.sum))

m.contexts.num.binary.sum.total <- rowSums(m.contexts.num.binary.sum, na.rm = TRUE)
m.contexts.num.tfbs <- sapply(m.contexts, nrow)
m.contexts.num.tfbs.total <- sum(m.contexts.num.tfbs, na.rm = TRUE)

#rm(m); m <- m.contexts[[22]]

count.cofactors.per.tfbs <- function(confirmation=NULL, TA.or.DN="any", threshold.num.reads=1 ) {
    sapply(m.contexts, function(m) {

        if (is.null(m)) {
            return(NA)
        }
        m.context.extra.check <- rep(TRUE, nrow(m))
        cols.Score.index <- grep("_NumInWindow$", colnames(m), value = FALSE)

        if ("promoter" %in% confirmation) {
            m.context.extra.check <- m.context.extra.check & m[, "InPromoter"]
        }

        if ("tp73" %in% confirmation) {

            if (TA.or.DN == "TA") {
                m.context.extra.check <- m.context.extra.check & m$"tp73_skmel29_2_TA" >= threshold.num.reads
            } else if (TA.or.DN == "DN") {
                m.context.extra.check <- m.context.extra.check & m$"tp73_skmel29_2_DN" >= threshold.num.reads
            } else if (TA.or.DN == "GFP") {
                m.context.extra.check <- m.context.extra.check & m$"tp73_skmel29_2_GFP" >= threshold.num.reads
            } else if (TA.or.DN == "any") {
                m.context.extra.check <- m.context.extra.check & (
                    m$"tp73_skmel29_2_TA" >= threshold.num.reads |
                    m$"tp73_skmel29_2_DN" >= threshold.num.reads |
                    m$"tp73_skmel29_2_GFP" >= threshold.num.reads
                )
            } else if (TA.or.DN == "none") {
                # doing nothing, as we want to check all rows
            } else {
                stop("E: Invalid TA.or.DN method specified.\n")
            }
        }

        if ("pos" %in% confirmation) {
            if (TA.or.DN == "TA") {
                m.context.extra.check <- m.context.extra.check & ( m$"pos_skmel29_2_TA" >= threshold.num.reads )
            } else if (TA.or.DN == "DN") {
                m.context.extra.check <- m.context.extra.check & ( m$"pos_skmel29_2_DN" >= threshold.num.reads )
            } else if (TA.or.DN == "GFP") {
                m.context.extra.check <- m.context.extra.check & ( m$"pos_skmel29_2_GFP" >= threshold.num.reads )
            } else if (TA.or.DN == "any") {
                m.context.extra.check <- m.context.extra.check & (
                    m$"pos_skmel29_2_DN" >= threshold.num.reads |
                    m$"pos_skmel29_2_TA" >= threshold.num.reads |
                    m$"pos_skmel29_2_GFP" >= threshold.num.reads
                    )
            } else if (TA.or.DN == "none") {
                # doing nothing, as we want to check all rows
            } else {
                stop("E: Invalid TA.or.DN method specified.\n")
            }
        }

        if ("pos-up" %in% confirmation) {
            if (TA.or.DN == "TA") {
                m.context.extra.check <- m.context.extra.check & ( m$"pos_skmel29_2_TA" >= m$"pos_skmel29_2_GFP" )
            } else if (TA.or.DN == "DN") {
                m.context.extra.check <- m.context.extra.check & ( m$"pos_skmel29_2_DN" >= m$"pos_skmel29_2_GFP" )
            } else if (TA.or.DN == "any") {
                m.context.extra.check <- m.context.extra.check & (
                    m$"pos_skmel29_2_DN" >= m$"pos_skmel29_2_GFP" |
                    m$"pos_skmel29_2_TA" >= m$"pos_skmel29_2_GFP"
                    )
            } else if (TA.or.DN == "none") {
                # doing nothing, as we want to check all rows
            } else {
                stop("E: Invalid TA.or.DN method specified.\n")
            }
        }

        m <- m[m.context.extra.check, , drop = FALSE]
        m <- m[, ..cols.Score.index, drop = FALSE]
        colSums(m > 0, na.rm = TRUE) # count the number of TFBS that are confirmed by TA, GFP or DN
    })
}

#m.contexts.num.binary.confirmed.taOrDn.sum <- count.cofactors.per.tfbs(confirmation=c("tp73"), TA.or.DN="any", threshold.num.reads=threshold.num.reads)
#m.contexts.num.binary.confirmed.taOrDn.sum.total <- rowSums(m.contexts.num.binary.confirmed.taOrDn.sum, na.rm = TRUE)

count.tfbs <- function(confirmation=NULL, TA.or.DN="any", threshold.num.reads=1 ) {
    sapply(m.contexts, function(m) {

        if (is.null(m)) {
            return(NA)
        }
        m.context.extra.check <- rep(TRUE, nrow(m))
        cols.Score.index <- grep("_NumInWindow$", colnames(m), value = FALSE)

        if ("promoter" %in% confirmation) {
            m.context.extra.check <- m.context.extra.check & m[, "InPromoter"]
        }

        if ("tp73" %in% confirmation) {

            if (TA.or.DN == "TA") {
                m.context.extra.check <- m.context.extra.check & ( m$"tp73_skmel29_2_TA" >= threshold.num.reads )
            } else if (TA.or.DN == "DN") {
                m.context.extra.check <- m.context.extra.check & ( m$"tp73_skmel29_2_DN" >= threshold.num.reads )
            } else if (TA.or.DN == "GFP") {
                m.context.extra.check <- m.context.extra.check & ( m$"tp73_skmel29_2_GFP" >= threshold.num.reads )
            } else if (TA.or.DN == "any") {
                m.context.extra.check <- m.context.extra.check & (
                    m$"tp73_skmel29_2_TA" >= threshold.num.reads |
                    m$"tp73_skmel29_2_DN" >= threshold.num.reads |
                    m$"tp73_skmel29_2_GFP" >= threshold.num.reads
                )
            } else if (TA.or.DN == "none") {
                # doing nothing, as we want to check all rows
            } else {
                stop("E: Invalid TA.or.DN method specified.\n")
            }
        }

        if ("pos" %in% confirmation) {
            if (TA.or.DN == "TA") {
                m.context.extra.check <- m.context.extra.check & ( m$"pos_skmel29_2_TA" >= threshold.num.reads )
            } else if (TA.or.DN == "DN") {
                m.context.extra.check <- m.context.extra.check & ( m$"pos_skmel29_2_DN" >= threshold.num.reads )
            } else if (TA.or.DN == "GFP") {
                m.context.extra.check <- m.context.extra.check & ( m$"pos_skmel29_2_GFP" >= threshold.num.reads )
            } else if (TA.or.DN == "any") {
                m.context.extra.check <- m.context.extra.check & (
                    m$"pos_skmel29_2_DN" >= threshold.num.reads |
                    m$"pos_skmel29_2_TA" >= threshold.num.reads |
                    m$"pos_skmel29_2_GFP" >= threshold.num.reads
                    )
            } else if (TA.or.DN == "none") {
                # doing nothing, as we want to check all rows
            } else {
                stop("E: Invalid TA.or.DN method specified.\n")
            }
        }

        if ("pos-up" %in% confirmation) {
            if (TA.or.DN == "TA") {
                m.context.extra.check <- m.context.extra.check & ( m$"pos_skmel29_2_TA" >= m$"pos_skmel29_2_GFP" )
            } else if (TA.or.DN == "DN") {
                m.context.extra.check <- m.context.extra.check & ( m$"pos_skmel29_2_DN" >= m$"pos_skmel29_2_GFP" )
            } else if (TA.or.DN == "any") {
                m.context.extra.check <- m.context.extra.check & (
                    m$"pos_skmel29_2_DN" >= m$"pos_skmel29_2_GFP" |
                    m$"pos_skmel29_2_TA" >= m$"pos_skmel29_2_GFP"
                    )
            } else if (TA.or.DN == "none") {
                # doing nothing, as we want to check all rows
            } else {
                stop("E: Invalid TA.or.DN method specified.\n")
            }
        }

        sum(m.context.extra.check, na.rm = TRUE) # count the number of TFBS that are confirmed by TA, GFP or DN - not the number of cofactors
    })
}

#m.contexts.num.tfbs.confirmed.taOrDn2 <- count.tfbs(confirmation=c("tp73"), TA.or.DN="any", threshold.num.reads=threshold.num.reads)
#m.contexts.num.tfbs.confirmed.taOrDn.total <- sum(m.contexts.num.tfbs.confirmed.taOrDn, na.rm = TRUE)
#frequency.cofactors.taOrDn <- m.contexts.num.binary.confirmed.taOrDn.sum.total/m.contexts.num.tfbs.confirmed.taOrDn.total
#names(frequency.cofactors.taOrDn) <- prettyIdentifierJaspar(names(frequency.cofactors.taOrDn))

m.contexts.num.binary.confirmed.ta.sum.total <- rowSums(count.cofactors.per.tfbs(confirmation=c("tp73"), TA.or.DN="TA", threshold.num.reads=threshold.num.reads), na.rm = TRUE)
m.contexts.num.tfbs.confirmed.ta.total <- sum(count.tfbs(confirmation=c("tp73"), TA.or.DN="TA", threshold.num.reads=threshold.num.reads))
frequency.cofactors.ta <- m.contexts.num.binary.confirmed.ta.sum.total/m.contexts.num.tfbs.confirmed.ta.total
names(frequency.cofactors.ta) <- prettyIdentifierJaspar(names(frequency.cofactors.ta))


m.contexts.num.binary.confirmed.dn.sum.total <- rowSums(count.cofactors.per.tfbs(confirmation=c("tp73"), TA.or.DN="DN", threshold.num.reads=threshold.num.reads), na.rm = TRUE)
m.contexts.num.tfbs.confirmed.dn.total <- sum(count.tfbs(confirmation=c("tp73"), TA.or.DN="DN", threshold.num.reads=threshold.num.reads))
frequency.cofactors.dn <- m.contexts.num.binary.confirmed.dn.sum.total/m.contexts.num.tfbs.confirmed.dn.total
names(frequency.cofactors.dn) <- prettyIdentifierJaspar(names(frequency.cofactors.dn))

log.ratio.confirmed.ta.vs.dn <- log2(frequency.cofactors.ta / frequency.cofactors.dn)
names(log.ratio.confirmed.ta.vs.dn) <- prettyIdentifierJaspar(names(log.ratio.confirmed.ta.vs.dn))
log.ratio.confirmed.ta.vs.dn.human <- log.ratio.confirmed.ta.vs.dn[is.human.jaspar.id(names(log.ratio.confirmed.ta.vs.dn))]

frequency.cofactors.initial <- m.contexts.num.binary.sum.total/m.contexts.num.tfbs.total

#frequency.cofactors.taOrDn.human <- frequency.cofactors.taOrDn[is.human.jaspar.id(names(frequency.cofactors.taOrDn))]
frequency.cofactors.ta.human <- frequency.cofactors.ta[is.human.jaspar.id(names(frequency.cofactors.ta))]
frequency.cofactors.dn.human <- frequency.cofactors.dn[is.human.jaspar.id(names(frequency.cofactors.dn))]

#enrichment.taOrDn <- log2(frequency.cofactors.taOrDn / frequency.cofactors.initial)
## enrichment.taOrDn <- log2((m.contexts.num.binary.confirmed.taOrDn.sum.total/m.contexts.num.tfbs.confirmed.taOrDn.total) / (m.contexts.num.binary.sum.total/m.contexts.num.tfbs.total))
#names(enrichment.taOrDn) <- prettyIdentifierJaspar(names(enrichment.taOrDn))
#enrichment.taOrDn.human <- enrichment.taOrDn[is.human.jaspar.id(names(enrichment.taOrDn))]

enrichment.ta <- log2(frequency.cofactors.ta / frequency.cofactors.initial)
names(enrichment.ta) <- prettyIdentifierJaspar(names(enrichment.ta))
enrichment.ta.human <- enrichment.ta[is.human.jaspar.id(names(enrichment.ta))]

enrichment.dn <- log2(frequency.cofactors.dn / frequency.cofactors.initial)
names(enrichment.dn) <- prettyIdentifierJaspar(names(enrichment.dn))
enrichment.dn.human <- enrichment.dn[is.human.jaspar.id(names(enrichment.dn))]

# Enrichments for methylation for CUT&RUN-confirmed p73 binding sites
#m.contexts.num.binary.confirmed.taOrDn.methylation.sum.total <- rowSums(count.cofactors.per.tfbs(confirmation=c("tp73","pos"), TA.or.DN="any", threshold.num.reads=threshold.num.reads), na.rm = TRUE))
m.contexts.num.binary.confirmed.ta.methylation.sum.total <- rowSums(count.cofactors.per.tfbs(confirmation=c("tp73","pos"), TA.or.DN="TA", threshold.num.reads=threshold.num.reads), na.rm = TRUE)
m.contexts.num.binary.confirmed.dn.methylation.sum.total <- rowSums(count.cofactors.per.tfbs(confirmation=c("tp73","pos"), TA.or.DN="DN", threshold.num.reads=threshold.num.reads), na.rm = TRUE)

#m.contexts.num.tfbs.confirmed.taOrDn.methylation <- sum(count.tfbs(confirmation=c("tp73","pos"), TA.or.DN="any", threshold.num.reads=threshold.num.reads), na.rm = TRUE)
m.contexts.num.tfbs.confirmed.ta.methylation.total <- sum(count.tfbs(confirmation=c("tp73","pos"), TA.or.DN="TA", threshold.num.reads=threshold.num.reads), na.rm = TRUE)
m.contexts.num.tfbs.confirmed.dn.methylation.total <- sum(count.tfbs(confirmation=c("tp73","pos"), TA.or.DN="DN", threshold.num.reads=threshold.num.reads), na.rm = TRUE)

#frequency.cofactors.taOrDn.methylation <- m.contexts.num.binary.confirmed.taOrDn.methylation.sum.total / m.contexts.num.tfbs.confirmed.taOrDn.methylation.total
#names(frequency.cofactors.taOrDn.methylation) <- prettyIdentifierJaspar(names(frequency.cofactors.taOrDn.methylation))
frequency.cofactors.ta.methylation <- m.contexts.num.binary.confirmed.ta.methylation.sum.total / m.contexts.num.tfbs.confirmed.ta.methylation.total
names(frequency.cofactors.ta.methylation) <- prettyIdentifierJaspar(names(frequency.cofactors.ta.methylation))
frequency.cofactors.dn.methylation <- m.contexts.num.binary.confirmed.dn.methylation.sum.total / m.contexts.num.tfbs.confirmed.dn.methylation.total
names(frequency.cofactors.dn.methylation) <- prettyIdentifierJaspar(names(frequency.cofactors.dn.methylation))

#enrichment.taOrDn.methylation.vs.taOrDn.tp73 <- log2(frequency.cofactors.taOrDn.methylation / frequency.cofactors.taOrDn)
enrichment.ta.methylation.vs.ta.tp73 <- log2(frequency.cofactors.ta.methylation / frequency.cofactors.ta)
enrichment.dn.methylation.vs.dn.tp73 <- log2(frequency.cofactors.dn.methylation / frequency.cofactors.dn)
#enrichment.taOrDn.methylation.vs.taOrDn.tp73.human <- enrichment.taOrDn.methylation.vs.taOrDn.tp73[is.human.jaspar.id(names(enrichment.taOrDn.methylation.vs.taOrDn.tp73))]
enrichment.ta.methylation.vs.ta.tp73.human <- enrichment.ta.methylation.vs.ta.tp73[is.human.jaspar.id(names(enrichment.ta.methylation.vs.ta.tp73))]
enrichment.dn.methylation.vs.dn.tp73.human <- enrichment.dn.methylation.vs.dn.tp73[is.human.jaspar.id(names(enrichment.dn.methylation.vs.dn.tp73))]

log.ratio.confirmed.ta.vs.dn.methylation <- log2(enrichment.ta.methylation.vs.ta.tp73 / enrichment.dn.methylation.vs.dn.tp73)
names(log.ratio.confirmed.ta.vs.dn.methylation) <- prettyIdentifierJaspar(names(log.ratio.confirmed.ta.vs.dn.methylation))
log.ratio.confirmed.ta.vs.dn.methylation.human <- log.ratio.confirmed.ta.vs.dn.methylation[is.human.jaspar.id(names(log.ratio.confirmed.ta.vs.dn.methylation))]

# Again for methylation for CUT&RUN-confirmed p73 binding sites, but only accepting TFBS with up-regulated methylation

#m.contexts.num.binary.confirmed.taOrDn.methylation.is.up.sum.total <- rowSums(count.cofactors.per.tfbs(confirmation=c("tp73","pos"), TA.or.DN="any", threshold.num.reads=threshold.num.reads), na.rm = TRUE)
m.contexts.num.binary.confirmed.ta.methylation.is.up.sum.total <- rowSums(count.cofactors.per.tfbs(confirmation=c("tp73","pos"), TA.or.DN="TA", threshold.num.reads=threshold.num.reads), na.rm = TRUE)
m.contexts.num.binary.confirmed.dn.methylation.is.up.sum.total <- rowSums(count.cofactors.per.tfbs(confirmation=c("tp73","pos"), TA.or.DN="DN", threshold.num.reads=threshold.num.reads), na.rm = TRUE)
#m.contexts.num.tfbs.confirmed.taOrDn.methylation.is.up.total <- sum(count.tfbs(confirmation=c("tp73","pos"), TA.or.DN="any", threshold.num.reads=threshold.num.reads), na.rm=T)
m.contexts.num.tfbs.confirmed.ta.methylation.is.up.total <- sum(count.tfbs(confirmation=c("tp73","pos"), TA.or.DN="TA", threshold.num.reads=threshold.num.reads), na.rm=T)
m.contexts.num.tfbs.confirmed.dn.methylation.is.up.total <- sum(count.tfbs(confirmation=c("tp73","pos"), TA.or.DN="DN", threshold.num.reads=threshold.num.reads), na.rm=T)

#frequency.cofactors.taOrDn.methylation.is.up <- m.contexts.num.binary.confirmed.taOrDn.methylation.is.up.sum.total / m.contexts.num.tfbs.confirmed.taOrDn.methylation.is.up.total
#names(frequency.cofactors.taOrDn.methylation.is.up) <- prettyIdentifierJaspar(names(frequency.cofactors.taOrDn.methylation.is.up))
frequency.cofactors.ta.methylation.is.up <- m.contexts.num.binary.confirmed.ta.methylation.is.up.sum.total / m.contexts.num.tfbs.confirmed.ta.methylation.is.up.total
names(frequency.cofactors.ta.methylation.is.up) <- prettyIdentifierJaspar(names(frequency.cofactors.ta.methylation.is.up))
frequency.cofactors.dn.methylation.is.up <- m.contexts.num.binary.confirmed.dn.methylation.is.up.sum.total / m.contexts.num.tfbs.confirmed.dn.methylation.is.up.total
names(frequency.cofactors.dn.methylation.is.up) <- prettyIdentifierJaspar(names(frequency.cofactors.dn.methylation.is.up))
#frequency.cofactors.taOrDn.methylation.is.up.human <- frequency.cofactors.taOrDn.methylation.is.up[is.human.jaspar.id(names(frequency.cofactors.taOrDn.methylation.is.up))]
frequency.cofactors.ta.methylation.is.up.human <- frequency.cofactors.ta.methylation.is.up[is.human.jaspar.id(names(frequency.cofactors.ta.methylation.is.up))]
frequency.cofactors.dn.methylation.is.up.human <- frequency.cofactors.dn.methylation.is.up[is.human.jaspar.id(names(frequency.cofactors.dn.methylation.is.up))]


#enrichment.taOrDn.methylation.is.up.vs.taOrDn.tp73 <- log2(frequency.cofactors.taOrDn.methylation.is.up / frequency.cofactors.taOrDn)
enrichment.ta.methylation.is.up.vs.ta.tp73 <- log2(frequency.cofactors.ta.methylation.is.up / frequency.cofactors.ta)
enrichment.dn.methylation.is.up.vs.dn.tp73 <- log2(frequency.cofactors.dn.methylation.is.up / frequency.cofactors.dn)
#enrichment.taOrDn.methylation.is.up.vs.taOrDn.tp73.human <- enrichment.taOrDn.methylation.is.up.vs.taOrDn.tp73[is.human.jaspar.id(names(enrichment.taOrDn.methylation.is.up.vs.taOrDn.tp73))]
enrichment.ta.methylation.is.up.vs.ta.tp73.human <- enrichment.ta.methylation.is.up.vs.ta.tp73[is.human.jaspar.id(names(enrichment.ta.methylation.is.up.vs.ta.tp73))]
enrichment.dn.methylation.is.up.vs.dn.tp73.human <- enrichment.dn.methylation.is.up.vs.dn.tp73[is.human.jaspar.id(names(enrichment.dn.methylation.is.up.vs.dn.tp73))]

log.ratio.confirmed.ta.vs.dn.methylation.is.up <- enrichment.ta.methylation.is.up.vs.ta.tp73 - enrichment.dn.methylation.is.up.vs.dn.tp73
log.ratio.confirmed.ta.vs.dn.methylation.is.up.human <- log.ratio.confirmed.ta.vs.dn.methylation.is.up[is.human.jaspar.id(names(log.ratio.confirmed.ta.vs.dn.methylation.is.up))]


# Volcano-like scatter plot: enrichment.taOrDn (X) vs frequency.cofactors.taOrDn (Y)
require(ggplot2)
require(ggrepel)
#log.ratio.confirmed.ta.vs.dn.human.ranked <- rank(log.ratio.confirmed.ta.vs.dn.human, ties.method = "first")
#log.ratio.confirmed.ta.vs.dn.human.ranked <- log.ratio.confirmed.ta.vs.dn.human.ranked - mean(log.ratio.confirmed.ta.vs.dn.human.ranked)

my.volcano.plot <- function(data,title="Volcano Plot: Enrichment vs Frequency",
                            xlab="log2 Enrichment (confirmed / initial)",
                            ylab="Frequency (confirmed / initial)",n.tfbs.prior=NA,n.tfbs.post=NA) {

    ggplot(data, aes(x=Enrichment, y=Frequency, label=TF)) +
        geom_point(#aes(fill=TAvsDN),
                color="black", size=1.5, shape=21, alpha=0.75) +
  #      scale_fill_gradient2(
  #          low = "blue", mid = "white", high = "red",
  #          midpoint = 0, name = "log2(TA/DN)"
  #      ) +
        geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "red") +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
        ggrepel::geom_text_repel(
            data = subset(data, Enrichment < -0.1 | Enrichment > 0.1 | TF=="TP73"),
            size = 2.5, max.overlaps = Inf, min.segment.length = 0,
            box.padding = 0.3, point.padding = 0.2, segment.color = "grey90"
        ) +
        labs(
            title=title,
            x=xlab,
            y=ylab
        ) +
        theme_minimal()
    
}

#enrichment.ta.human>0 & enrichment.dn.human<0

pdf.volcano.filename <- paste("volcano_plot_enrichment_vs_frequency_treshold_",threshold.num.reads,".pdf",sep="")
pdf(pdf.volcano.filename, width=11, height=8)

# p73 binding enrichment plots

if (FALSE) {
    volcano_data.taOrDn <- data.frame(
        TF = sapply(strsplit(x=names(enrichment.taOrDn.human),split="[_ ]",fixed=FALSE), function(x) x[1]),
        Enrichment = enrichment.taOrDn.human,
        Frequency = frequency.cofactors.taOrDn.human,
        TAvsDN = log.ratio.confirmed.ta.vs.dn.human
    )
    my.volcano.plot(volcano_data.taOrDn,title=paste0("Enrichment for p73 binding (TA, DN or GFP confirmed by ",threshold.num.reads,"+ reads)"))
}

volcano_data.ta <- data.frame(
    TF = sapply(strsplit(x=names(enrichment.ta.human),split="[_ ]",fixed=FALSE), function(x) x[1]),
    Enrichment = enrichment.ta.human,
    Frequency = frequency.cofactors.ta.human,
    TAvsDN = log.ratio.confirmed.ta.vs.dn.human
)
my.volcano.plot(volcano_data.ta,title=paste0("Enrichment for p73 binding (TA confirmed by ",threshold.num.reads,"+ reads)"))

relevant.cofactors.ta <- enrichment.ta.human > 0.1 & frequency.cofactors.ta.human > 0.05
relevant.cofactors <- TRUE & relevant.cofactors.ta

volcano_data.dn <- data.frame(
    TF = sapply(strsplit(x=names(enrichment.dn.human),split="[_ ]",fixed=FALSE), function(x) x[1]),
    Enrichment = enrichment.dn.human,
    Frequency = frequency.cofactors.dn.human,
    TAvsDN = log.ratio.confirmed.ta.vs.dn.human
)
my.volcano.plot(volcano_data.dn,title=paste0("Enrichment for p73 binding (DN confirmed by ",threshold.num.reads,"+ reads)"))

relevant.cofactors.dn <- enrichment.dn.human > 0.1 & frequency.cofactors.dn.human > 0.05
relevant.cofactors <- relevant.cofactors | relevant.cofactors.dn

# p73 binding and Methylation enrichment plots

if (FALSE) {
    if (FALSE) {
        volcano_data.methylation.is.up.taOrDn <- data.frame(
            TF = sapply(strsplit(x=names(enrichment.taOrDn.methylation.is.up.vs.taOrDn.tp73.human),split="[_ ]",fixed=FALSE), function(x) x[1]),
            Enrichment = enrichment.taOrDn.methylation.is.up.vs.taOrDn.tp73.human,
            Frequency = frequency.cofactors.taOrDn.methylation.is.up.human,
            TAvsDN = log.ratio.confirmed.ta.vs.dn.methylation.is.up.human
        )
        my.volcano.plot(volcano_data.methylation.is.up.taOrDn,title=paste0("Enrichment for p73 binding & increased methylation (TA, DN or GFP confirmed by ",threshold.num.reads,"+ reads)"))
    }
    volcano_data.methylation.is.up.ta <- data.frame(
        TF = sapply(strsplit(x=names(enrichment.ta.methylation.is.up.vs.ta.tp73.human),split="[_ ]",fixed=FALSE), function(x) x[1]),
        Enrichment = enrichment.ta.methylation.is.up.vs.ta.tp73.human,
        Frequency = frequency.cofactors.ta.methylation.is.up.human,
        TAvsDN = log.ratio.confirmed.ta.vs.dn.methylation.is.up.human
    )
    my.volcano.plot(volcano_data.methylation.is.up.ta,title=paste0("Enrichment for p73 binding & increased methylation (TA confirmed by ",threshold.num.reads,"+ reads)"))

    relevant.cofactors.methylation.is.up.ta <- enrichment.ta.methylation.is.up.vs.ta.tp73.human > 0.1 & frequency.cofactors.ta.methylation.is.up.human > 0.05
    relevant.cofactors.methylation.is.up <- TRUE & relevant.cofactors.methylation.is.up.ta

    volcano_data.methylation.is.up.dn <- data.frame(
        TF = sapply(strsplit(x=names(enrichment.dn.methylation.is.up.vs.dn.tp73.human),split="[_ ]",fixed=FALSE), function(x) x[1]),
        Enrichment = enrichment.dn.methylation.is.up.vs.dn.tp73.human,
        Frequency = frequency.cofactors.dn.methylation.is.up.human,
        TAvsDN = log.ratio.confirmed.ta.vs.dn.methylation.is.up.human
    )
    my.volcano.plot(volcano_data.methylation.is.up.dn,title=paste0("Enrichment for p73 binding & increased methylation (DN confirmed by ",threshold.num.reads,"+ reads)"))

    relevant.cofactors.methylation.is.up.dn <- enrichment.dn.methylation.is.up.vs.dn.tp73.human > 0.1 & frequency.cofactors.dn.methylation.is.up.human > 0.05
    relevant.cofactors.methylation.is.up <- relevant.cofactors.methylation.is.up | relevant.cofactors.methylation.is.up.dn
}

dev.off()
cat("I: Volcano plot saved to '",pdf.volcano.filename,"'.\n",sep="")

write.xlsx(list(TA=volcano_data.ta[,-4],DN=volcano_data.dn[,-4]), file="volcano_data.xlsx", row.names=TRUE)


relevant.cofactors.selection.ta <- names(which(relevant.cofactors.methylation.is.up.ta & relevant.cofactors))
relevant.cofactors.selection.dn <- names(which(relevant.cofactors.methylation.is.up.dn & relevant.cofactors))
relevant.cofactors.selection <- names(which(relevant.cofactors.methylation.is.up & relevant.cofactors))


