
# 
# Preparation of Vulcano Plot
#

# Also defines
#  relevant.cofactors.selection.ta
#  relevant.cofactors.selection.dn
#  relevant.cofactors.selection

require(ggplot2)

threshold.num.reads <- 1

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
m.contexts.num.tfbs.total <- sum(sapply(m.contexts, nrow),  na.rm = TRUE)
frequency.cofactors.initial <- m.contexts.num.binary.sum.total/m.contexts.num.tfbs.total

# 4446155

#rm(m); m <- m.contexts[[22]]

count.cofactors.per.tfbs <- function(confirmation=NULL, TA.or.DN, cell.lines, threshold.num.reads=1 ) {

    if (! all (cell.lines %in% c("skmel29_2","skmel29_1","saos2"))) {
        stop("E: Invalid cell line specified. Only 'skmel29_2' or 'saos2' are supported.\n")
    }

    if (! all (cell.lines %in% c("TA","DN","any"))) {
        stop("E: Invalid variant specified. Only 'TA', 'DN', 'GFP' or 'any' are supported.\n")
    }

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
            prefix <- paste0("tp73","_",cell.line,"_")
            for (cell.line in cell.lines) {
                if (TA.or.DN %in% c("TA","DN","GFP")) {
                    v <- paste0(prefix, TA.or.DN)
                    m.context.extra.check <- m.context.extra.check &  as.vector(m[,..v,drop=T] >= threshold.num.reads )
                } else if (TA.or.DN == "any") {
                    v.ta <- paste0(prefix, "TA")
                    v.dn <- paste0(prefix, "DN")
                    v.gfp <- paste0(prefix, "GFP")
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.ta, drop=T] >= threshold.num.reads |
                        m[,..v.dn, drop=T] >= threshold.num.reads |
                        m[,..v.gfp,drop=T] >= threshold.num.reads
                    )
                } else if (TA.or.DN == "ignored") {
                    # doing nothing, as we want to check all rows
                } else {
                    stop("E: Invalid TA.or.DN method specified.\n")
                }
            }
        }

        if ("pos" %in% confirmation) {
            for (cell.line in cell.lines) {
                prefix <- paste0("pos","_",cell.line,"_")
                if (TA.or.DN %in% c("TA","DN","GFP")) {
                    v <- paste0(prefix, TA.or.DN)
                    m.context.extra.check <- m.context.extra.check & as.vector( m[,..v,drop=T] >= threshold.num.reads )
                } else if (TA.or.DN == "any") {
                    v.ta <- paste0(prefix, "TA")
                    v.dn <- paste0(prefix, "DN")
                    v.gfp <- paste0(prefix, "GFP")
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.ta, drop=T] >= threshold.num.reads |
                        m[,..v.dn, drop=T] >= threshold.num.reads |
                        m[,..v.gfp,drop=T] >= threshold.num.reads
                    )
                } else if (TA.or.DN == "ignored") {
                    # doing nothing, as we want to check all rows
                } else {
                    stop("E: Invalid TA.or.DN method specified.\n")
                }
            }
        }

        if ("pos-up" %in% confirmation) {
            for (cell.line in cell.lines) {
                prefix <- paste0("pos","_",cell.line,"_")
                v.ta <- paste0(prefix, "TA")
                v.dn <- paste0(prefix, "DN")
                v.gfp <- paste0(prefix, "GFP")
                if (TA.or.DN == "TA") {
                    m.context.extra.check <- m.context.extra.check & as.vector( m[,..v.ta,drop=T] >= m[,..v.gfp,drop=T] )
                } else if (TA.or.DN == "DN") {
                    m.context.extra.check <- m.context.extra.check & as.vector( m[,..v.dn,drop=T] >= m[,..v.gfp,drop=T] )
                } else if (TA.or.DN == "any") {
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.dn,drop=T] >= m[,..v.gfp,drop=T] |
                        m[,..v.ta,drop=T] >= m[,..v.gfp,drop=T]
                    )
                } else if (TA.or.DN == "ignored") {
                    # doing nothing, as we want to check all rows
                } else {
                    stop("E: Invalid TA.or.DN method specified.\n")
                }
            }
        }

        m <- m[m.context.extra.check, , drop = FALSE]
        m <- m[, ..cols.Score.index, drop = FALSE]
        colSums(m > 0, na.rm = TRUE) # count the number of TFBS that are confirmed by TA, GFP or DN
    })
}

#m.contexts.num.binary.confirmed.taOrDn.sum <- count.cofactors.per.tfbs(confirmation=c("tp73"), TA.or.DN="any", threshold.num.reads=threshold.num.reads)
#m.contexts.num.binary.confirmed.taOrDn.sum.total <- rowSums(m.contexts.num.binary.confirmed.taOrDn.sum, na.rm = TRUE)

count.tfbs <- function(confirmation=NULL, TA.or.DN, cell.lines, threshold.num.reads=1 ) {

    if (! all (cell.lines %in% c("skmel29_2","skmel29_1","saos2"))) {
        stop("E: Invalid cell line specified. Only 'skmel29_2' or 'saos2' are supported.\n")
    }

    if (! TA.or.DN %in% c("TA","DN","GFP","any","TAandDN","all","ignored")) {
        stop("E: Invalid variant specified. Only 'TA', 'DN', 'GFP', 'any', 'TAandDN', 'all' or 'ignored' are supported, found '",TA.or.DN,"'.\n")
    }

    sapply(m.contexts, function(m) {
        # m <- m.contexts[[1]]
        if (is.null(m)) {
            return(NA)
        }
        m.context.extra.check <- rep(TRUE, nrow(m))
        #cat("D: Initial sum(m.context.extra.check): ",sum(m.context.extra.check, na.rm = TRUE),"\n") 
        cols.Score.index <- grep("_NumInWindow$", colnames(m), value = FALSE)

        if ("promoter" %in% confirmation) {
            m.context.extra.check <- m.context.extra.check & m[, "InPromoter"]
        }

        if ("tp73" %in% confirmation) {
            for (cell.line in cell.lines) {
                prefix <- paste0("tp73","_",cell.line,"_")
                v.ta <- paste0(prefix, "TA")
                v.dn <- paste0(prefix, "DN")
                v.gfp <- paste0(prefix, "GFP")
                if (TA.or.DN %in% c("TA","DN","GFP")) {
                    v <- paste0(prefix, TA.or.DN)
                    m.context.extra.check <- m.context.extra.check & as.vector( m[,..v, drop=T] >= threshold.num.reads )
                } else if (TA.or.DN == "any") {
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.ta, drop=T] >= threshold.num.reads |
                        m[,..v.dn, drop=T] >= threshold.num.reads |
                        m[,..v.gfp, drop=T] >= threshold.num.reads
                    )
                } else if (TA.or.DN == "TAandDN") {
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.ta, drop=T] >= threshold.num.reads &
                        m[,..v.dn, drop=T] >= threshold.num.reads
                    )
                } else if (TA.or.DN == "all") {
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.ta, drop=T] >= threshold.num.reads &
                        m[,..v.dn, drop=T] >= threshold.num.reads &
                        m[,..v.gfp, drop=T] >= threshold.num.reads
                    )
                } else if (TA.or.DN == "ignored") {
                    # doing nothing, as we want to check all rows
                } else {
                    stop("E: Invalid TA.or.DN method specified, found '",ta.or.dn,"'.\n")
                }
            }
        }

        if ("pos" %in% confirmation) {
            for (cell.line in cell.lines) {
                prefix <- paste0("pos","_",cell.line,"_")
                v.ta <- paste0(prefix, "TA")
                v.dn <- paste0(prefix, "DN")
                v.gfp <- paste0(prefix, "GFP")
                if (TA.or.DN %in% c("TA","DN","GFP")) {
                    v <- paste0(prefix, TA.or.DN)
                    m.context.extra.check <- m.context.extra.check & as.vector( m[,..v, drop=T] >= threshold.num.reads )
                } else if (TA.or.DN == "any") {
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.ta, drop=T] >= threshold.num.reads |
                        m[,..v.dn, drop=T] >= threshold.num.reads |
                        m[,..v.gfp, drop=T] >= threshold.num.reads
                    )
                } else if (TA.or.DN == "TAandDN") {
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.ta, drop=T] >= threshold.num.reads &
                        m[,..v.dn, drop=T] >= threshold.num.reads
                    )
                } else if (TA.or.DN == "all") {
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.ta, drop=T] >= threshold.num.reads &
                        m[,..v.dn, drop=T] >= threshold.num.reads &
                        m[,..v.gfp, drop=T] >= threshold.num.reads
                    )
                } else if (TA.or.DN == "ignored") {
                    # doing nothing, as we want to check all rows
                } else {
                    stop("E: Invalid TA.or.DN method specified, found '",ta.or.dn,"'.\n")
                }
            }
        }

        if ("pos-up" %in% confirmation) {
            for (cell.line in cell.lines) {
                prefix <- paste0("pos","_",cell.line,"_")
                v.ta <- paste0(prefix, "TA")
                v.dn <- paste0(prefix, "DN")
                v.gfp <- paste0(prefix, "GFP")
                if (TA.or.DN == "TA") {
                    m.context.extra.check <- m.context.extra.check & as.vector( m[,..v.ta, drop=T] >= m[,..v.gfp, drop=T] )
                } else if (TA.or.DN == "DN") {
                    m.context.extra.check <- m.context.extra.check & as.vector( m[,..v.dn, drop=T] >= m[,..v.gfp, drop=T] )
                } else if (TA.or.DN == "any") {
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.dn, drop=T] >= m[,..v.gfp, drop=T] |
                        m[,..v.ta, drop=T] >= m[,..v.gfp, drop=T]
                    )
                } else if (TA.or.DN == "all" || TA.or.DN == "TAandDN") {
                    m.context.extra.check <- m.context.extra.check & as.vector(
                        m[,..v.dn, drop=T] >= m[,..v.gfp, drop=T] &
                        m[,..v.ta, drop=T] >= m[,..v.gfp, drop=T]
                    )
                } else if (TA.or.DN == "ignored") {
                    # doing nothing, as we want to check all rows
                } else {
                    stop("E: Invalid TA.or.DN method specified - only TA or DN valid for pos-up, found '",TA.or.DN,"'.\n")
                }
            }
        }

        sum(m.context.extra.check, na.rm = TRUE) # count the number of TFBS that are confirmed by TA, GFP or DN - not the number of cofactors
    })
}

# 
# count.tfbs.cache
#
#
count.tfbs.cache <- matrix(rep(NA,2*7*2),nrow=2, ncol=7*2)
dim(count.tfbs.cache) <- c(2,7,2)
dimnames(count.tfbs.cache) <- list(cell.line=c("skmel29_2","saos2"),isoform=c("ignored","TA","DN","GFP","any","TAandDN","all"),inPromoter=c("maybe","yes"))
for( cell.line in c("skmel29_2","saos2") ) {
    for( ta.or.dn in c("ignored","TA","DN","GFP","any","TAandDN","all")) {
        cat("D: ta.or.dn: ",ta.or.dn," cell.line: ",cell.line,"\n",sep="")
        count.tfbs.cache[cell.line,ta.or.dn,"maybe"] <- sum(count.tfbs(confirmation=c("tp73"), TA.or.DN=ta.or.dn, cell.lines=cell.line, threshold.num.reads=threshold.num.reads))
        count.tfbs.cache[cell.line,ta.or.dn,"yes"] <- sum(count.tfbs(confirmation=c("tp73","promoter"), TA.or.DN=ta.or.dn, cell.lines=cell.line, threshold.num.reads=threshold.num.reads))
    }
}

#m.contexts.num.tfbs.confirmed.taOrDn2 <- count.tfbs(confirmation=c("tp73"), TA.or.DN="any", threshold.num.reads=threshold.num.reads)
#m.contexts.num.tfbs.confirmed.taOrDn.total <- sum(m.contexts.num.tfbs.confirmed.taOrDn, na.rm = TRUE)
#frequency.cofactors.taOrDn <- m.contexts.num.binary.confirmed.taOrDn.sum.total/m.contexts.num.tfbs.confirmed.taOrDn.total
#names(frequency.cofactors.taOrDn) <- prettyIdentifierJaspar(names(frequency.cofactors.taOrDn))


rm(list=c("cell.line","m.context.num.binary.confirmed.sum.total","m.context.num.tfbs.confirmed.total","frequency.cofactors","frequency.cofactors.human",
   "log.ratio.confirmed.ta.vs.dn","log.ratio.confirmed.ta.vs.dn.human","enrichment,enrichment.human","a","ta.or.dn","TA.or.DN","enrichment.mean"))

volcano.plot.for.cell.line <- function(cell.line, threshold.num.reads=1, debug=TRUE,
                                       m.context.num.binary.confirmed.sum.total=m.context.num.binary.confirmed.sum.total,
                                       m.context.num.tfbs.confirmed.total=m.context.num.tfbs.confirmed.total) {

    # cell.line <- "skmel29_2"

    if (! cell.line %in% c("skmel29_2","skmel29_1","saos2")) {
        stop("E: Invalid cell line specified. Only 'skmel29_2', 'skmel29_1' or 'saos2' are supported.\n")
    }

    if (debug) {
        cat("D: Preparing data for Volcano plot of cell line '",cell.line,"' with threshold ",threshold.num.reads,".\n",sep="")
    }

    m.context.num.binary.confirmed.sum.total <- sapply(c("TA","DN","GFP"), function(ta.or.dn) {
        rowSums(count.cofactors.per.tfbs(confirmation=c("tp73"), TA.or.DN=ta.or.dn, cell.lines=cell.line, threshold.num.reads=threshold.num.reads), na.rm = TRUE)
    })

    print(apply(m.context.num.binary.confirmed.sum.total, 2, quantile, na.rm = TRUE))

    if (debug) {
        cat("D: 1\n")
    }

    m.context.num.tfbs.confirmed.total <- sapply(c("TA","DN","GFP"), function(ta.or.dn) {
        sum(count.tfbs(confirmation=c("tp73"), TA.or.DN=ta.or.dn, cell.lines=cell.line, threshold.num.reads=threshold.num.reads))
    })

    if (debug) {
        cat("D: 2\n")
    }

    frequency.cofactors <- sapply(c("TA","DN","GFP"), function(ta.or.dn,m.context.num.binary.confirmed.sum.tota=m.context.num.binary.confirmed.sum.total,m.context.num.tfbs.confirmed.total=m.context.num.tfbs.confirmed.total) {

        total.confirmed <- m.context.num.binary.confirmed.sum.total[,ta.or.dn]
        if (is.na(total.confirmed)) {
            stop(paste0("E: Values for 'm.context.num.binary.confirmed.sum.total[\"",ta.or.dn,"\"]' have not been computed."))
        }

        total.predicted.jaspar <- m.context.num.tfbs.confirmed.total[ta.or.dn]
        if(is.na(m.context.num.tfbs.confirmed.total)) {
            stop(paste0("E: Values for 'm.context.num.tfbs.confirmed.total[\"",ta.or.dn,"\"]' have not been computed."))
        }

        a <- total.confirmed / m.context.num.tfbs.confirmed.total
        names(a) <- prettyIdentifierJaspar(names(a))
        a
    })

    print(apply(frequency.cofactors, 2, quantile, na.rm = TRUE))

    if (debug) {
        cat("D: 3\n")
    }

    frequency.cofactors.human <- frequency.cofactors[is.human.jaspar.id(rownames(frequency.cofactors)),]

    if (debug) {
        cat("D: 4\n")
    }

    log.ratio.confirmed.ta.vs.dn <- log2(frequency.cofactors[,"TA"] / frequency.cofactors[,"DN"])
    log.ratio.confirmed.ta.vs.dn.human <- log.ratio.confirmed.ta.vs.dn[is.human.jaspar.id(names(log.ratio.confirmed.ta.vs.dn))]
 
    if (debug) {
        cat("D: 5\n")
    }

    cat("All names in sync: "); print(all(rownames(frequency.cofactors) == prettyIdentifierJaspar(names(frequency.cofactors.initial))))

    enrichment <- apply(frequency.cofactors,2,function(X) {
        log2(X / frequency.cofactors.initial)
    })

    enrichment.mean <- apply(enrichment, 2, mean, na.rm = TRUE)
    if (any(abs(enrichment.mean) > 0.50)) {
        print(apply(enrichment, 2, quantile,na.rm=T))
        stop("E: Enrichment mean values are too high, check the data or threshold.")
    }

    if (debug) {
        cat("D: 6\n")
    }

    enrichment.human <- enrichment[is.human.jaspar.id(rownames(enrichment)),]

    list(
        enrichment = enrichment,
        enrichment.human = enrichment.human,
        m.context.num.tfbs.confirmed.total = m.context.num.tfbs.confirmed.total,
        m.context.num.binary.confirmed.sum.total = m.context.num.binary.confirmed.sum.total,
        frequency.cofactors = frequency.cofactors,
        frequency.cofactors.human = frequency.cofactors.human,
        threshold.num.reads = threshold.num.reads,
        log.ratio.confirmed.ta.vs.dn = log.ratio.confirmed.ta.vs.dn,
        log.ratio.confirmed.ta.vs.dn.human = log.ratio.confirmed.ta.vs.dn.human,
        cell.line = cell.line
        )
}

if (FALSE) {
    m.contexts.num.binary.confirmed.ta.sum.total <- rowSums(count.cofactors.per.tfbs(confirmation=c("tp73"), TA.or.DN="TA", threshold.num.reads=threshold.num.reads), na.rm = TRUE)
    m.contexts.num.tfbs.confirmed.ta.total <- sum(count.tfbs(confirmation=c("tp73"), TA.or.DN="TA", threshold.num.reads=threshold.num.reads))
    # 948129
    frequency.cofactors.ta <- m.contexts.num.binary.confirmed.ta.sum.total/m.contexts.num.tfbs.confirmed.ta.total
    names(frequency.cofactors.ta) <- prettyIdentifierJaspar(names(frequency.cofactors.ta))

    m.contexts.num.binary.confirmed.dn.sum.total <- rowSums(count.cofactors.per.tfbs(confirmation=c("tp73"), TA.or.DN="DN", threshold.num.reads=threshold.num.reads), na.rm = TRUE)
    m.contexts.num.tfbs.confirmed.dn.total <- sum(count.tfbs(confirmation=c("tp73"), TA.or.DN="DN", threshold.num.reads=threshold.num.reads))
    frequency.cofactors.dn <- m.contexts.num.binary.confirmed.dn.sum.total/m.contexts.num.tfbs.confirmed.dn.total
    names(frequency.cofactors.dn) <- prettyIdentifierJaspar(names(frequency.cofactors.dn))

    log.ratio.confirmed.ta.vs.dn <- log2(frequency.cofactors.ta / frequency.cofactors.dn)
    names(log.ratio.confirmed.ta.vs.dn) <- prettyIdentifierJaspar(names(log.ratio.confirmed.ta.vs.dn))
    log.ratio.confirmed.ta.vs.dn.human <- log.ratio.confirmed.ta.vs.dn[is.human.jaspar.id(names(log.ratio.confirmed.ta.vs.dn))]

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
}


# Volcano-like scatter plot: enrichment.taOrDn (X) vs frequency.cofactors.taOrDn (Y)
require(ggplot2)
require(ggrepel)
#log.ratio.confirmed.ta.vs.dn.human.ranked <- rank(log.ratio.confirmed.ta.vs.dn.human, ties.method = "first")
#log.ratio.confirmed.ta.vs.dn.human.ranked <- log.ratio.confirmed.ta.vs.dn.human.ranked - mean(log.ratio.confirmed.ta.vs.dn.human.ranked)

my.volcano.plot <- function(data,
                            title="Volcano Plot: Enrichment vs Frequency",
                            xlab="log2 Enrichment (confirmed / initial)",
                            ylab="Frequency (confirmed / initial)",n.tfbs.prior=NA,n.tfbs.post=NA,
                            enrichment.min=NA, enrichment.max=NA) {

    cat("D: Plotting Volcano plot with title '",title,"'.\n",sep="")

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
        xlim(enrichment.min, enrichment.max) +
        labs(
            title=title,
            x=xlab,
            y=ylab
        ) +
        theme_minimal()
    
}

#enrichment.ta.human>0 & enrichment.dn.human<0


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

if (FALSE) {

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

}

d.skmel29_2 <- volcano.plot.for.cell.line(cell.line="skmel29_2", threshold.num.reads=threshold.num.reads, debug=TRUE)

pdf.volcano.filename <- paste("volcano_plot_enrichment_vs_frequency_treshold_",threshold.num.reads,"_",format(Sys.time(), "%Y%m%d"),"_skmel29.pdf",sep="")
pdf(pdf.volcano.filename, width=11, height=8)

    # SKMel29

    ta.or.dn <- "TA" # "TA" or "DN"

    enrichment.human.min <- min(d.skmel29_2$enrichment.human, na.rm=TRUE)
    enrichment.human.max <- max(d.skmel29_2$enrichment.human, na.rm=TRUE)

    
    if (FALSE) {
        ta.or.dn <- "any" # "TA" or "DN"
        volcano_data <- data.frame(
            TF = sapply(strsplit(x=names(d.skmel29_2$enrichment[,ta.or.dn]),split="[_ ]",fixed=FALSE), function(x) x[1]),
            Enrichment = d.skmel29_2$enrichment[,ta.or.dn],
            Frequency = d.skmel29_2$frequency.cofactors[,ta.or.dn],
            TAvsDN=d.skmel29_2$enrichment[,"TA"] - d.skmel29_2$enrichment[,"DN"]
        )
        my.volcano.plot(volcano_data,title=paste0("Enrichment for p73 binding in SKMel-29 (",ta.or.dn," confirmed by ",threshold.num.reads,"+ reads) - all of JASPAR"),
                        enrichment.min=enrichment.human.min,enrichment.max=enrichment.human.max)
    }
    volcano_data.ta <- volcano_data <- data.frame(
        TF = sapply(strsplit(x=names(d.skmel29_2$enrichment.human[,ta.or.dn]),split="[_ ]",fixed=FALSE), function(x) x[1]),
        Enrichment = d.skmel29_2$enrichment.human[,ta.or.dn],
        Frequency = d.skmel29_2$frequency.cofactors.human[,ta.or.dn],
        TAvsDN=d.skmel29_2$enrichment.human[,"TA"] - d.skmel29_2$enrichment.human[,"DN"]
    )
    my.volcano.plot(volcano_data,title=paste0("Enrichment for p73 binding in SKMel-29 (",ta.or.dn," confirmed by ",threshold.num.reads,"+ reads) - Human subset of JASPAR"),
                        enrichment.min=enrichment.human.min,enrichment.max=enrichment.human.max)

    ta.or.dn <- "DN" # "TA" or "DN"

    volcano_data.dn <- volcano_data <- data.frame(
        TF = sapply(strsplit(x=names(d.skmel29_2$enrichment.human[,ta.or.dn]),split="[_ ]",fixed=FALSE), function(x) x[1]),
        Enrichment = d.skmel29_2$enrichment.human[,ta.or.dn],
        Frequency = d.skmel29_2$frequency.cofactors.human[,ta.or.dn],
        TAvsDN=d.skmel29_2$enrichment.human[,"TA"] - d.skmel29_2$enrichment.human[,"DN"]
    )
    my.volcano.plot(volcano_data,title=paste0("Enrichment for p73 binding in SKMel-29 (",ta.or.dn," confirmed by ",threshold.num.reads,"+ reads) - Human subset of JASPAR"),
                        enrichment.min=enrichment.human.min,enrichment.max=enrichment.human.max)

    ta.or.dn <- "GFP" # "TA" or "DN"

    volcano_data.gfp <- volcano_data <- data.frame(
        TF = sapply(strsplit(x=names(d.skmel29_2$enrichment.human[,ta.or.dn]),split="[_ ]",fixed=FALSE), function(x) x[1]),
        Enrichment = d.skmel29_2$enrichment.human[,ta.or.dn],
        Frequency = d.skmel29_2$frequency.cofactors.human[,ta.or.dn],
        TAvsDN=d.skmel29_2$enrichment.human[,"TA"] - d.skmel29_2$enrichment.human[,"DN"]
    )
    my.volcano.plot(volcano_data,title=paste0("Enrichment for p73 binding in SKMel-29 (",ta.or.dn," confirmed by ",threshold.num.reads,"+ reads) - Human subset of JASPAR"),
                        enrichment.min=enrichment.human.min,enrichment.max=enrichment.human.max)

dev.off()
cat("I: Volcano plot saved to '",pdf.volcano.filename,"'.\n",sep="")

write.table(list(SKMel29.TA=volcano_data.ta[,-4],SKMel29.DN=volcano_data.dn[,-4],SKMel29.GFP=volcano_data.gfp[,-4]), file=paste0("volcano_data_skmel29_",format(Sys.time(), "%Y%m%d"),".tsv"), row.names=TRUE,na="",col.names=NA,quote=FALSE,sep="\t", dec=",")



d.saos2 <- volcano.plot.for.cell.line(cell.line="saos2", threshold.num.reads=threshold.num.reads, debug=TRUE)

pdf.volcano.filename <- paste("volcano_plot_enrichment_vs_frequency_treshold_",threshold.num.reads,"_",format(Sys.time(), "%Y%m%d"),"_saos2.pdf",sep="")
pdf(pdf.volcano.filename, width=11, height=8)

    # Saos2

    ta.or.dn <- "TA" # "TA" or "DN"

    if (FALSE) {
        volcano_data <- data.frame(
            TF = sapply(strsplit(x=names(d.saos2$enrichment[,ta.or.dn]),split="[_ ]",fixed=FALSE), function(x) x[1]),
            Enrichment = d.saos2$enrichment[,ta.or.dn],
            Frequency = d.saos2$frequency.cofactors[,ta.or.dn]
        )
        my.volcano.plot(volcano_data,title=paste0("Enrichment for p73 binding in Saos2 (",ta.or.dn," confirmed by ",threshold.num.reads,"+ reads) - all of JASPAR"))
    }

    volcano_data <- data.frame(
        TF = sapply(strsplit(x=names(d.saos2$enrichment.human[,ta.or.dn]),split="[_ ]",fixed=FALSE), function(x) x[1]),
        Enrichment = d.saos2$enrichment.human[,ta.or.dn],
        Frequency = d.saos2$frequency.cofactors.human[,ta.or.dn],
        TAvsDN=d.saos2$enrichment.human[,"TA"] - d.saos2$enrichment.human[,"DN"]
    )
    my.volcano.plot(volcano_data,title=paste0("Enrichment for p73 binding in Saos2 (",ta.or.dn," confirmed by ",threshold.num.reads,"+ reads) - Human subset of JASPAR"))

    ta.or.dn <- "DN" # "TA" or "DN"

    volcano_data <- data.frame(
        TF = sapply(strsplit(x=names(d.saos2$enrichment.human[,ta.or.dn]),split="[_ ]",fixed=FALSE), function(x) x[1]),
        Enrichment = d.saos2$enrichment.human[,ta.or.dn],
        Frequency = d.saos2$frequency.cofactors.human[,ta.or.dn],
        TAvsDN=d.saos2$enrichment.human[,"TA"] - d.saos2$enrichment.human[,"DN"]
    )
    my.volcano.plot(volcano_data,title=paste0("Enrichment for p73 binding in Saos2 (",ta.or.dn," confirmed by ",threshold.num.reads,"+ reads) - Human subset of JASPAR"))


dev.off()


volcano_data.ta.human <- data.frame(
    TF = sapply(strsplit(x=names(d.skmel29_2$enrichment.human),split="[_ ]",fixed=FALSE), function(x) x[1]),
    Enrichment = d.skmel29_2$enrichment.human[,"TA"],
    Frequency = d.skmel29_2$frequency.cofactors.human[,"TA"],
    TAvsDN = log.ratio.confirmed.ta.vs.dn.human
)
my.volcano.plot(volcano_data.ta,title=paste0("Enrichment for p73 binding (TA confirmed by ",threshold.num.reads,"+ reads)"))



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

write.xlsx(list(TA=volcano_data.ta[,-4],DN=volcano_data.dn[,-4]), file=paste0("volcano_data_",format(Sys.time(), "%Y%m%d"),".xlsx"), row.names=TRUE)
write.table(list(TA=volcano_data.ta[,-4],DN=volcano_data.dn[,-4]), file=paste0("volcano_data_",format(Sys.time(), "%Y%m%d"),".tsv"), row.names=TRUE,na="",col.names=NA,quote=FALSE,sep="\t", dec=",")


relevant.cofactors.selection.ta <- names(which(relevant.cofactors.methylation.is.up.ta & relevant.cofactors))
relevant.cofactors.selection.dn <- names(which(relevant.cofactors.methylation.is.up.dn & relevant.cofactors))
relevant.cofactors.selection <- names(which(relevant.cofactors.methylation.is.up & relevant.cofactors))


