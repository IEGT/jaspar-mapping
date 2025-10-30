#!/usr/bin/R

# Script to interpret TF binding sites for their association with CUT&RUN data.

options(width=180)

# Remove prior data on ratios or sums
rm(list=grep(ls(),pattern="^sum*",value=T))
rm(list=grep(ls(),pattern="^ratio.*",value=T))
rm(list=grep(ls(),pattern="^quantiles.*",value=T))
rm(list=grep(ls(),pattern="^mean.NumInWindow*",value=T))

require(ggplot2)
require(corrplot)
require(xlsx)
require(gplots)
require(ggplot2)
require(data.table)

jaspar.human <- read.delim("jaspar_homo.tsv",row.names=1,col.names=FALSE)

# Import series of functions and basic data
source("analyze_matrix_parameters.R")
source("analyze_matrix_function_lists.R")
source("analyze_matrix_function_distances.R")
source("analyze_matrix_function_promoters.R")
source("analyze_matrix_function_helper.R")

m.findings.filename <- "m.findings.RData"
m.context.filename <- "m.contexts.RData"

#meta.from.scratch <- TRUE
meta.from.scratch <- FALSE

if (!meta.from.scratch && file.exists(m.context.filename) && file.exists(m.findings.filename)) {

    cat("Loading '",m.context.filename,"'\n",sep="")
    load(file=m.context.filename, verbose=TRUE)
    cat("Loading '",m.findings.filename,"'\n",sep="")
    load(file=m.findings.filename,verbose=TRUE)

} else {

    m.findings <- vector("list", length(chromosomes))
    m.contexts <- vector("list", length(chromosomes))
    names(m.findings) <- names(m.contexts) <- chromosomes

    for(i in chromosomes) {
        cat("I: processing chromosome ",i,"...\n",sep="")
        m <- read.data.table.for.chromosome(i)

        if (is.null(m)) {
            cat("E: No data for chromosome ",i," - skipping\n",sep="")
            next
        }

        if ("Feature" %in% colnames(m)) {
            cat("W: Found 'Feature' for chromosome ",i," - renaming to 'Name'\n",sep="")
            m.colnames.which <- which(colnames(m)=="Feature")
            colnames(m)[m.colnames.which] <- "Name"
        }

        # Store context
        m.contexts[[i]] <- m

        m <- create.lists.for.chromosome(m)
        m.findings[[i]] <- attributes(m)

        # Clean up
        rm(m)
        gc()
    }

    # Save the findings for each chromosome
    cat("I: saving p73 contexts for each chromosome...\n")
    save(m.contexts, file="m.contexts.RData")
    cat("I: ... saved contexts, now saving findings...\n")
    save(m.findings, file="m.findings.RData")
    cat("I: ... saved findings successfully.\n")

}

gc()

expressionData.dir <- "."
expressionData.comparison.filename <- "SkMel29_GFP_TAa_DNb_2x3x2_ohne_Filter_20.05.2025.tsv"

num.skipped.because.of.gene.name.ambiguity <- NA
num.skipped.because.of.gene.name.unknown <- NA

if (!meta.from.scratch && file.exists("combined.expression.data.RData")) {

    cat("I: Loading existing combined expression data from 'combined.expression.data.RData'...\n")
    load("combined.expression.data.RData", verbose=TRUE)
    if (!exists("combined.expression.data")) {
        stop("E: 'combined.expression.data' not found in 'combined.expression.data.RData'.")
    }
    cat("I: Found ",nrow(combined.expression.data)," rows in combined expression data.\n",sep="")
    print(head(combined.expression.data))


    load("max.binding.for.gene.RData", verbose=TRUE)
    # Check if max.binding.for.gene is available
    if (exists("max.binding.for.gene")) {
        cat("I: Found max.binding.for.gene with ",nrow(max.binding.for.gene)," rows.\n",sep="")
        print(head(max.binding.for.gene))
    } else {
        cat("W: max.binding.for.gene not found in 'combined.expression.data.RData'.\n")
    }

} else {

    require(openxlsx) # at the end to write the results to an Excel file

    num.skipped.because.of.gene.name.ambiguity <- 0
    num.skipped.because.of.gene.name.unknown <- 0

    # Load the expression data
    expressionData <- fread(file.path(expressionData.dir,expressionData.comparison.filename), sep="\t", header=TRUE, stringsAsFactors=FALSE)
    # Caveat: Expression data have duplicated gene symbols
    expressionData.symbols <- expressionData$"Gene Symbol"

    # Iteratate over genes in expression data to retrieve promoter locations and find max binding
    list.of.all.promoters <- promoterBedTables[["all.promoter.sorted.cutandrun.tp73bs.bed"]]
    if (is.null(list.of.all.promoters)) {
        stop("E: 'all.promoter.cutandrun.tp73bs.bed' not found in promoterBedTables.")
    }
    list.of.all.utr <- utrBedTables[["all.utr.sorted.cutandrun.tp73bs.bed"]]
    if (is.null(list.of.all.utr)) {
        stop("E: 'all.utr.cutandrun.tp73bs.bed' not found in utrBedTables.")
    }
    # Check if the promoter and UTR data have the same number of columns
    if (ncol(list.of.all.promoters) != ncol(list.of.all.utr)) {
        stop("E: 'all.promoter.cutandrun.tp73bs.bed' and 'all.utr.cutandrun.tp73bs.bed' have different number of columns.")
    }
    if (nrow(list.of.all.promoters) != nrow(list.of.all.utr)) {
        stop("E: 'all.promoter.cutandrun.tp73bs.bed' and 'all.utr.cutandrun.tp73bs.bed' have different number of rows")
    }

    list.sum <- list.of.all.utr[,-(1:6)] + list.of.all.promoters[,-(1:6)]
    colnames(list.of.all.promoters) <- paste("promoter", colnames(list.of.all.promoters), sep=".")
    colnames(list.of.all.utr) <- paste("utr", colnames(list.of.all.utr), sep=".")
    colnames(list.sum) <- paste("sum", colnames(list.sum), sep=".")

    list.of.all.promoters.plus.utr <- cbind(list.of.all.promoters[,1],
                                            apply(cbind(list.of.all.promoters[,2],list.of.all.utr[,2]), 1, min),
                                            apply(cbind(list.of.all.promoters[,3],list.of.all.utr[,3]), 1, max),
                                            list.of.all.promoters[,4:6],
                                            list.sum,
                                            list.of.all.promoters[,-(1:6)],
                                            list.of.all.utr[,-(1:6)])
    colnames(list.of.all.promoters.plus.utr)[1:6] <- c("Chr","From","To","Gene","Score","Strand")
    #write.xlsx(list.of.all.promoters.plus.utr, file="all.promoter.and.utr.cutandrun.tp73bs.xlsx", row.names=FALSE)
    write.table(list.of.all.promoters.plus.utr, file="all.genes.all.promoters.and.utr.cutandrun.tp73bs.tsv", row.names=FALSE,col.names=TRUE, sep="\t", quote=FALSE,na="")

    list.of.all.promoters.plus.utr.aggregated = aggregate(list.of.all.promoters.plus.utr[,-c(1:6)], by=list("Gene"=list.of.all.promoters.plus.utr$Gene), FUN = mean, na.rm = TRUE)
    rownames(list.of.all.promoters.plus.utr.aggregated) <- list.of.all.promoters.plus.utr.aggregated$Gene
    list.of.all.promoters.plus.utr.aggregated[,-1] <- round(list.of.all.promoters.plus.utr.aggregated[,-1], digits=2)
    write.table(list.of.all.promoters.plus.utr.aggregated, file="all.genes.all.promoters.and.utr.cutandrun.tp73bs.aggregated.tsv", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE,na="")

    max.binding.for.gene <- as.data.frame(matrix(NA, nrow=length(expressionData.symbols), ncol=1+17+3+1))
    colnames(max.binding.for.gene) <- c(colnames(m.contexts[[1]])[1:18],"DN Max Methylation 500 + 500", "GFP Max Methylation 500 + 500", "TA Max Methylation 500 + 500", "PromoterOfWhichGene")
    colnames(max.binding.for.gene)[4:6] <- c("TF","TF.Score","TF.Strand")

    debug <- FALSE
    columns.with.methylation.of.interest <- c("pos_skmel29_2_DN","pos_skmel29_2_GFP","pos_skmel29_2_TA")
    columns.with.methylation.of.interest.promoters <- paste("promoter", columns.with.methylation.of.interest,sep=".")
    columns.with.methylation.of.interest.utr <- paste("utr", columns.with.methylation.of.interest,sep=".")


    #ed.rowno <- 57
    for(ed.rowno in 1:nrow(expressionData)) {

        #if (ed.rowno>500) {
        #    cat("I: Stopping after 500 genes processed.\n")
        #    break
        #}

        gene <- expressionData.symbols[ed.rowno]
        cat(ed.rowno,": ",gene,"\n",sep="")

        # Check if promoter region of gene is known in promoters
        gene.in.promoter.rows <- (gene == list.of.all.promoters$"promoter.Gene")

        if (!any(gene.in.promoter.rows)) {
            cat("E: Gene ",gene," not found in promoter data\n",sep="")
            max.binding.for.gene[ed.rowno,] <- c(NA,rep(NA,17+3),gene)
            num.skipped.because.of.gene.name.unknown <- num.skipped.because.of.gene.name.unknown + 1
            next
        }

        # Get the promoter location
        gene.promoter.location.promoters <- list.of.all.promoters[gene.in.promoter.rows,c("promoter.Chr","promoter.From","promoter.To","promoter.Gene","promoter.Strand"),drop=F]
        colnames(gene.promoter.location.promoters) <- c("Chr","From","To","Gene","Strand")
        gene.promoter.location.utr <- list.of.all.utr[gene.in.promoter.rows,c("utr.Chr","utr.From","utr.To","utr.Gene","utr.Strand"),drop=F]
        colnames(gene.promoter.location.utr) <- c("Chr","From","To","Gene","Strand")

        stopifnot(all(gene.promoter.location.promoters$promoter.Chr == gene.promoter.location.utr$utr.Chr))

        if (debug) print(gene.promoter.location)

        chr <- unique(gene.promoter.location.promoters$Chr)
        stopifnot(!is.null(chr))

        if (length(chr) != 1) {
            cat("E: Found ",length(chr)," chromosomes for gene ",gene," - skipping\n",sep="")
            max.binding.for.gene[ed.rowno,] <- c(paste(chr,collapse="+",sep=""),rep(NA,17+3),gene)
            num.skipped.because.of.gene.name.ambiguity <- num.skipped.because.of.gene.name.ambiguity + 1
            next
        }

        if (is.null(m.contexts[[chr]])) {
            cat("W: No data for chromosome ",chr," - skipping\n",sep="")
            max.binding.for.gene[ed.rowno,] <- c(NA,rep(NA,17+3),gene)
            num.skipped.because.of.gene.name.unknown <- num.skipped.because.of.gene.name.unknown + 1
            next
        }

        cat("Chromosome: ",chr,"\n",sep="")
        overlap.index.of.m.chr.promoters <- checkBedOverlaps(gene.promoter.location.promoters[,1:3],m.contexts[[chr]][,1:3])
        overlap.index.of.m.chr.utr <- checkBedOverlaps(gene.promoter.location.utr[,1:3],m.contexts[[chr]][,1:3])


        if (0 < length(overlap.index.of.m.chr.promoters)) {

            cat("I: Found ",length(overlap.index.of.m.chr.promoters)," overlaps for gene ",gene," in chromosome ",chr,"\n",sep="")
            # Get the promoter location
            gene.promoter.location.promoters <- m.contexts[[chr]][overlap.index.of.m.chr.promoters,1:18]
            if (debug) print(gene.promoter.location.promoters)
            val <- gene.promoter.location.promoters$"tp73_skmel29_2_DN" + gene.promoter.location.promoters$"tp73_skmel29_2_TA"
            stopifnot(!is.null(val))
            val.max <- max(val)
            which.max.val <- which(val==val.max)
            select.line <- which.max.val
            if (length(select.line) > 1) {
                val.pos <- gene.promoter.location.promoters$"pos_skmel29_2_DN" + gene.promoter.location.promoters$"pos_skmel29_2_TA"
                stopifnot(!is.null(val.pos))
                val.pos[val != val.max] <- NA
                val.pos.max <- max(val.pos, na.rm=TRUE)
                select.line <- which(val.pos==val.pos.max)
                if (length(select.line) > 1) {
                    val.score <- gene.promoter.location.promoters$Score
                    val.score[val.pos != val.pos.max] <- NA
                    val.score[val != val.max] <- NA
                    val.score.max <- max(val.score, na.rm=TRUE)
                    select.line <- which.max(val.score)
                }
            }

            cat("I: Found max binding for gene ",gene," in chromosome ",chr," at line ",select.line," with value ",val.max,".\n",sep="")
            methylation.data.of.interest.promoter <- list.of.all.promoters[gene.in.promoter.rows,,drop=F][select.line,..columns.with.methylation.of.interest.promoters,drop=F]
            methylation.data.of.interest.utr <- list.of.all.utr[gene.in.promoter.rows,,drop=F][select.line,..columns.with.methylation.of.interest.utr,drop=F]
            methylation.data.of.interest.sum <- methylation.data.of.interest.promoter + methylation.data.of.interest.utr
            max.binding.for.gene[ed.rowno,] <- c(gene.promoter.location.promoters[select.line,],methylation.data.of.interest.sum,PromotorOfWhichGene=gene)
            #break;
        } else {
            if (sum(list.of.all.promoters[gene.in.promoter.rows,"promoter.num.tfbs"])>0) {
                stop("E: Discrepancy - gene ",gene," has promoter regions in all.promoter.sorted.cutandrun.tp73bs.bed but no overlaps found in m.contexts for chromosome ",chr,".\n",sep="")
            }
            mean.methylation.for.gene.promoter <- colSums(list.of.all.promoters[gene.in.promoter.rows,..columns.with.methylation.of.interest.promoters,drop=F], na.rm=TRUE)/length(gene.in.promoter.rows)
            names(mean.methylation.for.gene.promoter) <- paste0("mean.500bp.promoter.",names(mean.methylation.for.gene.promoter))
            mean.methylation.for.gene.utr <- colSums(list.of.all.promoters[gene.in.promoter.rows,..columns.with.methylation.of.interest.promoters,drop=F], na.rm=TRUE)/length(gene.in.promoter.rows)
            names(mean.methylation.for.gene.utr) <- paste0("mean.500bp.utr.",names(mean.methylation.for.gene.utr))
            mean.methylation.for.gene.joint <- mean.methylation.for.gene.utr+mean.methylation.for.gene.promoter
            names(mean.methylation.for.gene.joint) <- paste0("mean.1000bp.",names(mean.methylation.for.gene.utr))
            max.binding.for.gene[ed.rowno,] <- c(chr,rep(NA,17),round(mean.methylation.for.gene.joint,2),"PromoterOfWhichGene"=gene)
            cat("E: No overlaps found for gene ",gene," in chromosome ",chr,"\n",sep="")
        }
    }
    rm(ed.rowno)

    save(max.binding.for.gene, num.skipped.because.of.gene.name.ambiguity, num.skipped.because.of.gene.name.unknown,
         file="max.binding.for.gene.RData")

    # Check the dimensions of the max.binding.for.gene
    if(!all(max.binding.for.gene$"Gene" == expressionData.symbols)) {
        cat("E: Found ",sum(max.binding.for.gene$"Gene" != expressionData.symbols)," mismatches in gene names\n",sep="")
        print(max.binding.for.gene[max.binding.for.gene$"Gene" != expressionData.symbols,])
    } else {
        cat("I: Found all ",length(expressionData.symbols)," genes in max.binding.for.gene\n",sep="")
    }

    combined.expression.data <- cbind(max.binding.for.gene,expressionData)
    save(combined.expression.data, file="combined.expression.data.RData")
    require(xlsx)
    xlsx::write.xlsx(combined.expression.data, file="combined.expression.data.xlsx", row.names=FALSE)
    write.table(combined.expression.data, file="combined.expression.data.tsv", sep="\t", row.names=FALSE, col.names=TRUE, na="", quote=FALSE, dec=",")
    cat("I: Saved combined expression data to 'combined.expression.data.RData' and 'combined.expression.data.xlsx'.\n")
}

gc()

# Classify all TFBS in m.contexts for being in promoter regions
cat("I: Classifying TFBS in promoter regions...\n")

if (!meta.from.scratch && file.exists("tfbs.in.promoter.list.RData")) {

    cat("I: Loading existing TFBS in promoter regions from 'tfbs.in.promoter.list.RData'...\n")
    load("tfbs.in.promoter.list.RData")
    if (!exists("tfbs.in.promoter.list")) {
        stop("E: 'tfbs.in.promoter.list' not found in 'tfbs.in.promoter.list.RData'.")
    }
    cat("I: Found ",length(tfbs.in.promoter.list)," chromosomes in tfbs.in.promoter.list.\n",sep="")

} else {

    tfbs.in.promoter.list <- lapply(names(m.contexts), function(i) {
        cat("I: processing chromosome ",i,"...\n",sep="")
        if (is.null(m.contexts[[i]])) {
            cat("E: No data for chromosome ",i," - skipping\n",sep="")
            return(NULL)
        }
        # Process overlaps

        m.promoters.index <- lapply(promoterBedTables, function(x) {
            checkBedOverlaps(x,m.contexts[[i]])
        })

        tfbs.is.in.promoter.list<-lapply(names(m.promoters.index),function(x) {
            cat("Analyzing TFBS in promoter regions for ",x," in chromosome ",i,"...\n",sep="")
            tfbs.in.promoter <- rep(FALSE, nrow(m.contexts[[i]]))
            for(j in m.promoters.index[[x]]) {
                if (is.na(j) || j < 1 || j > nrow(m.contexts[[i]])) {
                    stop("E: Invalid index ",j," for chromosome ",i," - skipping\n",sep="")
                }
                tfbs.in.promoter[j] <- TRUE
            }
            cat("I: Found ",sum(tfbs.in.promoter)," TFBS in promoter regions for ",x," in chromosome ",i,".\n",sep="")
            return(tfbs.in.promoter)
        })
        names(tfbs.is.in.promoter.list) <- names(m.promoters.index)

        cat("I: Found ",sum(m.contexts[[i]]$InPromoter)," TFBS in promoter regions for chromosome ",i,".\n",sep="")
        return(tfbs.is.in.promoter.list)
    })
    names(tfbs.in.promoter.list) <- names(m.contexts)

    save(tfbs.in.promoter.list, file="tfbs.in.promoter.list.RData")
}

gc()

# Assigning the promoter classification to all TFBS in m.contexts

for (i in names(m.contexts)) {
    cat("I: processing chromosome ",i,"...\n",sep="")
    if (is.null(m.contexts[[i]])) {
        cat("E: No data for chromosome ",i," - skipping\n",sep="")
        next
    }
    if (!"InPromoter" %in% colnames(m.contexts[[i]])) {
        m.contexts[[i]]$InPromoter <- NA
    }
    if (is.null(tfbs.in.promoter.list[[i]])) {
        cat("E: No TFBS in promoter regions found for chromosome ",i," - skipping\n",sep="")
        next
    }
    if(is.null(tfbs.in.promoter.list[[i]][["all.promoter.bed"]])) {
        cat("E: List of all promoter regions not found for chromosome ",i," - skipping\n",sep="")
        next
    }
    stopifnot(nrow(m.contexts[[i]]) == length(tfbs.in.promoter.list[[i]][["all.promoter.bed"]]))
    m.contexts[[i]]$InPromoter <- tfbs.in.promoter.list[[i]][["all.promoter.bed"]]
}

gc()

m.contexts.annotated.filename.RData <- "m.contexts.with.promoter.info.RData"
cat("I: Saving annotated contexts with promoter info to '",m.contexts.annotated.filename.RData,"'...\n",sep="")
save(m.contexts, file=m.contexts.annotated.filename.RData)
cat("I: ... saved annotated contexts successfully.\n")
m.contexts.annotated.dirname.tsv <-  gsub(x=m.contexts.annotated.filename.RData,"RData","tsv")
dir.create(m.contexts.annotated.dirname.tsv, showWarnings=FALSE)
for(chr in names(m.contexts)) {
    f <- file.path(m.contexts.annotated.dirname.tsv, paste0("Chr=",chr,".tsv.gz"))
    cat("I: Writing annotated context for chromosome ",chr," to TSV at ",f,"...\n",sep="")
    write.table(m.contexts[[chr]], file=gzfile(f), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE, na="")
}
gc()



source("analyze_matrix_function_retrieve_context.R")  # Needs Context with InPromoter column



# Series of analyses across chromosomes on absolute and relative abundances of TFBS with/without CUT&RUN data confirmation

analyze.context.for.chromosome <- function(m.context.i, chromosome) {
    cat("I: Analyzing context for chromosome ",chromosome,"...\n",sep="")

    # Check if m.context is NULL
    if (is.null(m.context.i)) {
        cat("E: No data for chromosome ",chromosome," - skipping\n",sep="")
        return(NULL)
    }

    # Calculate the number of TFBS in each context
    num.TAa <- sum(m.context.i$"tp73_skmel29_2_TA", na.rm=TRUE)
    num.DNb <- sum(m.context.i$"tp73_skmel29_2_DN", na.rm=TRUE)
    num.GFP <- sum(m.context.i$"tp73_skmel29_2_GFP", na.rm=TRUE)

    # Calculate the ratio of TAa to DNb
    ratio.TAa.to.DNb <- num.TAa / num.DNb

    # Calculate the quantiles for TAa and DNb
    quantile.TAa <- quantile(m.context.i$"tp73_skmel29_2_TA", probs=c(0, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm=TRUE)
    quantile.DNb <- quantile(m.context.i$"tp73_skmel29_2_DN", probs=c(0, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm=TRUE)
    quantile.GFP <- quantile(m.context.i$"tp73_skmel29_2_GFP", probs=c(0, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm=TRUE)

    quantile.TAa.promoter <- quantile(m.context.i$"tp73_skmel29_2_TA"[m.context.i$InPromoter], probs=c(0, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm=TRUE)
    quantile.DNb.promoter <- quantile(m.context.i$"tp73_skmel29_2_DN"[m.context.i$InPromoter], probs=c(0, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm=TRUE)
    quantile.GFP.promoter <- quantile(m.context.i$"tp73_skmel29_2_GFP"[m.context.i$InPromoter], probs=c(0, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm=TRUE)

    p73.bs.confirmed.TAa.percent <- sum(m.context.i$"tp73_skmel29_2_TA" > 0) / nrow(m.context.i) * 100
    p73.bs.confirmed.DNb.percent <- sum(m.context.i$"tp73_skmel29_2_DN" > 0) / nrow(m.context.i) * 100
    p73.bs.confirmed.GFP.percent <- sum(m.context.i$"tp73_skmel29_2_GFP" > 0) / nrow(m.context.i) * 100

    p73.bs.confirmed.TAa.promoter.percent <- sum(m.context.i$"tp73_skmel29_2_TA"[m.context.i$InPromoter] > 0) / sum(m.context.i$InPromoter) * 100
    p73.bs.confirmed.DNb.promoter.percent <- sum(m.context.i$"tp73_skmel29_2_DN"[m.context.i$InPromoter] > 0) / sum(m.context.i$InPromoter) * 100
    p73.bs.confirmed.GFP.promoter.percent <- sum(m.context.i$"tp73_skmel29_2_GFP"[m.context.i$InPromoter] > 0) / sum(m.context.i$InPromoter) * 100

    p73.bs.with.methylation.TAa.percent <- sum(m.context.i$"pos_skmel29_2_TA" > 0) / nrow(m.context.i) * 100
    p73.bs.with.methylation.DNb.percent <- sum(m.context.i$"pos_skmel29_2_DN" > 0) / nrow(m.context.i) * 100
    p73.bs.with.methylation.GFP.percent <- sum(m.context.i$"pos_skmel29_2_GFP" > 0) / nrow(m.context.i) * 100

    p73.bs.with.methylation.TAa.promoter.percent <- sum(m.context.i$"pos_skmel29_2_TA"[m.context.i$InPromoter] > 0) / sum(m.context.i$InPromoter) * 100
    p73.bs.with.methylation.DNb.promoter.percent <- sum(m.context.i$"pos_skmel29_2_DN"[m.context.i$InPromoter] > 0) / sum(m.context.i$InPromoter) * 100
    p73.bs.with.methylation.GFP.promoter.percent <- sum(m.context.i$"pos_skmel29_2_GFP"[m.context.i$InPromoter] > 0) / sum(m.context.i$InPromoter) * 100

    p73.bs.confirmed.TAa.with.methylation.TAa.percent <- sum(m.context.i$"tp73_skmel29_2_TA" > 0 & m.context.i$"pos_skmel29_2_TA" > 0) / nrow(m.context.i) * 100
    p73.bs.confirmed.DNb.with.methylation.DNb.percent <- sum(m.context.i$"tp73_skmel29_2_DN" > 0 & m.context.i$"pos_skmel29_2_DN" > 0) / nrow(m.context.i) * 100
    p73.bs.confirmed.GFP.with.methylation.GFP.percent <- sum(m.context.i$"tp73_skmel29_2_GFP" > 0 & m.context.i$"pos_skmel29_2_GFP" > 0) / nrow(m.context.i) * 100

    p73.bs.confirmed.TAa.with.methylation.TAa.promoter.percent <- sum(m.context.i$"tp73_skmel29_2_TA" > 0 & m.context.i$"pos_skmel29_2_TA">0 &m.context.i$InPromoter > 0) / sum(m.context.i$InPromoter) * 100
    p73.bs.confirmed.DNb.with.methylation.DNb.promoter.percent <- sum(m.context.i$"tp73_skmel29_2_DN" > 0 & m.context.i$"pos_skmel29_2_DN">0 & m.context.i$InPromoter > 0) / sum(m.context.i$InPromoter) * 100
    p73.bs.confirmed.GFP.with.methylation.GFP.promoter.percent <- sum(m.context.i$"tp73_skmel29_2_GFP" > 0 & m.context.i$"pos_skmel29_2_GFP">0 & m.context.i$InPromoter > 0) / sum(m.context.i$InPromoter) * 100

    p73.bs.with.methylation.TAa.if.confirmed.TAa.promoter.percent <-
       sum(m.context.i$"pos_skmel29_2_TA"[m.context.i$"tp73_skmel29_2_TA">0 & m.context.i$InPromoter]>0 ) / sum(m.context.i$InPromoter & m.context.i$"tp73_skmel29_2_TA">0) * 100
    p73.bs.with.methylation.DNb.if.confirmed.DNb.promoter.percent <-
       sum(m.context.i$"pos_skmel29_2_DN"[m.context.i$"tp73_skmel29_2_DN">0 & m.context.i$InPromoter]>0 ) / sum(m.context.i$InPromoter & m.context.i$"tp73_skmel29_2_DN">0) * 100
    p73.bs.with.methylation.GFP.if.confirmed.GFP.promoter.percent <-
       sum(m.context.i$"pos_skmel29_2_GFP"[m.context.i$"tp73_skmel29_2_GFP">0 & m.context.i$InPromoter]>0 ) / sum(m.context.i$InPromoter & m.context.i$"tp73_skmel29_2_GFP">0) * 100

    # Determine the number of TFBS that are confirmed by TA but not by DN
    p73.bs.with.confirmation.TAa.only.percent <- sum(m.context.i$"tp73_skmel29_2_TA" > 0 & m.context.i$"tp73_skmel29_2_DN" == 0)/nrow(m.context.i) * 100
    p73.bs.with.confirmation.DNb.only.percent <- sum(m.context.i$"tp73_skmel29_2_DN" > 0 & m.context.i$"tp73_skmel29_2_TA" == 0)/nrow(m.context.i) * 100
    p73.bs.with.confirmation.GFP.only.percent <- sum(m.context.i$"tp73_skmel29_2_TA" == 0 & m.context.i$"tp73_skmel29_2_DN" == 0 & m.context.i$"tp73_skmel29_2_GFP" > 0)/nrow(m.context.i) * 100
    p73.bs.with.confirmation.TAa.only.promoter.percent <-
       sum(m.context.i$"tp73_skmel29_2_TA" > 0 & m.context.i$"tp73_skmel29_2_DN" == 0 & m.context.i$InPromoter)/sum(m.context.i$InPromoter) * 100
    p73.bs.with.confirmation.DNb.only.promoter.percent <-
       sum(m.context.i$"tp73_skmel29_2_DN" > 0 & m.context.i$"tp73_skmel29_2_TA" == 0 & m.context.i$InPromoter)/sum(m.context.i$InPromoter) * 100
    p73.bs.with.confirmation.GFP.only.promoter.percent <-
       sum(m.context.i$"tp73_skmel29_2_TA" == 0 & m.context.i$"tp73_skmel29_2_DN" == 0 & m.context.i$"tp73_skmel29_2_GFP">0 & m.context.i$InPromoter)/sum(m.context.i$InPromoter) * 100

    p73.bs.with.methylation.TAa.only.confirmed.TAa.promoter.percent <-
       sum(m.context.i$"tp73_skmel29_2_TA">0 & m.context.i$"pos_skmel29_2_TA" > 0 & m.context.i$"pos_skmel29_2_DN" & m.context.i$InPromoter)/sum(m.context.i$InPromoter) * 100
    p73.bs.with.methylation.DNb.only.confirmed.DNb.promoter.percent <-
       sum(m.context.i$"tp73_skmel29_2_DN">0 & m.context.i$"pos_skmel29_2_DN" > 0 & m.context.i$"pos_skmel29_2_TA" & m.context.i$InPromoter)/sum(m.context.i$InPromoter) * 100
    p73.bs.with.methylation.GFP.only.confirmed.GFP.promoter.percent <-
       sum(m.context.i$"tp73_skmel29_2_GFP">0 & m.context.i$"pos_skmel29_2_TA" == 0 & m.context.i$"pos_skmel29_2_DN" == 0 & m.context.i$"pos_skmel29_2_GFP" > 0 & m.context.i$InPromoter)/sum(m.context.i$InPromoter) * 100

       #Num.TAa = num.TAa,
        #Num.DNb = num.DNb,
        #Num.GFP = num.GFP,

    # Create a summary table
    summary.table <- data.frame(
        Chromosome = chromosome,
        Num.p73bs = nrow(m.context.i),
        Num.p73bs.InPromoter.percent = sum(m.context.i$InPromoter, na.rm=TRUE)/nrow(m.context.i) * 100,
        Ratio.TAa.to.DNb = round(ratio.TAa.to.DNb,2),
        p73.bs.confirmed.TAa.percent = round(p73.bs.confirmed.TAa.percent,2),
        p73.bs.confirmed.DNb.percent = round(p73.bs.confirmed.DNb.percent,2),
        p73.bs.confirmed.GFP.percent = round(p73.bs.confirmed.GFP.percent,2),
        p73.bs.confirmed.TAa.promoter.percent = round(p73.bs.confirmed.TAa.promoter.percent,2),
        p73.bs.confirmed.DNb.promoter.percent = round(p73.bs.confirmed.DNb.promoter.percent,2),
        p73.bs.confirmed.GFP.promoter.percent = round(p73.bs.confirmed.GFP.promoter.percent,2),
        p73.bs.with.methylation.TAa.percent = round(p73.bs.with.methylation.TAa.percent,2),
        p73.bs.with.methylation.DNb.percent = round(p73.bs.with.methylation.DNb.percent,2),
        p73.bs.with.methylation.GFP.percent = round(p73.bs.with.methylation.GFP.percent,2),
        p73.bs.with.methylation.TAa.promoter.percent = round(p73.bs.with.methylation.TAa.promoter.percent,2),
        p73.bs.with.methylation.DNb.promoter.percent = round(p73.bs.with.methylation.DNb.promoter.percent,2),
        p73.bs.with.methylation.GFP.promoter.percent = round(p73.bs.with.methylation.GFP.promoter.percent,2),
        p73.bs.confirmed.TAa.with.methylation.TAa.percent = round(p73.bs.confirmed.TAa.with.methylation.TAa.percent,2),
        p73.bs.confirmed.DNb.with.methylation.DNb.percent = round(p73.bs.confirmed.DNb.with.methylation.DNb.percent,2),
        p73.bs.confirmed.GFP.with.methylation.GFP.percent = round(p73.bs.confirmed.GFP.with.methylation.GFP.percent,2),
        p73.bs.confirmed.TAa.with.methylation.TAa.promoter.percent = round(p73.bs.confirmed.TAa.with.methylation.TAa.promoter.percent,2),
        p73.bs.confirmed.DNb.with.methylation.DNb.promoter.percent = round(p73.bs.confirmed.DNb.with.methylation.DNb.promoter.percent,2),
        p73.bs.confirmed.GFP.with.methylation.GFP.promoter.percent = round(p73.bs.confirmed.GFP.with.methylation.GFP.promoter.percent,2),
        p73.bs.with.methylation.TAa.if.confirmed.TAa.promoter.percent = round(p73.bs.with.methylation.TAa.if.confirmed.TAa.promoter.percent,2),
        p73.bs.with.methylation.DNb.if.confirmed.DNb.promoter.percent = round(p73.bs.with.methylation.DNb.if.confirmed.DNb.promoter.percent,2),
        p73.bs.with.methylation.GFP.if.confirmed.GFP.promoter.percent = round(p73.bs.with.methylation.GFP.if.confirmed.GFP.promoter.percent,2),
        p73.bs.with.confirmation.TAa.only.percent = round(p73.bs.with.confirmation.TAa.only.percent,2),
        p73.bs.with.confirmation.DNb.only.percent = round(p73.bs.with.confirmation.DNb.only.percent,2),
        p73.bs.with.confirmation.GFP.only.percent = round(p73.bs.with.confirmation.GFP.only.percent,2),
        p73.bs.with.confirmation.TAa.only.promoter.percent = round(p73.bs.with.confirmation.TAa.only.promoter.percent,2),
        p73.bs.with.confirmation.DNb.only.promoter.percent = round(p73.bs.with.confirmation.DNb.only.promoter.percent,2),
        p73.bs.with.confirmation.GFP.only.promoter.percent = round(p73.bs.with.confirmation.GFP.only.promoter.percent,2),
        p73.bs.with.methylation.TAa.only.confirmed.TAa.promoter.percent = round(p73.bs.with.methylation.TAa.only.confirmed.TAa.promoter.percent,2),
        p73.bs.with.methylation.DNb.only.confirmed.DNb.promoter.percent = round(p73.bs.with.methylation.DNb.only.confirmed.DNb.promoter.percent,2),
        p73.bs.with.methylation.GFP.only.confirmed.GFP.promoter.percent = round(p73.bs.with.methylation.GFP.only.confirmed.GFP.promoter.percent,2),

        Quantile.TAa.0   = quantile.TAa[1],
        Quantile.TAa.50  = quantile.TAa[2],
        Quantile.TAa.75  = quantile.TAa[3],
        Quantile.TAa.90  = quantile.TAa[4],
        Quantile.TAa.95  = quantile.TAa[5],
        Quantile.TAa.99  = quantile.TAa[6],
        Quantile.TAa.100 = quantile.TAa[7],
        Quantile.TAa.95.promoter  = quantile.TAa.promoter[5],
        Quantile.TAa.99.promoter  = quantile.TAa.promoter[6],
        Quantile.TAa.100.promoter = quantile.TAa.promoter[7],
        Quantile.DNb.0   = quantile.DNb[1],
        Quantile.DNb.50  = quantile.DNb[2],
        Quantile.DNb.75  = quantile.DNb[3],
        Quantile.DNb.90  = quantile.DNb[4],
        Quantile.DNb.95  = quantile.DNb[5],
        Quantile.DNb.99  = quantile.DNb[6],
        Quantile.DNb.100 = quantile.DNb[7],
        Quantile.DNb.90.promoter  = quantile.DNb.promoter[4],
        Quantile.DNb.95.promoter  = quantile.DNb.promoter[5],
        Quantile.DNb.99.promoter  = quantile.DNb.promoter[6],
        Quantile.DNb.100.promoter = quantile.DNb.promoter[7]
    )

    return(summary.table)
}

distribution.result.table.filename <- "distribution.result.table.RData"

if (!meta.from.scratch && file.exists(distribution.result.table.filename)) {
    cat("I: Loading existing distribution result table from '", distribution.result.table.filename, "'...\n", sep="")
    load(distribution.result.table.filename)
    if (!exists("distribution.result.table")) {
        stop("E: 'distribution.result.table' not found in '", distribution.result.table.filename, "'.")
    }
    cat("I: Found ", ncol(distribution.result.table), " columns in distribution result table.\n", sep="")
} else {
    cat("I: No existing distribution result table found, creating a new one...\n")

    distribution.result.table <- sapply(names(m.contexts), function(i) {
        cat("I: processing chromosome ", i, "...\n", sep = "")
        if (is.null(m.contexts[[i]])) {
            cat("E: No data for chromosome ", i, " - skipping\n", sep = "")
            return(NULL)
        }
        # Analyze the context for this chromosome and return properties
        analyze.context.for.chromosome(m.contexts[[i]], i)
    })

    distribution.result.table.summary <- (t(apply(distribution.result.table, 1, function(x) {
        round(c(summary(as.numeric(x)),sd=sd(as.numeric(x))),1)
    })))

    save(distribution.result.table, distribution.result.table.summary, file=distribution.result.table.filename)
}

#
# Prepare Volcano plots
#
source("analyze_matrix_function_volcano.R")


# All chromosomal results taken together, what fraction of TFBS are confirmed by CUT&RUN data in dependency of the score of the TFBS binding site?
cat("I: Analyzing all chromosomes for TFBS confirmation by CUT&RUN data...\n")

require(ggplot2)

# FIXME: Find better name and possibly split into multiple files
source("analyze_matrix_function_lists_genomewide.R")

# Figure 1C can be recreated with data collected at this point
# source("analyze_matrix_figure1C.R")

#
# Determine genes 1.5fold change upon DN-overexpression but not for TA-overexpression

#              TA                                  DN
DN.enriched <- combined.expression.data[,23]<=1.25 & combined.expression.data[,28]>=1.5
#              TA                                  DN
TA.enriched <- combined.expression.data[,23]>=1.5 & combined.expression.data[,28]<=1.25

combined.expression.data.valid <- !is.na(combined.expression.data$From) & !is.na(combined.expression.data$To) & !is.na(combined.expression.data$Gene)

print(paste0("I: Found ",sum(DN.enriched)," genes enriched in DN-overexpression but not in TA-overexpression."))
print(paste0("I: Found ",sum(TA.enriched)," genes enriched in TA-overexpression but not in DA-overexpression."))

DN.enriched.valid <- DN.enriched & combined.expression.data.valid
TA.enriched.valid <- TA.enriched & combined.expression.data.valid

print(paste0("I: Found ",sum(DN.enriched.valid)," genes enriched in DN-overexpression but not in TA-overexpression (valid genes)."))
print(paste0("I: Found ",sum(TA.enriched.valid)," genes enriched in TA-overexpression but not in DN-overexpression (valid genes)."))

source("analyze_matrix_function_retrieve_context.R")

tp73_noConfirmationWhatsoever <- retrieve_context_data_by_chromosome(NULL,confirmation="none",TA.or.DN="any")
tp73_inPromoter <- retrieve_context_data_by_chromosome(NULL,confirmation="promoter",TA.or.DN="any")


tp73_tp73ConfirmTA <- retrieve_context_data_by_chromosome(NULL,confirmation="tp73",TA.or.DN="TA"); gc()
tp73_tp73ConfirmDN <- retrieve_context_data_by_chromosome(NULL,confirmation="tp73",TA.or.DN="DN"); gc()

ta_vs_dn_tp73Confirm <- cbind(mean.TA=tp73_tp73ConfirmTA$mean_total,
                              mean.DN=tp73_tp73ConfirmDN$mean_total,
                              ratio=tp73_tp73ConfirmTA$mean_total / tp73_tp73ConfirmDN$mean_total)
ta_vs_dn_tp73Confirm.sorted <- ta_vs_dn_tp73Confirm[order(ta_vs_dn_tp73Confirm[, "ratio"]), ]
rownames(ta_vs_dn_tp73Confirm.sorted) <- prettyIdentifierJaspar(rownames(ta_vs_dn_tp73Confirm.sorted))

tp73_tp73ConfirmTA_inPromoter <- retrieve_context_data_by_chromosome(NULL,confirmation=c("tp73","promoter"),TA.or.DN="TA")
tp73_tp73ConfirmDN_inPromoter <- retrieve_context_data_by_chromosome(NULL,confirmation=c("tp73","promoter"),TA.or.DN="DN")
ta_vs_dn_tp73Confirm_inPromoter <- cbind(mean.TA=tp73_tp73ConfirmTA_inPromoter$mean_total_binary,
                                         mean.DN=tp73_tp73ConfirmDN_inPromoter$mean_total_binary,
                                         ratio=tp73_tp73ConfirmTA_inPromoter$mean_total_binary / tp73_tp73ConfirmDN_inPromoter$mean_total_binary)
ta_vs_dn_tp73Confirm_inPromoter.sorted <- ta_vs_dn_tp73Confirm_inPromoter[order(ta_vs_dn_tp73Confirm_inPromoter[, "ratio"]), ]
rownames(ta_vs_dn_tp73Confirm_inPromoter.sorted) <- prettyIdentifierJaspar(rownames(ta_vs_dn_tp73Confirm_inPromoter.sorted))


tp73_tp73ConfirmTA_inPromoter_posConfirmTA <- retrieve_context_data_by_chromosome(NULL,confirmation=c("tp73","pos","promoter"),TA.or.DN="TA")
tp73_tp73ConfirmDN_inPromoter_posConfirmDN <- retrieve_context_data_by_chromosome(NULL,confirmation=c("tp73","pos","promoter"),TA.or.DN="DN")
ta_vs_dn_tp73Confirm_inPromoter_posConfirm <- cbind(
                                        mean.TA=tp73_tp73ConfirmTA_inPromoter_posConfirmTA$mean_total,
                                        mean.DN=tp73_tp73ConfirmDN_inPromoter_posConfirmDN$mean_total,
                                        ratio=tp73_tp73ConfirmTA_inPromoter_posConfirmTA$mean_total / tp73_tp73ConfirmDN_inPromoter_posConfirmDN$mean_total)
ta_vs_dn_tp73Confirm_inPromoter_posConfirm.sorted <- ta_vs_dn_tp73Confirm_inPromoter_posConfirm[order(ta_vs_dn_tp73Confirm_inPromoter_posConfirm[, "ratio"]), ]
rownames(ta_vs_dn_tp73Confirm_inPromoter_posConfirm.sorted) <- prettyIdentifierJaspar(rownames(ta_vs_dn_tp73Confirm_inPromoter_posConfirm.sorted))


dn_results_tp73ConfirmAny <- retrieve_context_data_by_chromosome(DN.enriched.valid,confirmation="tp73",TA.or.DN="any")
ta_results_tp73ConfirmAny <- retrieve_context_data_by_chromosome(TA.enriched.valid,confirmation="tp73",TA.or.DN="any")
#all(dn_results_tp73ConfirmAny$mean_by_chromosome == dn_results$mean_by_chromosome)

#dn_results <- retrieve_context_data_by_chromosome(DN.enriched.valid)
#ta_results <- retrieve_context_data_by_chromosome(TA.enriched.valid)

ta_effects_tp73ConfirmTA_posConfirmTA <- retrieve_context_data_by_chromosome(TA.enriched.valid,confirmation=c("tp73","pos"),TA.or.DN="TA")
dn_effects_tp73ConfirmDN_posConfirmDN <- retrieve_context_data_by_chromosome(DN.enriched.valid,confirmation=c("tp73","pos"),TA.or.DN="DN")
ta_vs_effects_tp73Confirm_posConfirm <- cbind(
    mean.TA=ta_effects_tp73ConfirmTA_posConfirmTA$mean_total,
    mean.DN=dn_effects_tp73ConfirmDN_posConfirmDN$mean_total,
    ratio=ta_effects_tp73ConfirmTA_posConfirmTA$mean_total / dn_effects_tp73ConfirmDN_posConfirmDN$mean_total
)
ta_vs_effects_tp73Confirm_posConfirm.sorted <- ta_vs_effects_tp73Confirm_posConfirm[order(ta_vs_effects_tp73Confirm_posConfirm[, "ratio"]), ]
rownames(ta_vs_effects_tp73Confirm_posConfirm.sorted) <- prettyIdentifierJaspar(rownames(ta_vs_effects_tp73Confirm_posConfirm.sorted))

cat("I: Retrieved context data, calculated colSums, number of matches, and mean for DN.enriched.valid and TA.enriched.valid rows successfully.\n")

combined.expression.in.promoters.of.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <- combined.expression.data[,"Gene Symbol"] %in% promoterBedTables$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.promoter.bed$Gene
genesetEMT_tp73ConfirmTA <- retrieve_context_data_by_chromosome(combined.expression.in.promoters.of.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,confirmation=c("tp73"),TA.or.DN="TA")
genesetEMT_tp73ConfirmDN <- retrieve_context_data_by_chromosome(combined.expression.in.promoters.of.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,confirmation=c("tp73"),TA.or.DN="DN")
ta_vs_dn_genesetEMT_tp73Confirm <- cbind(
    mean.TA=genesetEMT_tp73ConfirmTA$mean_total,
    mean.DN=genesetEMT_tp73ConfirmDN$mean_total,
    ratio=genesetEMT_tp73ConfirmTA$mean_total / genesetEMT_tp73ConfirmDN$mean_total
)
ta_vs_dn_genesetEMT_tp73Confirm.sorted <- ta_vs_dn_genesetEMT_tp73Confirm[order(ta_vs_dn_genesetEMT_tp73Confirm[, "ratio"]), ]
rownames(ta_vs_dn_genesetEMT_tp73Confirm.sorted) <- prettyIdentifierJaspar(rownames(ta_vs_dn_genesetEMT_tp73Confirm.sorted))


combined.expression.in.promoters.of.EMT_TA_induced <- combined.expression.data[,"Gene Symbol"] %in% promoterBedTables$Nico_Analysis_TA_20250508.promoter.bed$Gene
genesetEMT.TA_tp73ConfirmTA_posConfirmTA <- retrieve_context_data_by_chromosome(combined.expression.in.promoters.of.EMT_TA_induced,confirmation=c("tp73","pos"),TA.or.DN="TA")
genesetEMT.TA_tp73ConfirmDN_posConfirmDN <- retrieve_context_data_by_chromosome(combined.expression.in.promoters.of.EMT_TA_induced,confirmation=c("tp73","pos"),TA.or.DN="DN")
print(rbind(genesetEMT.TA_tp73ConfirmTA_posConfirmTA$num_matches_by_chromosome, genesetEMT.TA_tp73ConfirmDN_posConfirmDN$num_matches_by_chromosome))
ta_vs_dn_genesetEMT.TA_tp73Confirm_posConfirm <- cbind(
    mean.TA=genesetEMT.TA_tp73ConfirmTA_posConfirmTA$mean_total,
    mean.DN=genesetEMT.TA_tp73ConfirmDN_posConfirmDN$mean_total,
    ratio=genesetEMT.TA_tp73ConfirmTA_posConfirmTA$mean_total / genesetEMT.TA_tp73ConfirmDN_posConfirmDN$mean_total
)
ta_vs_dn_genesetEMT.TA_tp73Confirm_posConfirm.sorted <- ta_vs_dn_genesetEMT.TA_tp73Confirm_posConfirm[order(ta_vs_dn_genesetEMT.TA_tp73Confirm_posConfirm[, "ratio"]), ]
rownames(ta_vs_dn_genesetEMT.TA_tp73Confirm_posConfirm.sorted) <- prettyIdentifierJaspar(rownames(ta_vs_dn_genesetEMT.TA_tp73Confirm_posConfirm.sorted))


# Not used - external EMT gene list
genesetEMT_tp73ConfirmTA_posConfirmTA <- retrieve_context_data_by_chromosome(combined.expression.in.promoters.of.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,confirmation=c("tp73","pos"),TA.or.DN="TA")
genesetEMT_tp73ConfirmDN_posConfirmDN <- retrieve_context_data_by_chromosome(combined.expression.in.promoters.of.HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,confirmation=c("tp73","pos"),TA.or.DN="DN")
ta_vs_dn_genesetEMT_tp73Confirm_posConfirm <- cbind(
    mean.TA=genesetEMT_tp73ConfirmTA_posConfirmTA$mean_total,
    mean.DN=genesetEMT_tp73ConfirmDN_posConfirmDN$mean_total,
    ratio=genesetEMT_tp73ConfirmTA_posConfirmTA$mean_total / genesetEMT_tp73ConfirmDN_posConfirmDN$mean_total
)
ta_vs_dn_genesetEMT_tp73Confirm_posConfirm.sorted <- ta_vs_dn_genesetEMT_tp73Confirm_posConfirm[order(ta_vs_dn_genesetEMT_tp73Confirm_posConfirm[, "ratio"]), ]
rownames(ta_vs_dn_genesetEMT_tp73Confirm_posConfirm.sorted) <- prettyIdentifierJaspar(rownames(ta_vs_dn_genesetEMT_tp73Confirm_posConfirm.sorted))


# TA-/DN-specific gene selection for EMT

#gene.selection.EMT.TAonly <- c("EDIL3","MMP2","LRRC15","LAMA3","BGN","FBLN2","ACTA2","TGFBR3","TNFRSF12A","VCAN","MYLK","SERPINH1","TPM1","COL8A2","CRLF1","LGALS1","COL5A2","SLIT2","MATN2")
#gene.selection.EMT.DNonly <- c("LAMA1","CDH11","SPARC","EMP3","DAB2","FN1","THBS1","DKK1","COLGALT1","PCOLCE","VIM","SNAI2","MCM7","ITGAV","PMP22","PLOD1")
#gene.selection.EMT.TAandDN <- c("CD44","PMEPA1","SFRP1","GADD45B","CADM1","SLC6A8","ID2","SDC1")
gene.selection.EMT.TA <- promoterBedTables$Nico_Analysis_TA_20250618.promoter.bed$Gene
gene.selection.EMT.DN <- promoterBedTables$Nico_Analysis_DN_20250618.promoter.bed$Gene
gene.selection.EMT.TAandDN <- intersect(gene.selection.EMT.TA, gene.selection.EMT.DN)
gene.selection.EMT.TAorDN <- unique(c(gene.selection.EMT.TA, gene.selection.EMT.DN, gene.selection.EMT.TAandDN))
gene.selection.EMT.TAonly <- gene.selection.EMT.TA[!gene.selection.EMT.TA %in% gene.selection.EMT.DN]
gene.selection.EMT.DNonly <- gene.selection.EMT.DN[!gene.selection.EMT.DN %in% gene.selection.EMT.TA]

genesetEMT.TA_tp73ConfirmTA_posConfirmTA <- retrieve_context_data_by_chromosome(combined.expression.data[,"Gene Symbol"] %in% gene.selection.EMT.TA, confirmation=c("tp73","pos"),TA.or.DN="TA")
genesetEMT.DN_tp73ConfirmDN_posConfirmDN <- retrieve_context_data_by_chromosome(combined.expression.data[,"Gene Symbol"] %in%gene.selection.EMT.DN, confirmation=c("tp73","pos"),TA.or.DN="DN")
ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm <- cbind(
    mean.TA=genesetEMT.TA_tp73ConfirmTA_posConfirmTA$mean_total,
    mean.DN=genesetEMT.DN_tp73ConfirmDN_posConfirmDN$mean_total,
    ratio=genesetEMT.TA_tp73ConfirmTA_posConfirmTA$mean_total / genesetEMT.DN_tp73ConfirmDN_posConfirmDN$mean_total
)
ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm.sorted <- ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm[order(ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm[, "ratio"]), ]
rownames(ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm.sorted) <- prettyIdentifierJaspar(rownames(ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm.sorted))

# Prepare data for vertical plot: top 5, bottom 5, and genes of interest

if (FALSE) {
    # To help with debugging, we can set these variables to NULL
    sorted.matrix.to.display=NULL
    num.from.extrema <- 10
    tf.patterns.of.interest <- NULL
}


source("analyze_matrix_plots_slots.R")


# Retrieve all transcripts associated with a given genomic region
transcript.associated.with <- function(chr, start, end) {
   if (!"PromoterOfWhichGene" %in% colnames(combined.expression.data)) {
     stop("E: Column 'PromoterOfWhichGene' not found in combined.expression.data.")
   }
   chr.vec <- as.character(combined.expression.data$Chr)
   chr.query <- as.character(chr)
   stopifnot(length(chr.query)==1)
   lines.of.interest <- rep(TRUE,length(chr.vec))
   lines.of.interest <- lines.of.interest & lines.of.interest & chr.vec == chr.query & !is.na(chr.vec)
   lines.of.interest <- lines.of.interest & as.numeric(combined.expression.data$From)-500 <= as.numeric(start)
   lines.of.interest <- lines.of.interest & as.numeric(combined.expression.data$To)+500 >= as.numeric(end)
    lines.of.interest.which <- which(lines.of.interest)
   if (length(lines.of.interest) > 0) {
     cat("I: Found ", length(lines.of.interest.which), " transcripts\n", sep = "")
     if (length(lines.of.interest.which)<10) print(lines.of.interest.which)
     return(combined.expression.data[lines.of.interest.which, "PromoterOfWhichGene"])
   } else {
     cat("I: Found no transcripts\n", sep = "")
     return(NA)
   }
}
#   combined.expression.data[combined.expression.data$Chr==chr & combined.expression.data$From>=start & combined.expression.data$To<=end, ]


# Genes with more than 10 p73 binding sites in 500 bp promoter
for(i in chromosomes) {
    cat("Chr ",i,": ",sep="")
    i.max <- max(m.contexts[[i]][m.contexts[[i]]$InPromoter, ]$"TP73_MA0861.1_NumInWindow")
    cat("Max: ", i.max, "\n")
    if (i.max >= 10) {
        i.m <- m.contexts[[i]][m.contexts[[i]]$InPromoter & m.contexts[[i]]$TP73_MA0861.1_NumInWindow>=10,
            c("Chr","From","To","Score",
            "pos_saos2_DN","pos_saos2_GFP","pos_saos2_TA",
            "pos_skmel29_2_DN","pos_skmel29_2_GFP","pos_skmel29_2_TA",
            "tp73_saos2_DN","tp73_saos2_GFP","tp73_saos2_TA",
            "tp73_skmel29_2_DN","tp73_skmel29_2_GFP","tp73_skmel29_2_TA","TP73_MA0861.1_NumInWindow") ]
        genes <- c()
        for(i.n in 1:nrow(i.m)) {
            g <-transcript.associated.with(i.m[i.n,"Chr"],i.m[i.n,"From"],i.m[i.n,"To"])
            if (0 == length(g)) {
                print(i.m[i.n,])
            }
            genes <- c(genes, g)
        }
        print(genes)
        print(i.m)
    }
}


plot.distance.distribution.pre.post.intern <- function(d,name, scaling="logarithmic") {

    if (! scaling %in% c("linear","logarithmic")) {
        stop("E: scaling must be 'linear' or 'logarithmic'.")
    }
    stopifnot(!any(is.na(d$Distance)) & !any(is.na(d$Coverage)))

    m.rect.y <- min(100,max(d$Coverage, na.rm = TRUE))
    m.rect <- matrix(NA, nrow=100+1+100, ncol=m.rect.y+1)
    for(i in 1:nrow(d)) {
        x <- d[i,"Distance"]
        stopifnot(-100 <= x & x <= 100)
        y <- d[i,"Coverage"]
        stopifnot(!is.na(y) & y>=0)

        if (y>m.rect.y) {
            y <- m.rect.y
        }
        #x.bin <- round((x - min(df_ta$Distance, na.rm = TRUE)) / (max(df_ta$Distance, na.rm = TRUE)-min(df_ta$Distance, na.rm = TRUE)) * 100) + 100 + 1
        x.bin <- x + 100 + 1
        y.bin <- y + 1
        # Ensure indices are within bounds
        if (x.bin >= 1 && x.bin <= nrow(m.rect) && y.bin >= 1 && y.bin <= ncol(m.rect) ) {
            if (is.na(m.rect[x.bin, y.bin])) {
                m.rect[x.bin, y.bin] <- 1
            } else {
                m.rect[x.bin, y.bin] <- m.rect[x.bin, y.bin] + 1
            }
        } else {
            cat("W: Skipping out-of-bounds index: x.bin=", x.bin, ", y=", y, "\n")
        }
    }

    for(j in (ncol(m.rect)-1):1) {
        m.rect[,j] <- apply(cbind(m.rect[,j],m.rect[,j+1]), 1, sum, na.rm=TRUE)
    }

    max.rect.rows <- apply(m.rect, 2, max)
    sum.rect.rows <- apply(m.rect, 2, sum)

    if (scaling=="logarithmic") {
        m.rect <- log2(m.rect + 1)
        m.rect <- sweep(m.rect, 2, log2(max.rect.rows+1), FUN="/")

    } else if (scaling=="linear") {

        #for(j in ncol(m.rect):1) {
        #    m.rect[,j] <- m.rect[,j] / max.rect.rows[j]
        #}
        sweep(m.rect, 2, max.rect.rows, FUN="/")
    }

    # Plot m.rect as an image with values as intensities
    # Add margin to the right for axis labels
    old.mar <- par("mar")
    par(mar = c(5, 4, 4, 8) + 0.1)  # increase right margin

    image(
        z = m.rect * 256,
        #col = c(0, colorRampPalette(c("blue", "red"))(255)),
        col = gray.colors(256, start = 1, end = 0),
        axes = FALSE,
        main = paste(name, " Coverage Heatmap for ", prettyIdentifierJaspar(tf), sep = ""),
        xlab = "Distance (bp)", ylab = paste(name,"Coverage (# reads)",sep=" ")
    )
    axis(1, at = seq(0, 1, length.out = 5), labels = seq(-100, 100, length.out = 5))
    axis(2, at = seq(0, 1, length.out = 5), labels = seq(0, m.rect.y, length.out = 5))
    axis(4, at = seq(0, 1, length.out = ncol(m.rect)), labels = sum.rect.rows, las = 2)
    mtext("# p73 binding sites", side = 4, line = 5, cex = 1)
    box()

    par(mar = old.mar)  # restore original margins
}


## Show distribution of distances prior and after filtering for confirmed p73 binding sites

plot.distance.distribution.pre.post <- function(tf,chr=chromosomes) {
    tf.colname <- paste0(tf,"_Shift")

    distances.pre <- distances.post.taOrDn <- distances.post.ta <- distances.post.dn <- distances.post.taOrDn.InPromoter <- numeric()
    coverage.ta <- coverage.dn <- numeric()
    inPromoter <- logical()

    for(i in chr) {
        if (!tf.colname %in% colnames(m.contexts[[i]])) {
            cat("E: TF '",tf,"' not found in chromosome '",i,"' - skipping\n",sep="")
            break;
        }
        if (is.null(m.contexts[[i]])) {
            stop("E: Chromosome '",i,"' provides not content.\n",sep="")
        }

        coverage.ta <- c(coverage.ta, m.contexts[[i]][["tp73_skmel29_2_TA"]])
        coverage.dn <- c(coverage.dn, m.contexts[[i]][["tp73_skmel29_2_DN"]])
        distances.pre <- c(distances.pre, m.contexts[[i]][[tf.colname]])
        inPromoter <- c(inPromoter, m.contexts[[i]]$InPromoter)


        distances.post.taOrDn <- c(distances.post.taOrDn,
                m.contexts[[i]][m.contexts[[i]]$"tp73_skmel29_2_TA">0 | m.contexts[[i]]$"tp73_skmel29_2_DN">0,][[tf.colname]])
        distances.post.taOrDn.InPromoter <- c(distances.post.taOrDn.InPromoter,
                                       m.contexts[[i]][
                                           (m.contexts[[i]]$"tp73_skmel29_2_TA">0 | m.contexts[[i]]$"tp73_skmel29_2_DN">0) & m.contexts[[i]]$InPromoter,][[tf.colname]])

    }
    stopifnot(length(coverage.ta) == length(coverage.dn))
    stopifnot(length(coverage.ta) == length(distances.pre))
    stopifnot(length(coverage.ta) == length(inPromoter))

    # Now filter for non-NA values to constrain on TF's presence
    distances.pre.nonNA <- distances.pre[!is.na(distances.pre)]
    coverage.ta.nonNA <- coverage.ta[!is.na(distances.pre)]
    coverage.dn.nonNA <- coverage.dn[!is.na(distances.pre)]
    inPromoter.nonNA <- inPromoter[!is.na(distances.pre)]

    distances.post.taOrDn2.nonNA <- distances.pre.nonNA[coverage.ta.nonNA>0 | coverage.dn.nonNA>0]
    distances.post.taOrDn.InPromoter2.nonNA <- distances.pre.nonNA[(coverage.ta.nonNA>0 | coverage.dn.nonNA>0) & inPromoter.nonNA]
    distances.post.ta <- distances.pre.nonNA[coverage.ta.nonNA>0]
    distances.post.dn <- distances.pre.nonNA[coverage.dn.nonNA>0]

    stopifnot(all(distances.post.taOrDn2.nonNA==distances.post.taOrDn[!is.na(distances.post.taOrDn)]))

    # Controlling reimplementation
    distances.post.taOrDn.InPromoter.nonNA <- distances.post.taOrDn.InPromoter[!is.na(distances.post.taOrDn.InPromoter)]
    stopifnot(all(distances.post.taOrDn.InPromoter2.nonNA==distances.post.taOrDn.InPromoter2.nonNA))


    # Plot side-by-side distribution of distances (pre vs post) for non-NA values

    if (length(distances.pre.nonNA) == 0 || length(distances.post.taOrDn2.nonNA) == 0) {
        cat("E: No valid distances found for TF '",tf,"'\n",sep="")
        return(NULL)
    }

    df <- data.frame(
        Distance = c(distances.pre.nonNA, distances.post.taOrDn2.nonNA, distances.post.taOrDn.InPromoter2.nonNA),
        Group = factor(rep(c(paste0("All TFBS (",length(distances.pre.nonNA),")"),
                            paste0("Confirmed TFBS (",length(distances.post.taOrDn2.nonNA),")"),
                            paste0("In Promoter (",length(distances.post.taOrDn.InPromoter2.nonNA),")")),
                           c(length(distances.pre.nonNA), length(distances.post.taOrDn2.nonNA), length(distances.post.taOrDn.InPromoter2.nonNA))))
    )

    p <- ggplot(df, aes(x=Distance, fill=Group)) +
        geom_density(alpha=0.5, position="identity") +
        theme_minimal() +
        labs(title=paste("Distribution of TF ",prettyIdentifierJaspar(tf)," Distances", sep=""), x="Distance", y="Density")
    print(p)

    # TA and DN coverage plots

    # Prepare data for TA coverage
    df_ta <- data.frame(
        Distance = distances.pre.nonNA,
        Coverage = coverage.ta.nonNA
    )

    stopifnot(!any(is.na(df_ta$Distance)) & !any(is.na(df_ta$Coverage)))

    plot.distance.distribution.pre.post.intern(df_ta,"TA")


    # Prepare data for DN coverage
    df_dn <- data.frame(
        Distance = distances.pre.nonNA,
        Coverage = coverage.dn.nonNA
    )

    stopifnot(!any(is.na(df_dn$Distance)) & !any(is.na(df_dn$Coverage)))

    plot.distance.distribution.pre.post.intern(df_dn,"DN")


    if (FALSE) {
        p_ta <- ggplot(df_ta, aes(x=Distance, y=TA_Coverage)) +
            geom_point(alpha=0.3, size=0.5) +
            theme_minimal() +
            labs(title=paste("TA CUT&RUN Coverage vs Distance for", prettyIdentifierJaspar(tf)),
                x="Distance", y="TA CUT&RUN Coverage")
        #print(p_ta)

        p_dn <- ggplot(df_dn, aes(x=Distance, y=DN_Coverage)) +
            geom_bin2d(bins = c(100, min(30,max(df_dn$DN_Coverage, na.rm = TRUE)))) +
            scale_fill_gradient(low="white", high="red") +
            theme_minimal() +
            labs(title=paste("DN CUT&RUN Coverage vs Distance for", prettyIdentifierJaspar(tf)),
                x="Distance", y="DN CUT&RUN Coverage")
        #print(p_dn)
    }
    #invisible(list(p, p_ta, p_dn))
    invisible(NULL)
}

tf.of.interest <- c(# JASPAR TFs from EMT heatmap
                    "TP73_MA0861.1","TP63_MA0525.2","TP53_MA0106.3",
                    "ATF7_MA0834.1",
                    "BACH2_MA1470.1",
                    "DBP_MA0639.1","DMRT3_MA0610.1",
                    "E2F7_MA0758.1","ETV2_MA0762.1","ELK4_MA0076.2","ESR1_MA0112.3",
                    "FOXF2_MA0030.1","FOXP2_MA0593.1","GFI1_MA0038.2","GATA5_MA0766.2",
                    "HES7_MA0822.1","HOXC9_MA0485.2","HOXC11_MA0651.2",
                    "JDP2_MA0656.1","JUND_MA0492.1",
                    "KLF13_MA0657.1","KLF15_MA1513.1",
                    "LMX1B_MA0703.2",
                    "MEF2B_MA0660.1","MEF2D_MA0773.1","MYBL2_MA0777.1",
                    "NFKB1_MA0105.4","NFKB2_MA0778.1",
                    "NFIC_MA0161.1","NFIC_MA0161.2","NFIC_MA1527.1",
                    "NFIX_MA1528.1","NRL_MA0842.2",
                    "ONECUT2_MA0756.2","OVOL2_MA1545.1",
                    "PATZ1_MA1961.1",
                    "PLAG1_MA0163.1",
                    "POU2F2_MA0507.2","POU1F1_MA0784.2","PPARG_MA0066.1","PROX1_MA0794.1",
                    "RELA_MA0107.1","REST_MA0138.2","RXRB_MA0855.1",
                    "SPDEF_MA0686.1",
                    "ZBTB32_MA1580.1",
                    # Prominent in global assessment
                    "CDX4_MA1473.1",
                    "HEY1_MA0823.1",
                    "HOXD11_MA0908.1",
                    "KLF1_MA0493.2",
                    "KLF7_MA1959.1",
                    "KLF12_MA0742.2",
                    "PLAGL2_MA1548.1",
                    "TEF_MA0843.1",
                    #"TFAP2A", # multiple, check ID
                    #"TFAP2B", # multiple, check ID
                    #"TFAP2C", # multiple, check ID
                    "TFAP2E_MA1569.1",
                    # Figure 3 - Enriched for TAp73a
                    "PAX1_MA0779.1",
                    #"CEBPB", # multiple, check ID
                    #"CEBPE", # multiple, check ID
                    #"GLI3", # multiple, check ID
                    #"SP3", # multiple, check ID
                    #"KLF10", # multiple, check ID
                    "MTF1_MA0863.1",
                    # Figure 3 - Enriched for DNp73a
                    #"NR1D2 ", # multiple, check ID
                    "MAZ_MA1522.1",
                    # Others
                    "IRF1_MA0050.1","IRF1_MA0050.2","IRF2_MA0051.1","IRF3_MA1418.1","IRF4_MA1419.1","IRF5_MA1420.1",
                    "IRF6_MA1509.1","IRF7_MA0772.1","IRF8_MA0652.1","IRF9_MA0653.1",
                    "OVOL1_MA1544.1",
                    "RORA_MA0071.1",
                    "RORA_MA0072.1",
                    "RUNX2_MA0511.2",
                    "SATB1_MA1963.1",
                    "STAT1_MA0137.1","STAT1_MA0137.2","STAT1_MA0137.3",
                    "SP1_MA0079.5",
                    "TBP_MA0108.2",
                    "YY1-2_MA1927.1",
                    "BCL6_MA0463.2",
                    "BCL6_MA0463.3",
                    "BCL6B_MA0731.1",
                    # Other species
                    "YAP1_MA0415.1"
)
f <- "distances.pdf"
pdf(f)
for(tf in tf.of.interest) {
    cat("I: Plotting distance distribution for TF '",tf,"'...\n",sep="")
    p.list <- plot.distance.distribution.pre.post(tf)
    if (is.null(p)) next();
    for(p in p.list) print(p)
}
dev.off()
cat("I: Image stored at '",Sys.info()["effective_user"],"@",Sys.info()["nodename"],":",getwd(),"/",f,"\n",sep="")



## manual ## source("analyze_matrix_plots_heatmap.R")

#
# OUTDATED (?) BELOW
#

# Load the findings if needed
# load("m.findings.RData")
# Check the structure of the findings
# str(m.findings)
# Check the names of the findings
# names(m.findings[[1]])
# Check the dimensions of the findings
# dim(m.findings[[1]][["ratio.skmel29_2.TAa.99.vs.ratio.skmel29_2.DNb.90"]])
# Check the first few rows of the findings
# head(m.findings[[1]][["ratio.skmel29_2.TAa.99.vs.ratio.skmel29_2.DNb.90"]])
# Check the first few rows and columns of the findings
# head(m.findings[[1]][["ratio.skmel29_2.TAa.99.vs.ratio.skmel29_2.DNb.90"]][1:20,1:10])

# Identify most promising transcription factors for 90% quantile
m.all.findings <- NULL
#lists <- c("ratio.skmel29_2.TAa.99.vs.ratio.skmel29_2.DNb.90","ratio.saos2.TAa.99.vs.ratio.saos2.DNb.90")
#lists <- c("ratio.skmel29_2.TAa.99.vs.ratio.skmel29_2.DNb.90","ratio.saos2.TAa.99.vs.ratio.saos2.DNb.90")
lists <- c("ratio.skmel29_2.TAa.quantile.90.vs.ratio.skmel29_2.DNb.quantile.90","ratio.saos2.TAa.99.vs.ratio.saos2.DNb.99")

#all(names(m.findings[[1]][["ratio.saos2.TAa.99.vs.ratio.saos2.DNb.90"]])==names(m.findings[[2]][["ratio.saos2.TAa.99.vs.ratio.saos2.DNb.90"]]))
for (l in lists) {
    if (! l %in% names(m.findings[[1]])) {
        cat("E: ",l," not found in m.findings[[1]] - skipping\n",sep="")
        next
    }
    cat("I: processing findings for ",l,"...\n",sep="")
    m.all.findings[[l]] <- NULL
    for(i in chromosomes) {
        m.all.findings[[l]] <- cbind(m.all.findings[[l]],m.findings[[i]][[l]])
    }

    if (any(is.na(m.all.findings[[l]]))) {
        cat("I: Found ",sum(is.na(m.all.findings[[l]]))," NA values in m.all.findings[[",l,"]], replacing with 100\n",sep="")
        m.all.findings[[l]][is.na(m.all.findings[[l]])] <- 100
    }

    if (any(is.na(m.all.findings[[l]]))) {
        stop(paste("Still finding NAs for ''",l,"''.",sep=""))
    }

    rownames(m.all.findings[[l]]) <- prettyIdentifierJaspar(rownames(m.all.findings[[l]]))
}


motifs.of.interest <- c("TP73_MA0861.1","TP63_MA0525.2","TP53_MA0106.3", "CTCF_MA0139.1",
            "E2F1_MA0024.3", "E2F2_MA0864.2", "E2F4_MA0470.2", "E2F6_MA0471.2",
            "FOS_MA1951.1","FOSL1--JUN_MA1129.1",
            "FOXO4_MA0848.1","Foxq1_MA0040.1", "KLF14_MA0740.2","Klf15_MA1890.1","REL_MA0101.1", "PLAG1_MA0163.1","POU2F1--SOX2_MA1962.1", "PPARG_MA0066.1",
            "RELA_MA0107.1","REST_MA0138.2", "RXRA--1VDR_MA0074.1", "SP1_MA0079.5", "TBP_MA0108.2","YY1-2_MA1927.1"
)

# Add additional motifs of interest - Nico's set of DN-associated TFs
motifs.of.interest <-  c(motifs.of.interest, "ARALYDRAFT_496250", "Ar_MA0007.3", "dve_MA0915.1", "Gsc_MA0190.1", "OLIG1_MA0826.1",
                                             "pnr_MA0536.1", "Ptx1_MA0201.1", "RORA_MA0071.1", "TCP13_MA1282.1", "TCP5_MA1067.1",
                                             "TDA9_MA0431.1")

# Add additional motifs of interest - Nico's set of TA-associated TFs
motifs.of.interest <-  c(motifs.of.interest, "ABF1_MA0570.2", "bHLH18_MA1361.1", "br_MA0010.1", "BRN2_MA1741.1", "Dmbx1_MA0883.1", "ETV2_MA0762.1",
                                             "ETV5_MA0765.3", "FEV_MA0156.3", "FOXD1_MA0031.1", "FOXO1--ELK3_MA1955.1", "GT-1_MA1020.1", "HEY1_MA0823.1",
                                             "HOXD9_MA0913.2", "MYB10_MA1762.1", "MYB17_MA1766.1", "MYB40_MA1770.1", "MYB41_MA1771.1", "MYB44_MA2027.1",
                                             "MYB51_MA1773.1", "MYB73_MA1394.2", "MYB80_MA1778.1", "MYB92_MA1780.1", "NAC010_MA2009.1", "NAC031_MA2045.1",
                                             "NAC075_MA2035.1", "OJ1058_F05.8", "RFX1_MA0365.1", "SPDEF_MA0686.1", "ZNF680_MA1729.1")

m.all.findings.summary <- NULL
for(l in lists) {
    cat("I: processing summary for ",l,"...\n",sep="")
    m.all.findings.summary[[l]] <- t(apply(m.all.findings[[l]],1,function(X) as.numeric(summary(X))))
    print(dim(m.all.findings.summary[[l]]))
    colnames(m.all.findings.summary[[l]]) <- c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
    rownames(m.all.findings.summary[[l]]) <- prettyIdentifierJaspar(rownames(m.all.findings.summary[[l]]))

    cat("\nConsistent findings for ",l,":\n",sep="")
    p <- m.all.findings.summary[[l]][order(m.all.findings.summary[[l]][,"Median"]),]
    cat("\n   top 20:\n")
    print(tail(p[p[,"1st Qu."]>1,],20))
    cat("\n   bottom 20:\n")
    print(head(p[p[,"3rd Qu."]<1,],20))
    cat("\n   motifs of interest:\n")
    print(p[prettyIdentifierJaspar(motifs.of.interest),])
    cat("\n   motifs of interest - details:\n")
    print(m.all.findings[[l]][prettyIdentifierJaspar(motifs.of.interest),])
}

# Plot distance of binding sites to CUT&RUN-confirmed p73 binding sites

l <- list("TAa"=sum.cutandrun.tp73.TAa>0, "DNb"=sum.cutandrun.tp73.DNb>0, "GFP"=sum.cutandrun.tp73.GFP>0,
          "TAa_without_DNb"=sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.90 &
                            sum.cutandrun.tp73.DNb<=sum.cutandrun.tp73.DNb.quantile.50,
          "DNb_without_TAa"=sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.90 &
                            sum.cutandrun.tp73.TAa<=sum.cutandrun.tp73.TAa.quantile.50)

for (l.name in names(l)) {
    for(tf in motifs.of.interest) {
        cat("I: plotting '",tf,"' (",l.name,") : ",sep="")
        plotShiftOfBinding(tf=tf,tfbs.selection=l[[l.name]],subset.name=l.name)
        cat("\n")
    }
}

