#!/usr/bin/R

# Script to interpret TF binding sites for their association with CUT&RUN data.

options(width=180)

# Remove prior data on ratios or sums
rm(list=grep(ls(),pattern="^sum*",value=T))
rm(list=grep(ls(),pattern="^ratio.*",value=T))
rm(list=grep(ls(),pattern="^quantiles.*",value=T))
rm(list=grep(ls(),pattern="^mean.NumInWindow*",value=T))

jaspar.human <- read.delim("jaspar_homo.tsv",row.names=1,col.names=FALSE)

# Import function to read and interpret the matrix for individual chromosomes
source("analyze_matrix_function_lists.R")
source("analyze_matrix_function_distances.R")
source("analyze_matrix_function_promoters.R")

chromosomes <- c(as.character(1:22),"X","Y")

m.findings.filename <- "m.findings.RData"
m.context.filename <- "m.contexts.RData"

meta.from.scratch <- TRUE

if (!meta.from.scratch && file.exists(m.context.filename) && file.exists(m.findings.filename)) {

    cat("Loading '",m.context.filename,"'\n")
    load(file=m.context.filename)
    cat("Loading '",m.findings.filename,"'\n")
    load(file=m.finding.filename)

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

expressionData.dir <- "."
expressionData.comparison.filename <- "SkMel29_GFP_TAa_DNb_2x3x2_ohne_Filter_20.05.2025.tsv"
num.skipped.because.of.gene.name.ambiguity <- 0
num.skipped.because.of.gene.name.unknown <- 0

if (!meta.from.scratch && file.exists("combined.expression.data.RData")) {

    cat("I: Loading existing combined expression data from 'combined.expression.data.RData'...\n")
    load("combined.expression.data.RData")
    if (!exists("combined.expression.data")) {
        stop("E: 'combined.expression.data' not found in 'combined.expression.data.RData'.")
    }
    cat("I: Found ",nrow(combined.expression.data)," rows in combined expression data.\n",sep="")
    print(head(combined.expression.data))

    # Check if max.binding.for.gene is available
    if (exists("max.binding.for.gene")) {
        cat("I: Found max.binding.for.gene with ",nrow(max.binding.for.gene)," rows.\n",sep="")
        print(head(max.binding.for.gene))
    } else {
        cat("W: max.binding.for.gene not found in 'combined.expression.data.RData'.\n")
    }

} else {

    require(openxlsx) # at the end to write the results to an Excel file

    # Load the expression data
    expressionData <- fread(file.path(expressionData.dir,expressionData.comparison.filename), sep="\t", header=TRUE, stringsAsFactors=FALSE)
    expressionData.symbols <- expressionData$"Gene Symbol"

    # Iteratate over genes in expression data to retrieve promoter locations and find max binding
    list.of.all.promoters <- promoterBedTables[["all.promoter.bed"]]
    max.binding.for.gene <- as.data.frame(matrix(NA, nrow=length(expressionData.symbols), ncol=19))
    colnames(max.binding.for.gene) <- c(colnames(m.contexts[[1]])[1:18],"PromoterOfWhichGene")
    colnames(max.binding.for.gene)[4] <- "TF"
    colnames(max.binding.for.gene)[5] <- "TF.Score"
    colnames(max.binding.for.gene)[6] <- "TF.Strand"

    for(ed.rowno in 1:nrow(expressionData)) {

        gene <- expressionData.symbols[ed.rowno]
        cat(ed.rowno,": ", gene, "\n",sep="")
        cat("I: processing gene ",gene,"...\n",sep="")
        # Check if promoter region of gene is known in promoters
        gene.in.promoter.rows <- (gene == list.of.all.promoters$"Gene")

        if (0 == sum(gene.in.promoter.rows)) {
            cat("E: Gene ",gene," not found in promoter data\n",sep="")
            max.binding.for.gene[ed.rowno,] <- c(NA,rep(NA,17),gene)
            num.skipped.because.of.gene.name.unknown <- num.skipped.because.of.gene.name.unknown + 1
            next
        }

        # Get the promoter location
        gene.promoter.location <- list.of.all.promoters[gene.in.promoter.rows,c("Chr","From","To","Gene","Strand"),drop=F]
        print(gene.promoter.location)

        chr <- unique(gene.promoter.location$Chr)
        if (length(chr) != 1) {
            cat("E: Found ",length(chr)," chromosomes for gene ",gene," - skipping\n",sep="")
            max.binding.for.gene[ed.rowno,] <- c(paste(chr,collapse="+",sep=""),rep(NA,17),gene)
            num.skipped.because.of.gene.name.ambiguity <- num.skipped.because.of.gene.name.ambiguity + 1
            next
        }

        if (is.null(m.contexts[[chr]])) {
            cat("W: No data for chromosome ",chr," - skipping\n",sep="")
            max.binding.for.gene[ed.rowno,] <- c(NA,rep(NA,17),gene)
            num.skipped.because.of.gene.name.unknown <- num.skipped.because.of.gene.name.unknown + 1
            next
        }

        cat("Chromosome: ",chr,"\n",sep="")
        overlap.index.of.m.chr <- checkBedOverlaps(gene.promoter.location[,1:3],m.contexts[[chr]][,1:3])

        if (0 < length(overlap.index.of.m.chr)) {

            cat("I: Found ",length(overlap.index.of.m.chr)," overlaps for gene ",gene," in chromosome ",chr,"\n",sep="")
            # Get the promoter location
            gene.promoter.location <- m.contexts[[chr]][overlap.index.of.m.chr,1:18]
            print(gene.promoter.location)
            val <- gene.promoter.location$"tp73_skmel29_2_DN" + gene.promoter.location$"tp73_skmel29_2_TA"
            val.max <- max(val)
            which.max.val <- which(val==val.max)
            select.line <- which.max.val
            if (length(select.line) > 1) {
                val.pos <- gene.promoter.location$"pos_skmel29_2_DN" + gene.promoter.location$"pos_skmel29_2_TA"
                val.pos[val != val.max] <- NA
                val.pos.max <- max(val.pos, na.rm=TRUE)
                select.line <- which(val.pos==val.pos.max)
                if (length(select.line) > 1) {
                    val.score <- gene.promoter.location$Score
                    val.score[val.pos != val.pos.max] <- NA
                    val.score[val != val.max] <- NA
                    val.score.max <- max(val.score, na.rm=TRUE)
                    select.line <- which.max(val.score)
                }
            }

            cat("I: Found max binding for gene ",gene," in chromosome ",chr," at line ",select.line," with value ",val.max,".\n",sep="")
            max.binding.for.gene[ed.rowno,] <- c(gene.promoter.location[select.line,],gene)
            #break;
        } else {
            max.binding.for.gene[ed.rowno,] <- c(chr,rep(NA,17),gene)
            cat("E: No overlaps found for gene ",gene," in chromosome ",chr,"\n",sep="")
        }
    }

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
    write.xlsx(combined.expression.data, file="combined.expression.data.xlsx", rowNames=FALSE)
    cat("I: Saved combined expression data to 'combined.expression.data.RData' and 'combined.expression.data.xlsx'.\n")
}


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
    m.contexts[[i]]$InPromoter <- tfbs.in.promoter.list[[i]][["all.promoter.bed"]]
}


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

svg("fractions_by_score.svg", width=800, height=400)
# Plot the fractions    
ggplot(fractions, aes(x = ScoreBin+bin.size/2)) +
    geom_bar(aes(y = TFBS_Count / max(TFBS_Count), fill = "p73 BS Count"), stat = "identity", alpha = 0.5) +
    geom_line(aes(y = Fraction_DN, color = "DN"), linewidth = 1) +
    geom_line(aes(y = Fraction_TA, color = "TA"), linewidth = 1) +
    scale_y_continuous(
        limits = c(0, 1),  # Ensure the y-axis covers 0 to 1 so 0.15 is visible
        name = "Fraction confirmed by CUT&RUN",
        sec.axis = sec_axis(~ . * max(fractions$TFBS_Count), name = "TFBS Count")
    ) +
    labs(
        title = "Fraction of Entries with Value > 0 by Score and TFBS Count",
        x = "Score",
        color = "Legend",
        fill = "Legend"
    ) +
    theme_minimal()
dev.off()


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

tp73_tp73ConfirmTA <- retrieve_context_data_by_chromosome(NULL,confirmation="tp73",TA.or.DN="TA")
tp73_tp73ConfirmDN <- retrieve_context_data_by_chromosome(NULL,confirmation="tp73",TA.or.DN="DN")

ta_vs_dn_tp73Confirm <- cbind(mean.TA=tp73_tp73ConfirmTA$mean_total,
                              mean.DN=tp73_tp73ConfirmDN$mean_total,
                              ratio=tp73_tp73ConfirmTA$mean_total / tp73_tp73ConfirmDN$mean_total)
ta_vs_dn_tp73Confirm.sorted <- ta_vs_dn_tp73Confirm[order(ta_vs_dn_tp73Confirm[, "ratio"]), ]
rownames(ta_vs_dn_tp73Confirm.sorted) <- prettyIdentifierJaspar(rownames(ta_vs_dn_tp73Confirm.sorted))

tp73_tp73ConfirmTA_inPromoter <- retrieve_context_data_by_chromosome(NULL,confirmation=c("tp73","promoter"),TA.or.DN="TA")
tp73_tp73ConfirmDN_inPromoter <- retrieve_context_data_by_chromosome(NULL,confirmation=c("tp73","promoter"),TA.or.DN="DN")
ta_vs_dn_tp73Confirm_inPromoter <- cbind(mean.TA=tp73_tp73ConfirmTA_inPromoter$mean_total,
                                         mean.DN=tp73_tp73ConfirmDN_inPromoter$mean_total,
                                         ratio=tp73_tp73ConfirmTA_inPromoter$mean_total / tp73_tp73ConfirmDN_inPromoter$mean_total)
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
gene.selection.EMT.TAonly <- promoterBedTables$Nico_Analysis_TA_20250618.promoter.bed$Gene
gene.selection.EMT.DNonly <- promoterBedTables$Nico_Analysis_DN_20250618.promoter.bed$Gene
gene.selection.EMT.TAandDN <- c()
gene.selection.EMT.TAorDN <- unique(c(gene.selection.EMT.TAonly, gene.selection.EMT.DNonly, gene.selection.EMT.TAandDN))
gene.selection.EMT.TA <- c(gene.selection.EMT.TAonly, gene.selection.EMT.TAandDN)
gene.selection.EMT.DN <- c(gene.selection.EMT.DNonly, gene.selection.EMT.TAandDN)

genesetEMT.TA_tp73ConfirmTA_posConfirmTA <- retrieve_context_data_by_chromosome(combined.expression.data[,"Gene Symbol"] %in% gene.selection.EMT.TA, confirmation=c("tp73","pos"),TA.or.DN="TA")
genesetEMT.DN_tp73ConfirmDN_posConfirmDN <- retrieve_context_data_by_chromosome(combined.expression.data[,"Gene Symbol"] %in%gene.selection.EMT.DN, confirmation=c("tp73","pos"),TA.or.DN="DN")
ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm <- cbind(
    mean.TA=genesetEMT.TA_tp73ConfirmTA_posConfirmTA$mean_total,
    mean.DN=genesetEMT.DN_tp73ConfirmDN_posConfirmDN$mean_total,
    ratio=genesetEMT.TA_tp73ConfirmTA_posConfirmTA$mean_total / genesetEMT.DN_tp73ConfirmDN_posConfirmDN$mean_total
)
ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm.sorted <- ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm[order(ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm[, "ratio"]), ]
rownames(ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm.sorted) <- prettyIdentifierJaspar(rownames(ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm.sorted))

# TA-/DN-specific gene expresssion change from Venn diagrams

require(readxl)
e <- read.xlsx("GeneLists/Venn Diagram (Gene mit hoher p73 Bindung) für Korrelationsgrafik.xlsx")
e.ta.up<-unlist(e[1,!is.na(e[1,]),drop=T])
e.dn.up<-unlist(e[2,!is.na(e[2,]),drop=T])

# Prepare data for vertical plot: top 5, bottom 5, and genes of interest
require(ggplot2)


if (FALSE) {
    # To help with debugging, we can set these variables to NULL
    sorted.matrix.to.display=NULL
    num.from.extrema <- 10
    tf.patterns.of.interest <- NULL
}


source("analyze_matrix_plots_slots.R")

source("analyze_matrix_plots_heatmap.R")

# Plot correlation matrix of all transcription factors (columns) against each other

require(corrplot)

gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant <- gene.selection_tp73ConfirmAny_context_interest_heatmap_input[!duplicated(gene.selection_tp73ConfirmAny_context_interest_heatmap_input[,"Gene"]),]
rownames(gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant) <- gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant$Gene


# Use the context data for all promoters with TP73 confirmation as the matrix
cor_matrix <- cor(gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant[,-1],
                     use="pairwise.complete.obs", method="pearson")

pdf("correlation_matrix_TFBS.pdf", width=12, height=10)
corrplot(cor_matrix, is.corr=TRUE, method="color", type="upper", order="hclust",
         tl.col="black", tl.cex=0.7, addCoef.col="black", number.cex=0.35, number.digits=1,insig="n",
         title="Correlation Matrix of TFBS (all promoters, TP73 confirmed)", mar=c(0,0,2,0))
dev.off()
cat("I: Correlation matrix plot saved to 'correlation_matrix_TFBS.pdf'.\n")

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

