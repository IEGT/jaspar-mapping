#!/usr/bin/R

# Script to interpret TF binding sites for their association with CUT&RUN data.

options(width=180)

# Remove prior data on ratios or sums
rm(list=grep(ls(),pattern="^sum*",value=T))
rm(list=grep(ls(),pattern="^ratio.*",value=T))
rm(list=grep(ls(),pattern="^quantiles.*",value=T))
rm(list=grep(ls(),pattern="^mean.NumInWindow*",value=T))

# Import function to read and interpret the matrix for individual chromosomes
source("analyze_matrix_function_lists.R")  
source("analyze_matrix_function_distances.R")  
source("analyze_matrix_function_promoters.R")


chromosomes <- c(as.character(1:22) ) # ,"X","Y")

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


if (file.exists("combined.expression.data.RData")) {

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

    expressionData.dir <- "."
    # Load the expression data
    expressionData <- fread(file.path(expressionData.dir,"SkMel29_GFP_TAa_DNb_2x3x2_ohne_Filter_20.05.2025.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
    expressionData.symbols <- expressionData$"Gene Symbol"

    # Iteratate over genes in expression data to retrieve promoter locations and find max binding
    list.of.all.promoters <- promoterBedTables[["all.promoter.bed"]]
    max.binding.for.gene <- as.data.frame(matrix(NA, nrow=length(expressionData.symbols), ncol=19))
    colnames(max.binding.for.gene) <- c(colnames(m.contexts[[chr]])[1:18],"PromoterOfWhichGene")
    colnames(max.binding.for.gene)[4] <- "TF"
    colnames(max.binding.for.gene)[5] <- "TF.Score"
    colnames(max.binding.for.gene)[6] <- "TF.Strand"

    for(ed.rowno in 1:nrow(expressionData)) {
        gene <- expressionData.symbols[ed.rowno]
        cat(ed.nrowno,": ", gene, "\n",sep="")
        
        cat("I: processing gene ",gene,"...\n",sep="")
        # Check if promoter region of gene is known in promoters
        gene.in.promoter.rows <- (gene == list.of.all.promoters$"Gene")
        if (0 == sum(gene.in.promoter.rows)) {
            cat("E: Gene ",gene," not found in promoter data\n",sep="")
            next
        }
        # Get the promoter location
        gene.promoter.location <- list.of.all.promoters[gene.in.promoter.rows,c("Chr","From","To","Gene","Strand"),drop=F]
        print(gene.promoter.location)

        chr <- unique(gene.promoter.location$Chr)
        if (length(chr) != 1) {
            cat("E: Found ",length(chr)," chromosomes for gene ",gene," - skipping\n",sep="")
            next
        }

        if (is.null(m.contexts[[chr]])) {
            cat("W: No data for chromosome ",chr," - skipping\n",sep="")
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

    save(max.binding.for.gene, file="max.binding.for.gene.RData")


    # Check the dimensions of the max.binding.for.gene
    if(!all(max.binding.for.gene$"Gene" == expressionData.symbols)) {
        cat("E: Found ",sum(max.binding.for.gene$"Gene" != expressionData.symbols)," mismatches in gene names\n",sep="")
        print(max.binding.for.gene[max.binding.for.gene$"Gene" != expressionData.symbols,])
    } else {
        cat("I: Found all ",length(expressionData.symbols)," genes in max.binding.for.gene\n",sep="")
    }

    combined.expression.data <- cbind(max.binding.for.gene,expressionData)
    save(combined.expression.data, file="combined.expression.data.RData")
    require(openxlsx)
    write.xlsx(combined.expression.data, file="combined.expression.data.xlsx", rowNames=FALSE)
    cat("I: Saved combined expression data to 'combined.expression.data.RData' and 'combined.expression.data.xlsx'.\n")
}


# Classify all TFBS in m.contexts for being in promoter regions
cat("I: Classifying TFBS in promoter regions...\n")

if (file.exists("tfbs.in.promoter.list.RData")) {
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


# Assigning the promoter classification to m.contexts

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

save(distribution.result.table, distribution.result.table.summary, file="distribution.result.table.RData")


# All chromosomal results taken together, what fraction of TFBS are confirmed by CUT&RUN data in dependency of the score of the TFBS binding site?
cat("I: Analyzing all chromosomes for TFBS confirmation by CUT&RUN data...\n")

m.contexts.all.OneTo18 <- do.call(rbind, lapply(m.contexts, function(x) x[, 1:18, drop = FALSE]))
# Calculate the fraction of entries with value > 0 for "tp73_skmel29_2_DN" and "tp73_skmel29_2_TA" based on "Score"
bin.size <- 3
score_bins <- seq(0, max(m.contexts.all.OneTo18$Score, na.rm = TRUE), by = bin.size)
fractions <- data.frame(
    ScoreBin = score_bins[-length(score_bins)],
    Fraction_DN = sapply(1:(length(score_bins) - 1), function(i) {
        bin <- m.contexts.all.OneTo18$Score >= score_bins[i] & m.contexts.all.OneTo18$Score < score_bins[i + 1]
        sum(m.contexts.all.OneTo18$"tp73_skmel29_2_DN"[bin] > 0, na.rm = TRUE) / sum(bin, na.rm = TRUE)
    }),
    Fraction_TA = sapply(1:(length(score_bins) - 1), function(i) {
        bin <- m.contexts.all.OneTo18$Score >= score_bins[i] & m.contexts.all.OneTo18$Score < score_bins[i + 1]
        sum(m.contexts.all.OneTo18$"tp73_skmel29_2_TA"[bin] > 0, na.rm = TRUE) / sum(bin, na.rm = TRUE)
    }),
      # Add bars to indicate the number of TFBS within each bin
    TFBS_Count = sapply(1:(length(score_bins) - 1), function(i) {
        bin <- m.contexts.all.OneTo18$Score >= score_bins[i] & m.contexts.all.OneTo18$Score < score_bins[i + 1]
        sum(bin, na.rm = TRUE)
    })
)

require(ggplot2)

png("fractions_by_score.png", width=800, height=400)
# Plot the fractions
ggplot(fractions, aes(x = ScoreBin+bin.size/2)) +
    geom_bar(aes(y = TFBS_Count / max(TFBS_Count), fill = "p73 BS Count"), stat = "identity", alpha = 0.5) +
    geom_line(aes(y = Fraction_DN, color = "DN"), linewidth = 1) +
    geom_line(aes(y = Fraction_TA, color = "TA"), linewidth = 1) +
    scale_y_continuous(
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
