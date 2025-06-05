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


# Iterate over chromosomes and retrieve context data for DN.enriched.valid and TA.enriched.valid rows
cat("I: Retrieving context data for DN.enriched.valid and TA.enriched.valid rows by iterating over chromosomes...\n")

retrieve_context_data_by_chromosome <- function(enriched_rows=NULL,confirmation=NULL,TA.or.DN=NULL) {
    context_data <- NULL
    cutandrun_data <- NULL
    coordinates <- NULL
    downstream_genes <- NULL
    context_matches <- c()

    if (length(setdiff(confirmation,c("none","tp73","pos","promoter")))>0) {
        stop("I: No method for confirmation specified: (",
                paste(setdiff(confirmation,c("none","tp73","pos","promoter")),collapse=",",sep=""),
                "), make it none, tp73, pos or c(tp73,pos).\n")
    } else {
        cat("I: Using method ", confirmation, ".\n", sep = "")
    }

    if (is.null(TA.or.DN) || ! TA.or.DN %in% c("TA","DN","any")) {
        stop("I: No enricment for TA or DN specified, make it TA or DN or any.\n")
    } else {
        cat("I: Using enrichment ", TA.or.DN, ".\n", sep = "")
    }

    # Create a matrix to store column sums for each chromosome
    mean_by_chromosome <- col_sums_by_chromosome <- matrix(
        NA,
        nrow = length(m.contexts),
        ncol = sum(cols.NumInWindow),
        dimnames = list(names(m.contexts), colnames(m.contexts[[1]][,..cols.NumInWindow]))
    )
    num_matches_by_chromosome <- setNames(vector("numeric", length(m.contexts)), names(m.contexts))

    for (chromosome in names(m.contexts)) {

        context_data_chromosome <- NULL
        cutandrun_data_chromosome <- NULL
        coordinates_chromosome <- NULL
        downstream_genes_chromosome <- NULL
        context_matches_chromosome <- c()
    

        if (is.null(m.contexts[[chromosome]])) {
            cat("E: No data for chromosome ", chromosome, " - skipping\n", sep = "")
            next
        }

        m.context.extra.check <- rep(TRUE, nrow(m.contexts[[chromosome]]))

        if ("promoter" %in% confirmation) {
            m.context.extra.check <- m.context.extra.check & m.contexts[[chromosome]][, "InPromoter"]
        }

        if ("tp73" %in% confirmation) {
            
            if (TA.or.DN == "TA") {
                m.context.extra.check <- m.context.extra.check & m.contexts[[chromosome]][, "tp73_skmel29_2_TA"] > 0
            } else if (TA.or.DN == "DN") {
                m.context.extra.check <- m.context.extra.check & m.contexts[[chromosome]][, "tp73_skmel29_2_DN"] > 0
            } else if (TA.or.DN == "any") {
                m.context.extra.check <- m.context.extra.check & (
                    m.contexts[[chromosome]][, "tp73_skmel29_2_TA"] > 0 |
                    m.contexts[[chromosome]][, "tp73_skmel29_2_DN"] > 0 |
                    m.contexts[[chromosome]][, "tp73_skmel29_2_GFP"] > 0
                )
            } else if (TA.or.DN == "none") {
                # doing nothing, as we want to check all rows
            } else {
                stop("E: Invalid TA.or.DN method specified.\n")
            }
        }

        if ("pos" %in% confirmation) {
            if (TA.or.DN == "TA") {
                m.context.extra.check <- m.context.extra.check & m.contexts[[chromosome]][, "pos_skmel29_2_TA"] > 0
            } else if (TA.or.DN == "DN") {
                m.context.extra.check <- m.context.extra.check & m.contexts[[chromosome]][, "pos_skmel29_2_DN"] > 0
            } else if (TA.or.DN == "any") { 
                m.context.extra.check <- m.context.extra.check & (
                    m.contexts[[chromosome]][, "pos_skmel29_2_DN"] > 0 |
                    m.contexts[[chromosome]][, "pos_skmel29_2_TA"] > 0 |
                    m.contexts[[chromosome]][, "pos_skmel29_2_GFP"] > 0
                    )
            } else if (TA.or.DN == "none") {
                # doing nothing, as we want to check all rows
            } else {
                stop("E: Invalid TA.or.DN method specified.\n")
            }                
        }

        if (is.null(enriched_rows)) {

            context_matches_chromosome <- which(m.context.extra.check)
            downstream_genes_chromosome <- rep(NA, length(context_matches_chromosome))
            cat("I: No enriched rows specified, using all rows for chromosome ", chromosome, ", found ",length(context_matches_chromosome)," hits.\n", sep = "")

        } else {

            rows_in_chromosome <- which(enriched_rows & combined.expression.data[, 1] == chromosome)
            if (length(rows_in_chromosome) == 0) {
                cat("E: No enriched rows found for chromosome ", chromosome, " - skipping\n", sep = "")
                next
            }
 
            for (row_idx in rows_in_chromosome) {

                start_pos <- combined.expression.data[row_idx, 2]
                end_pos <- combined.expression.data[row_idx, 3]
                gene <- combined.expression.data[row_idx, "Gene Symbol"]

                m.context.row.check <- m.context.extra.check & m.contexts[[chromosome]][, 2] == start_pos &
                                    m.contexts[[chromosome]][, 3] == end_pos

                match_idx <- which(m.context.row.check)
            
                if (length(match_idx) == 0) {
                    cat("E: No matching context found for row ", row_idx, " (chromosome: ", chromosome, ", start: ", start_pos, ", end: ", end_pos, ")\n", sep = "")
                    next
                }

                context_matches_chromosome <-  c(context_matches_chromosome, match_idx)
                downstream_genes_chromosome <- c(downstream_genes_chromosome, gene) 
            }
        }
        
        context_data_chromosome <- m.contexts[[chromosome]][context_matches_chromosome, ..cols.NumInWindow]
        cutandrun_data_chromosome <- m.contexts[[chromosome]][context_matches_chromosome, 7:18, drop = FALSE]
        coordinates_chromosome <- m.contexts[[chromosome]][context_matches_chromosome, 1:6, drop = FALSE]
        

        cat("I: Retrieved context for ", nrow(context_data_chromosome), " binding sites on chromosome: ", chromosome, "\n", sep = "")


        if (nrow(context_data_chromosome) == 0) {
            cat("E: No context data found for chromosome ", chromosome, " - skipping\n", sep = "")
            next
        }
        else {
            context_data <- rbind(context_data, context_data_chromosome)
            context_matches <- c(context_matches, context_matches_chromosome)
            cutandrun_data <- rbind(cutandrun_data, cutandrun_data_chromosome)
            coordinates <- rbind(coordinates, coordinates_chromosome)
            downstream_genes <- c(downstream_genes, downstream_genes_chromosome)
            col_sums_by_chromosome[chromosome,] <- colSums(context_data_chromosome, na.rm = TRUE)
            num_matches_by_chromosome[chromosome] <- nrow(context_data_chromosome)
            mean_by_chromosome[chromosome,] <- colMeans(context_data_chromosome, na.rm = TRUE)
            cat("I: Calculated colSums, number of matches, and mean for chromosome ", chromosome, "\n", sep = "")
        }
    }
    col_sums_total = colSums(col_sums_by_chromosome)
    return(list(
        context_data = context_data,
        context_matches = context_matches,
        downstream_genes = downstream_genes,
        cutandrun_data = cutandrun_data,
        coordinates = coordinates,
        col_sums_by_chromosome = col_sums_by_chromosome,
        col_sums_total = col_sums_total,
        num_matches_by_chromosome = num_matches_by_chromosome,
        num_matches_total = sum(num_matches_by_chromosome),
        mean_by_chromosome = mean_by_chromosome,
        mean_total = colMeans(context_data, na.rm = TRUE)
    ))
}

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

gene.selection.EMT.TAonly <- c("EDIL3","MMP2","LRRC15","LAMA3","BGN","FBLN2","ACTA2","TGFBR3","TNFRSF12A","VCAN","MYLK","SERPINH1","TPM1","COL8A2","CRLF1","LGALS1","COL5A2","SLIT2","MATN2")
gene.selection.EMT.DNonly <- c("LAMA1","CDH11","SPARC","EMP3","DAB2","FN1","THBS1","DKK1","COLGALT1","PCOLCE","VIM","SNAI2","MCM7","ITGAV","PMP22","PLOD1")
gene.selection.EMT.TAandDN <- c("CD44","PMEPA1","SFRP1","GADD45B","CADM1","SLC6A8","ID2","SDC1")
gene.selection.EMT.TAorDN <- c(gene.selection.EMT.TAonly, gene.selection.EMT.DNonly, gene.selection.EMT.TAandDN)
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


# Prepare data for vertical plot: top 5, bottom 5, and genes of interest
require(ggplot2)


if (FALSE) {
    # To help with debugging, we can set these variables to NULL
    sorted.matrix.to.display=NULL
    num.from.extrema <- 10
    tf.patterns.of.interest <- NULL
}

is.human.jaspar.id <- function(id) {
    if (0 == length(id)) {
        return(FALSE)
    }
    else if (length(id) == 1) {
        if ("" == id) {
            return(FALSE)
        }
        if (is.character(id) && length(id) == 0) {
            return(FALSE)
        }
        strsplit.id <- strsplit(id, split="[ _()]")[[1]]
        if (length(strsplit.id) == 1) {
            #cat("Found single ID: ", strsplit.id, "\n", sep="")
            return(strsplit.id %in% rownames(jaspar.human) || strsplit.id %in% jaspar.human)
        } else if (length(strsplit.id) == 2) {
            return(is.human.jaspar.id(strsplit.id[1]) || is.human.jaspar.id(strsplit.id[2]))
        } else if (length(strsplit.id) == 3) {
            return(is.human.jaspar.id(strsplit.id[1]) || is.human.jaspar.id(strsplit.id[2]) || is.human.jaspar.id(strsplit.id[3]))
        } else {
            stop("E: Unexpected ID format with ",length(id)," parts: ", id, "\n")
        }
    } else {
        return(sapply(id, is.human.jaspar.id))  
    }
}

#is.human.jaspar.id("(MA1466.1")

plot.series.of.ratio.matrices <- function(list.of.sorted.matrices.to.display=NULL,num.from.extrema=rep(10,length(list.of.sorted.matrices.to.display)),tf.patterns.of.interest=NULL,human.only=TRUE) {

    if (is.null(tf.patterns.of.interest)) {
        #tf.patterns.of.interest <- "^(jun|sp1|rest|yap1|yy1|e2f1|e2f7|tp53|tp63|tp73|plag1|rara|odd|irf2|bpc5|irf2|nr2f6|pdr3|pax1|rarg|esr1|odd|elt-3|nfkb1|rela|rora|erf6) "
        tf.patterns.of.interest <- paste0(
            #"^(jun|sp1|rest|yap1|yy1|e2f1|nfkb1|tp53|tp63|tp73|",
            #"BZIP43|IRF2|RXRA--VDR|pan|HAP1|Nr2F6|RORA|E2F3|E2F7|THI2|SPL15|TGIF1|MYB52|ZBED1|NAC083", #  paste(sapply(strsplit(rownames(head(ta_vs_effects_tp73Confirm_posConfirm.sorted,15)),split=" "),function(X) X[1]),collapse="|")
            #"NR3C1|TP63|RREB1|ATHB-9|GLIS1|PLAG1|E2FC|MYB15|DYT1|bZIP911|odd|RARA|BZIP42|BPC5|TCP14", # paste(sapply(strsplit(rownames(tail(ta_vs_effects_tp73Confirm_posConfirm.sorted[!is.na(ta_vs_effects_tp73Confirm_posConfirm.sorted[,"ratio"]),],15)),split=" "),function(X) X[1]),collapse="|")
            "RARA|GLIS1|PLAG1|RREB1|TP63|NR3C1|SRF|NR3C2|EWSR1-FLI1|THRB|TFAP2E|GLIS3|TFAP2C|RXRB|KLF14|",
            "IRF2|RXRA|RORA|E2F7|E2F3|TGIF1|ZBED1|IRF4|RFX3|MEF2D|STAT1|MEF2B|NFIC|IRF8|CREB3L1",
        ") ")
    }   

    if (is.null(list.of.sorted.matrices.to.display)) {
        cat("I: No sorted matrix provided, using default sorted matrices.\n")
        list.of.sorted.matrices.to.display <- list(
                effect_nonmethylated=ta_vs_dn_nonMethylated.sorted,
                effect_methylated=ta_vs_dn_w_methylation.sorted
        )
    } else if (!is.list(list.of.sorted.matrices.to.display)) {
        if (is.matrix(list.of.sorted.matrices.to.display) || is.data.frame(sorted.matrix.to.display)) {
            list.of.sorted.matrices.to.display <- list(unnamed=sorted.matrix.to.display)
        } else {
            cat("E: Invalid input for sorted.matrix.to.display, expected a matrix, data frame, or list.\n")
            stop("Exiting due to invalid input.")
        }
        cat("I: Using provided sorted matrix list.\n")
    } 

    if (0 == length(list.of.sorted.matrices.to.display)) {
        cat("E: No sorted matrices to display, exiting.\n")
        stop("Exiting due to no sorted matrices to display.")
    } else {
           cat("I: Found ", length(list.of.sorted.matrices.to.display), " sorted matrices to display:\n", sep="")
           cat("  ", paste(names(list.of.sorted.matrices.to.display), collapse=", "), "\n", sep="")
    }   

    pdf("vertical_plot.pdf", width=3*length(list.of.sorted.matrices.to.display)+1, height=8)

    for(i in 1:length(list.of.sorted.matrices.to.display)) {

        sorted.matrix.to.display.tag <- names(list.of.sorted.matrices.to.display)[i]
        add.to.existing.graph <- ( i > 1 )

        cat("I: Processing sorted matrix for tag: ", sorted.matrix.to.display.tag, "\n", sep="")

        print_legend <- (i == length(list.of.sorted.matrices.to.display))

        if (add.to.existing.graph) {
            cat("I: Adding to existing graph for tag: ", sorted.matrix.to.display.tag, "\n", sep="")
        } else {
            cat("I: Creating new graph for tag: ", sorted.matrix.to.display.tag, "\n", sep="")
        }

        sorted.matrix.to.display <- list.of.sorted.matrices.to.display[[sorted.matrix.to.display.tag]]
        print(dim(sorted.matrix.to.display))
        sorted.matrix.to.display <- sorted.matrix.to.display[which(is.human.jaspar.id(rownames(sorted.matrix.to.display))),,drop=FALSE]
        print(dim(sorted.matrix.to.display))

        sorted.matrix.to.display <- sorted.matrix.to.display[order(sorted.matrix.to.display[, "ratio"],decreasing=TRUE), , drop=FALSE]

        if (is.null(sorted.matrix.to.display)) {
            cat("E: No data for tag ", sorted.matrix.to.display.tag, " - skipping\n", sep="")
            next
        }

        sorted.matrix.to.display.valid <- rowSums(sorted.matrix.to.display[, c("mean.TA","mean.DN")])>0.0005

        # Get top and bottom
        top5 <- head(sorted.matrix.to.display[sorted.matrix.to.display.valid,], num.from.extrema[i])
        bottom5 <- tail(sorted.matrix.to.display[sorted.matrix.to.display.valid,], num.from.extrema[i])

        # Genes of interest (tf.of.interest) that are not in top5/bottom5
        tf_interest_names <- grep(rownames(ta_vs_dn_nonMethylated.sorted),pattern=tf.patterns.of.interest,ignore.case=T,value=T)
        print(tf_interest_names)

        tf_interest_set <- setdiff(tf_interest_names, c(rownames(top5), rownames(bottom5)))
        tf_interest <- sorted.matrix.to.display[rownames(sorted.matrix.to.display) %in% tf_interest_set, , drop=FALSE] # extra complicated to preserve order

        # Assign y positions: top5 (1:5), gap, tf_interest (by their rank), gap, bottom5 (last 5)
        n_top <- nrow(top5)
        n_bottom <- nrow(bottom5)
        gap_size <- 2

        # Get the ranks of tf_interest in the full sorted table
        tf_interest_ranks <- match(rownames(tf_interest), rownames(sorted.matrix.to.display))
        

        # Compose plotting data
        plot_data <- rbind(
            data.frame(Name=rownames(top5), mean.TA=top5[,1], mean.DN=top5[,2], ratio=top5[,3], group="Top", y=seq(1,n_top)),
            data.frame(Name=rownames(tf_interest), mean.TA=tf_interest[,1], mean.DN=tf_interest[,2], ratio=tf_interest[,3], group="Interest", y=tf_interest_ranks),
            data.frame(Name=rownames(bottom5), mean.TA=bottom5[,1], mean.DN=bottom5[,2], ratio=bottom5[,3], group="Bottom",
                        y=seq(sum(sorted.matrix.to.display.valid)-n_bottom+1,sum(sorted.matrix.to.display.valid)))
        )

        # Assign y for genes of interest: scale to fit in the gap
        gap_start <- 0
        gap_end <- -gap_size
        if(nrow(tf_interest)>0) {
            # Scale ranks to fit in the gap
            plot_data$y[plot_data$group=="Expression/Literature"] <- gap_start - (tf_interest_ranks - min(tf_interest_ranks)) / 
                (max(tf_interest_ranks)-min(tf_interest_ranks)+1) * gap_size
        }

        # Add a column for circle size (sum of means)
        plot_data$circle_size <- plot_data$mean.TA + plot_data$mean.DN

        # Assign y positions with reserved space: 20% for top5, 20% for bottom5, 60% for interest
        n_total <- sum(sorted.matrix.to.display.valid)
        n_top <- nrow(top5)
        n_bottom <- nrow(bottom5)
        n_interest <- nrow(tf_interest)

        # Define y-axis range
        y_min <- 0
        y_max <- 1

        # Allocate space
        top_space <- num.from.extrema[i]/50
        bottom_space <- num.from.extrema[i]/50
        interest_space <- 1-top_space-bottom_space

        # Calculate y positions for each group
        if (n_top > 1) {
            y_top <- seq(y_max, y_max - top_space, length.out = n_top + 1)[-1]
        } else {
            y_top <- y_max - top_space / 2
        }
        if (n_bottom > 1) {
            y_bottom <- seq(y_min + bottom_space, y_min , length.out = n_bottom+1)[-1]
        } else {
            y_bottom <- y_min + bottom_space / 2
        }
        if (n_interest > 0) {
            y_interest <- seq(y_max - top_space, y_min + bottom_space, length.out = n_interest + 2)[-c(1, n_interest + 2)]
        } else {
            y_interest <- numeric(0)
        }

        # Assign y to plot_data
        plot_data$y <- NA
        plot_data$y[plot_data$group == "Top"] <- y_top
        plot_data$y[plot_data$group == "Bottom"] <- y_bottom
        plot_data$y[plot_data$group == "Interest"] <- y_interest
        
        # Plot
        # Filter out rows with NA y values to avoid warnings and missing labels
        plot_data_valid <- plot_data[!is.na(plot_data$y), ]

        # Instead of using labs(title=...), add the tag as a text annotation at the correct x/y position
        # Place the label at the center top of the current column (x = (i-1)*2, y = max y + offset)
        p <- ggplot(plot_data_valid, aes(x=(i-1)*1.5, y=y)) +
            geom_point(aes(size=circle_size, fill=group), shape=21, color="black", alpha=0.8) +
            geom_text(aes(label=Name), hjust=0, nudge_x=0.12, size=2) +
            geom_text(aes(label=round(log2(ratio),2)), hjust=1, vjust=0.5, size=2, nudge_x=-0.12, fontface="bold") +
            scale_size_continuous(breaks = c(0.05,0.2,1,5), labels = c("0.05","0.2","1","5")) +
            scale_fill_manual(values=c("Top"="#1b9e77", "Interest"="#d95f02", "Bottom"="#7570b3")) +
            theme_minimal() +
            labs(size="Frequency", fill="Group") +
            theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
            panel.grid=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),
            legend.position=if (exists("print_legend") && print_legend) "right" else "none") +
            xlim(c(-0.75, (length(list.of.sorted.matrices.to.display)-1)*1.5+0.75)) +
            ylim(min(plot_data_valid$y, na.rm=TRUE), max(plot_data_valid$y, na.rm=TRUE)+0.2) +
            annotate("text", x=0.5+(i-1)*1.5, y=max(plot_data_valid$y, na.rm=TRUE)+0.1, na.rm=TRUE,
                 label=sorted.matrix.to.display.tag, fontface="bold", size=4, hjust=0.5, vjust=0)
 
        if (add.to.existing.graph) {
            print(p, newpage=FALSE)
        } else {
            print(p)
        }
        add.to.existing.graph <- TRUE

    }

    dev.off()
    cat("I: Vertical plot saved to 'vertical_plot.pdf'.\n")
}



list.of.matrices <- list(
    "p73 binding"=ta_vs_dn_tp73Confirm.sorted,
    "in promoter"= ta_vs_dn_tp73Confirm_inPromoter.sorted,
    "also methylated"=ta_vs_dn_tp73Confirm_inPromoter_posConfirm.sorted,
    # "effect nonmethylated"=ta_vs_dn_nonMethylated.sorted,
    #"effects expression"=ta_vs_dn_w_methylation.sorted
    "effects expression"=ta_vs_effects_tp73Confirm_posConfirm.sorted,
    #"geneset EMT"=ta_vs_dn_genesetEMT_tp73Confirm_posConfirm.sorted
    "TA/DN-induced EMT"=ta_vs_dn_genesetEMT.TA.DN_tp73Confirm_posConfirm.sorted
)

plot.series.of.ratio.matrices(
    list.of.sorted.matrices.to.display=list.of.matrices,
    num.from.extrema=c(10,10,10,15,15), # Adjusted to match the number of matrices
)


#
# Heatmap for gene selection
#


#gene.selection <- promoterBedTables$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.promoter.bed$Gene
gene.selection <- gene.selection.EMT.TAorDN
combined.expression.data.selected <- combined.expression.data[,"Gene Symbol"] %in% gene.selection
gene.selection_tp73ConfirmAny <- retrieve_context_data_by_chromosome(combined.expression.data.selected, confirmation=c("tp73"),TA.or.DN="any")
all_inPromoter_tp73ConfirmAny <- retrieve_context_data_by_chromosome(NULL, confirmation=c("tp73","promoter"),TA.or.DN="any")


cat("I: Distribution of number of matches by chromosome:\n")
print(gene.selection_tp73ConfirmAny$num_matches_by_chromosome)
cat("I: Total number of matches across all chromosomes: ",sum(gene.selection_tp73ConfirmAny$num_matches_by_chromosome),"\n",sep="")

gene.selection_tp73ConfirmAny.colSums <- colSums(gene.selection_tp73ConfirmAny$context_data,na.rm=T)
all_inPromoter_tp73ConfirmAny.colSums <- colSums(all_inPromoter_tp73ConfirmAny$context_data,na.rm=T)

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

gene.vs.all_inPromoter_tp73ConfirmAny.sorted <- gene.vs.all_inPromoter_tp73ConfirmAny[order(gene.vs.all_inPromoter_tp73ConfirmAny[, "ratio"],decreasing=T), ]
rownames(gene.vs.all_inPromoter_tp73ConfirmAny.sorted) <- prettyIdentifierJaspar(rownames(gene.vs.all_inPromoter_tp73ConfirmAny.sorted))
gene.vs.all_inPromoter_tp73ConfirmAny.rownames.interest <- c(rownames(gene.vs.all_inPromoter_tp73ConfirmAny.sorted)[is.human.jaspar.id(rownames(gene.vs.all_inPromoter_tp73ConfirmAny.sorted))][1:50],
    grep(x=prettyIdentifierJaspar(rownames(gene.vs.all_inPromoter_tp73ConfirmAny.sorted)),pattern=tf.patterns.of.interest,ignore.case=T,value=T)
    )
gene.vs.all_inPromoter_tp73ConfirmAny.rownames.interest <- unique(gene.vs.all_inPromoter_tp73ConfirmAny.rownames.interest)
    
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

pdf("heatmap_gene_selection_tp73ConfirmAny.pdf", width=18, height=8)
n <- gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant$Gene
heatmap.2(as.matrix(gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant[,-c(1)]),
          trace="none",
          col = colorRampPalette(c("white", "gray50", "black"))(100),
          scale="none", margins=c(10,10),
          key.title="Expression", key.xlab="scaled CUT&RUN and JASPAR binding scores", key.ylab="",
          dendrogram="none", Rowv=FALSE, Colv=FALSE,
          main="Heatmap of Gene Selection for EMT in TP73 Confirmed Contexts",
          labRow=n,
          labCol=colnames(gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant)[-c(1)],
          cexRow=0.75, cexCol=0.75
)
dev.off()

# Plot correlation matrix of all transcription factors (columns) against each other

require(corrplot)

# Use the context data for all promoters with TP73 confirmation as the matrix
cor_matrix <- cor(gene.selection_tp73ConfirmAny_context_interest_heatmap_input_nonRedundant[,-1], use="pairwise.complete.obs", method="pearson")

pdf("correlation_matrix_TFBS.pdf", width=12, height=10)
corrplot(cor_matrix, method="color", type="upper", order="hclust",
         tl.col="black", tl.cex=0.7, addCoef.col="black", number.cex=0.5,
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

