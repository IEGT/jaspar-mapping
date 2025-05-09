#!/usr/bin/R

require("xlsx")

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


m.findings <- list()
m.contexts <- list()
promoter.for.gene.groups.500 <- list()
chromosomes <- c(as.character(1:22)) #,"X","Y")
#chromosomes <- c(as.character(13:22),"X","Y")

keep.all.data <- FALSE

for(i in chromosomes) {
    
    cat("I: processing chromosome ",i,"...\n",sep="")
    m <- read.data.table.for.chromosome(i)

    if (is.null(m)) {
        cat("E: No data for chromosome ",i," - skipping\n",sep="")
        next
    }

    if ("Feature" %in% colnames(m)) {
        m.colnames.which <- which(colnames(m)=="Feature")
        colnames(m)[m.colnames.which] <- "Name"
    }

    promoter.for.gene.groups.500[[i]] <- m.promoters.index <- lapply(promoterBedTables, function(x) {
        checkBedOverlaps(x,m)
    })

    if (keep.all.data) {
        m.contexts[[i]] <- m
    
        for (j in names(m.promoters.index)) {
            cat("  ", j, ": ", length(m.promoters.index[[j]]), "\n", sep = "")
            print(m[m.promoters.index[[j]],c("Chr","From","To")])
        }
    } else {
        if (is.null(m.contexts[[i]])) {
            m.contexts[[i]] <- list()
        }
        for (j in names(m.promoters.index)) {
            if (is.null(m.contexts[[i]][[j]])) {
                m.contexts[[i]][[j]] <- m[m.promoters.index[[j]],]
            } else {
                m.contexts[[i]][[j]] <- rbind(m.contexts[[i]][[j]],m[m.promoters.index[[j]],])
            }
        }
    }
    gc()
}


# concatenate all gene-groups' contexts
m.contexts.all <- list()
for(i in chromosomes) {
    #for(groupname in  grep(x=names(promoterBedTables),pattern="^Nico",value=T))
    for(groupname in  names(promoterBedTables)) {
        if (is.null(m.contexts.all[[groupname]])) {
            m.contexts.all[[groupname]] <-list()
        }
        if (is.null(m.contexts[[i]][[groupname]]) || nrow(m.contexts[[i]][[groupname]])==0) {
            #cat("I: No data for ",groupname," in chromosome ",i," - skipping\n",sep="")
            next
        }

        if ("Feature" %in% colnames(m.contexts[[i]][[groupname]])) {
            cat("W: Found 'Feature' for chromosome ",i," - renaming to 'Name'\n",sep="")
            m.colnames.which <- which(colnames(m.contexts[[i]][[groupname]])=="Feature")
            colnames(m.contexts[[i]][[groupname]])[m.colnames.which] <- "Name"
        }

        if (is.null(m.contexts.all[[groupname]])) {
            m.contexts.all[[groupname]] <- m.contexts[[i]][[groupname]]
        }
        
        if (! all(colnames(m.contexts[[i]][[groupname]])==colnames(m.contexts.all[[groupname]]))) {
            cat("E: Column names do not match for ",groupname," in chromosome ",i," - skipping\n",sep="")
            next
        }
        #cat("I: concatenating ",groupname," for chromosome ",i,"...\n",sep="")
            #cat("I: ",dim(m.contexts.all[[groupname]]),"\n")
            #cat("I: ",dim(m.contexts[[i]][[groupname]]),"\n")
            #cat("I: ",dim(rbind(m.contexts.all[[groupname]],m.contexts[[i]][[groupname]])),"\n")
        cat("I: concatenating ",groupname," for chromosome ",i," - ", nrow(m.contexts[[i]][[groupname]])," rows...\n",sep="")
        m.contexts.all[[groupname]] <- rbind(m.contexts.all[[groupname]],m.contexts[[i]][[groupname]])
    }
}

#colnames.ref <- colnames(m.contexts.all[["Nico_Analysis_TA_20250508.promoter.bed"]])
#colnames.prob <- colnames(m.contexts[["10"]][["Nico_Analysis_TA_20250508.promoter.bed"]])
#colnames.prob[! colnames.prob%in%colnames.ref]
#colnames.ref[! colnames.ref%in%colnames.prob]

length(colnames(m.contexts[["10"]][["Nico_Analysis_TA_20250508.promoter.bed"]]))
length(colnames(m.contexts.all[["Nico_Analysis_TA_20250508.promoter.bed"]]))

cols.of.interest <- c(1,2,3,5,10,11,12,16,17,18,6855,6856) #,6857,6858)

for(groupname in names(m.contexts.all) ) {
    cat("I: ",groupname," :\n")
    print(m.contexts.all[[groupname]][,..cols.of.interest],"\n",sep="")
}
for(groupname in grep(x=names(m.contexts.all),"Nico",value=T)) {
    cat("I: ",groupname," :\n")
    for(i in chromosomes) {
        cat("Chromosome ",i," :\n")
        print(m.contexts[[i]][[groupname]][,1:18],"\n",sep="")
    }
}


scores <- sapply(m.contexts.all, function(X) {
    cols.Score <- grepl("_Score", colnames(X))
    apply(X[,..cols.Score],2,function(Y) {
        sum(!is.na(Y))
    })})

selective.for.TA <- scores[scores[,"Nico_Analysis_DN_20250508.promoter.bed"]==0 & scores[,"Nico_Analysis_TA_20250508.promoter.bed"]>=7,]
selective.for.DN <- scores[scores[,"Nico_Analysis_DN_20250508.promoter.bed"]>=2 & scores[,"Nico_Analysis_TA_20250508.promoter.bed"]<=1,]
write.xlsx(selective.for.TA, file="selective_for_TA.xlsx", sheetName="selective_for_TA", append=FALSE)
write.xlsx(selective.for.DN, file="selective_for_DN.xlsx", sheetName="selective_for_DN", append=FALSE)

for(groupname in 1:ncol(selective.for.TA) ) {
 cat(colnames(selective.for.TA)[groupname],"\n mean TA: ",
     mean(selective.for.TA[,groupname]),"  mean DN: ",
     mean(selective.for.DN[,groupname]),"\n",sep="")
 print(wilcox.test(selective.for.TA[,groupname], selective.for.DN[,groupname], alternative="two.sided", paired=FALSE,exact=FALSE))
}

candidate.cofactors <- c("HEY1_MA0823.1","ABF1_MA0570.2","br_MA0010.1","MYB17_MA1766.1","MYB80_MA1778.1",
                        "TDA9_MA0431.1","RORA_MA0071.1","OLIG1_MA0826.1","Ar_MA0007.3")
candidate.cofactors.shift <- paste(candidate.cofactors,"Shift",sep="_")


shifts.sd <- lapply(m.contexts.all, function(X) {
    cols.Shift <- grepl("_Shift", colnames(X))
    apply(X[,..cols.Shift],2,function(Y) {
        sd(Y,na.rm=TRUE)
    })})

shifts.median <- lapply(m.contexts.all, function(X) {
    cols.Shift <- grepl("_Shift", colnames(X))
    apply(X[,..cols.Shift],2,function(Y) {
        median(Y,na.rm=TRUE)
    })})

shifts.mean <- lapply(m.contexts.all, function(X) {
    cols.Shift <- grepl("_Shift", colnames(X))
    apply(X[,..cols.Shift],2,function(Y) {
        mean(Y,na.rm=TRUE)
    })})

shifts.coefficientOfVariatin <- lapply(m.contexts.all, function(X) {
    cols.Shift <- grepl("_Shift", colnames(X))
    apply(X[,..cols.Shift],2,function(Y) {
        sd(Y,na.rm=TRUE)/mean(Y,na.rm=TRUE)
    })})

for(groupname in names(m.contexts.all) ) {
    cat("I: ",groupname," :\n")
    cat("   sd\n")
    print(shifts.sd[[groupname]][candidate.cofactors.shift])
    cat("   coefficient of variation\n")
    print(shifts.coefficientOfVariatin[[groupname]][candidate.cofactors.shift])
    cat("   median\n")
    print(shifts.median[[groupname]][candidate.cofactors.shift])
    cat("   mean\n")
    print(shifts.mean[[groupname]][candidate.cofactors.shift])
}


if (false) 
    m <- create.lists.for.chromosome(m)
    m.findings[[i]] <- attributes(m)
    m.findings[[i]]$promoterBedTables <- m.promoters.index
    rm(m)
}

# Identify most promising transcription factors for 90% quantile
m.all.findings <- NULL
#lists <- c("ratio.skmel29_2.TAa.99.vs.ratio.skmel29_2.DNb.90","ratio.saos2.TAa.99.vs.ratio.saos2.DNb.90")
#lists <- c("ratio.skmel29_2.TAa.99.vs.ratio.skmel29_2.DNb.90","ratio.saos2.TAa.99.vs.ratio.saos2.DNb.90")
lists <- c("ratio.skmel29_2.TAa.quantile.90.vs.ratio.skmel29_2.DNb.quantile.90","ratio.saos2.TAa.99.vs.ratio.saos2.DNb.90")
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
            "E2F1_MA0024.3", "E2F2_MA0864.2", "E2F4_MA0470.2", "E2F6_MA0471.2"
            "FOS_MA1951.1","FOSL1--JUN_MA1129.1",
            "FOXO4_MA0848.1","Foxq1_MA0040.1", "KLF14_MA0740.2","Klf15_MA1890.1","REL_MA0101.1", "PLAG1_MA0163.1","POU2F1--SOX2_MA1962.1", "PPARG_MA0066.1",
            "RELA_MA0107.1","REST_MA0138.2", "RXRA--VDR_MA0074.1", "SP1_MA0079.5", "TBP_MA0108.2","YY1-2_MA1927.1",
)

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
