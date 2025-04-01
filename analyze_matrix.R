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

for(i in as.character(1:22)) {
    cat("I: processing chromosome ",i,"...\n",sep="")
    m <- read.data.table.for.chromosome(i)
    m <- create.lists.for.chromosome(m,i)
}

l <- list("TAa"=sum.cutandrun.tp73.TAa>0, "DNb"=sum.cutandrun.tp73.DNb>0, "GFP"=sum.cutandrun.tp73.GFP>0,
          "TAa_without_DNb"=sum.cutandrun.tp73.TAa>=sum.cutandrun.tp73.TAa.quantile.90 & 
                            sum.cutandrun.tp73.DNb<=sum.cutandrun.tp73.DNb.quantile.50,
          "DNb_without_TAa"=sum.cutandrun.tp73.DNb>=sum.cutandrun.tp73.DNb.quantile.90 & 
                            sum.cutandrun.tp73.TAa<=sum.cutandrun.tp73.TAa.quantile.50)

for (l.name in names(l)) {
    for(tf in c("TP73_MA0861.1","TP63_MA0525.2","TP53_MA0106.3",
            "CTCF_MA0139.1","E2F2_MA0864.2", "FOS_MA1951.1","FOSL1--JUN_MA1129.1",
            "FOXO4_MA0848.1","Foxq1_MA0040.1", "KLF14_MA0740.2","Klf15_MA1890.1","REL_MA0101.1",
            "PLAG1_MA0163.1","POU2F1--SOX2_MA1962.1", "PPARG_MA0066.1",
            "RELA_MA0107.1","REST_MA0138.2", "RXRA--VDR_MA0074.1", "SP1_MA0079.5",
            "TBP_MA0108.2","YY1-2_MA1927.1")) {
        cat("I: plotting '",tf,"' (",l.name,") : ",sep="")
        plotShiftOfBinding(tf=tf,tfbs.selection=l[[l.name]],subset.name=l.name)
        cat("\n")
    }
}
