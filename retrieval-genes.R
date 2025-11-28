# Retrieval of data from context for c(chromosome, min, max) triplet.
library(data.table)
gene.positions <- data.table(
    group      = c("Max",     "Max",    "ImmuneCheckpoint",  "ImmuneCheckpoint", "ImmuneCheckpoint", "Max",    "Max",    "Max",    "p53",   "p53",     "p53",   "p53",   "Interact","Yamanaka","Yamanaka","Yamanaka","Yamanaka","Target"),  
    gene       = c("IL10",    "TGFB1",  "CD274",             "PDCD1LG2",         "INCR1",            "SP1",    "PATZ1",   "CD109",  "IGFBP4", "LRRC32", "TP53",  "TP63",    "TP73", "MDM2",  "REST",    "POU5F1",  "SOX2",    "KLF4",    "MYC",     "GABBR2"),
    chromosome = c("1",       "19",     "9",                 "9",                "9",                "12",     "22",      "6",      "17",     "11",     "17",    "3",       "1",    "12",    "4",       "6",       "3",       "9",       "8",       "9"),
    start      = c(206767602, 41301587, 5450503,             5510491,            5457477,            53380176, 31325804,  73695785, 40443446, 76657524, 7661779, 189631389, 3652516, 68808175, 56907876,  31164337,  181711925, 107484852, 127735434, 98288096),
    end        = c(206774541, 41353961, 5470566,             5571282,            5629780,            53416446, 31346605,  73828316, 40457727, 76657868, 7687546, 189897276, 3736201, 68845544, 56966808,  31180731,  181714436, 107490482, 127742951, 98709739),
    strand     = c("-",       "-",       "+",                "+",                "-",                "+",       "-",       "+",   "+",      "-",      "-",     "+",       "+",     "+",      "+",       "-",       "+",       "-",       "+",       "-")
)

/* Convert gene region to promoter region.
   upstream: distance upstream of TSS to include in promoter region.
   intragenic: distance downstream of TSS to include in promoter region.
*/
gene_region_to_promoter_region  <-  function(gene_region, upstream = 2000, intragenic = 2000) {

    r <- list()
    r[["chromosome"]] <- gene_region$chromosome
    r[["strand"]] <- gene_region$strand
    r[["name"]] <- gene_region$name

    if (gene_region$strand == "+") {
        r[["gene_start"]] <- min(gene_region$start, gene_region$end)
        r[["start"]] <- r[["gene_start"]] - upstream
        r[["end"]] <- r[["gene_start"]] + intragenic
    } else if (gene_region$strand == "-") {
        r[["gene_start"]] <- max(gene_region$start, gene_region$end)
        r[["start"]] <- r[["gene_start"]] + upstream
        r[["end"]] <- r[["gene_start"]] - intragenic
    } else {
        stop("Strand information is invalid.")
    }

    return(r)
}   


/* Example usage: Retrieve gene regions for IL10, TGFB1, and CD274 */
get_gene_region <- function(gene_name) {
    gene_info <- gene.positions[gene == gene_name]
    if (nrow(gene_info) == 0) {
        stop("Gene not found in the database.")
    }
    #print(gene_info)
    return(list(
        name = gene_name,
        chromosome = gene_info$chromosome,
        start = gene_info$start,
        end = gene_info$end,
        strand = gene_info$strand
    ))
}

retrieve_p73_scores_from_contexts <- function(genomic_region, columns.of.interest=1:18) {

    # Variant presuming global variable m.contexts
    if (!exists("m.contexts")) {
        stop("Global variable m.contexts does not exist.")
    }
    if (! genomic_region$chromosome %in% names(m.contexts)) {
        stop(paste0("Chromosome '", genomic_region$chromosome, "' not found in m.contexts."))
    }

    p73.matches <- m.contexts[[genomic_region$chromosome]]
    cat("I: Retrieving p73 scores for chromosome:", genomic_region$chromosome, "with",nrow(p73.matches)," binding sites.\n")
    #p73.matches <- m.contexts[["1"]]
    eligible <- p73.matches$From >= min(genomic_region$start,genomic_region$end) & p73.matches$To <= max(genomic_region$start,genomic_region$end)
    cat("I: Found",sum(eligible)," binding sites in region",genomic_region$chromosome,":",genomic_region$start,"-",genomic_region$end,"\n")
    lines.of.interest <- p73.matches[eligible, ..columns.of.interest]
    return(lines.of.interest)

}

genomicRegion2description <- function(genomic_region) {
    paste0(genomic_region$name, " ", genomic_region$chromosome, " ", genomic_region$start, "-", genomic_region$end, " ", genomic_region$strand, "")
}

retrieve_p73_scores_for_gene <- function(gene_name, columns.of.interest=1:18, upstream=2000, intragenic=2000) {
    gene_region <- get_gene_region(gene_name)
    if (is.null(gene_region)) {
        stop(paste0("Gene '", gene_name, "' not found."))
    }
    if (is.na(intragenic) || "intragenic" == intragenic) {
        intragenic <- abs(gene_region$end - gene_region$start)+1
    }
    promoter_region <- gene_region_to_promoter_region(gene_region, upstream=upstream, intragenic=intragenic)
    r <- retrieve_p73_scores_from_contexts(promoter_region, columns.of.interest)
    location <- ifelse(r$From >= gene_region$start & r$To <= gene_region$end, "Intragenic", "Promoter")
    r <- cbind(Location=location, r)
    description <- genomicRegion2description(gene_region)
    list("description"=description,"gene"=gene_region,"promoter"=promoter_region,"details"=r)
}

time <- format(Sys.time(),"%Y%m%d")
for (group in unique(gene.positions$group)) {
    genes_in_group <- gene.positions$gene[gene.positions$group == group]
    regions_to_retrieve <- lapply(genes_in_group,
                                retrieve_p73_scores_for_gene, upstream=5000, intragenic="intragenic")

    l <- lapply(regions_to_retrieve, function(x) {
        x[["details"]]
    })
    l.desc <- sapply(regions_to_retrieve, function(x) {
        x[["description"]]
    })
    names(l) <- l.desc

    require(openxlsx)
    f<-paste0("p73_binding_sites_intragenic_group_",group,"_",time,".xlsx")
    openxlsx::write.xlsx(l, file=f, asTable=TRUE)
    cat(Sys.info()["user"],"@",Sys.info()["nodename"],":",getwd(),"/",f,sep="","\n",file=stdout())
}