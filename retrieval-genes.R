# Retrieval of data from context for c(chromosome, min, max) triplet.
library(data.table)
gene.positions <- data.table(
    gene       = c("IL10",    "TGFB1",  "CD274"),
    chromosome = c("1",       "19",     "9"),
    start      = c(206767602, 41301587, 5450503),
    end        = c(206774541, 41353961, 5470566),
    strand     = c("-",       "-",       "+")
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
    print(gene_region)
    promoter_region <- gene_region_to_promoter_region(gene_region, upstream=upstream, intragenic=intragenic)
    r <- retrieve_p73_scores_from_contexts(promoter_region, columns.of.interest)
    description <- genomicRegion2description(gene_region)
    list("description"=description,"gene"=gene_region,"promoter"=promoter_region,"details"=r)
}

regions_to_retrieve <- lapply(c("IL10", "TGFB1", "CD274"), retrieve_p73_scores_for_gene, upstream=5000, intragenic=5000)

l <- lapply(regions_to_retrieve, function(x) {
    x[["details"]]
})
l.desc <- sapply(regions_to_retrieve, function(x) {
    x[["description"]]
})
names(l) <- l.desc

#require(write.xlsx)
require(openxlsx)
openxlsx::write.xlsx(l, file="p73_binding_sites_IL10_TGFB1_CD274.xlsx", asTable=TRUE)