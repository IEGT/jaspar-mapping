
# Iterate over chromosomes and retrieve context data for DN.enriched.valid and TA.enriched.valid rows
cat("I: Retrieving context data for DN.enriched.valid and TA.enriched.valid rows by iterating over chromosomes...\n")

retrieve_context_data_by_chromosome <- function(enriched_rows=NULL,confirmation=NULL,TA.or.DN=NULL,verbose=TRUE) {
    context_data <- NULL
    context_data_binary <- NULL
    context_shifts <- NULL
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

    cols.NumInWindow <- grepl("_NumInWindow$", colnames(m.contexts[[1]]))
    cols.Shift <- grepl("__Shift$", colnames(m.contexts[[1]]))

    # Create a matrix to store column sums for each chromosome
    mean_by_chromosome <- col_sums_by_chromosome <- col_sums_by_chromosome_binary <- matrix(
        NA,
        nrow = length(m.contexts),
        ncol = sum(cols.NumInWindow),
        dimnames = list(names(m.contexts), colnames(m.contexts[[1]][,..cols.NumInWindow]))
    )
    num_matches_by_chromosome <- setNames(vector("numeric", length(m.contexts)), names(m.contexts))

    for (chromosome in names(m.contexts)) {

        context_data_chromosome <- NULL
        context_data_chromosome_binary <- NULL
        context_shifts_chromosome <- NULL
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
            # FIXME
            cat("W: Need to revise the logic for enriched_rows==NULL for downstream genes.\n")
            downstream_genes_chromosome <- c()
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
        context_data_chromosome_binary <- m.contexts[[chromosome]][context_matches_chromosome, ..cols.NumInWindow]>0
        context_shifts_chromosome <- m.contexts[[chromosome]][context_matches_chromosome, ..cols.Shift]
        cutandrun_data_chromosome <- m.contexts[[chromosome]][context_matches_chromosome, 7:18, drop = FALSE]
        coordinates_chromosome <- m.contexts[[chromosome]][context_matches_chromosome, 1:6, drop = FALSE]

        cat("I: Retrieved context for ", nrow(context_data_chromosome), " binding sites on chromosome: ", chromosome, "\n", sep = "")


        if (nrow(context_data_chromosome) == 0) {
            cat("E: No context data found for chromosome ", chromosome, " - skipping\n", sep = "")
            next
        }
        else {
            context_data <- rbind(context_data, context_data_chromosome)
            context_data_binary <- rbind(context_data_binary, context_data_chromosome_binary)
            context_shifts <- rbind(context_shifts, context_shifts_chromosome)
            context_matches <- c(context_matches, context_matches_chromosome)
            cutandrun_data <- rbind(cutandrun_data, cutandrun_data_chromosome)
            coordinates <- rbind(coordinates, coordinates_chromosome)
            downstream_genes <- c(downstream_genes, downstream_genes_chromosome)
            col_sums_by_chromosome[chromosome,] <- colSums(context_data_chromosome, na.rm = TRUE)
            col_sums_by_chromosome_binary[chromosome,] <- colSums(context_data_chromosome_binary, na.rm = TRUE)
            num_matches_by_chromosome[chromosome] <- nrow(context_data_chromosome)
            mean_by_chromosome[chromosome,] <- colMeans(context_data_chromosome, na.rm = TRUE)
            #mean_shifts_by_chromosome[chromosome,] <- colMeans(context_shifts_chromosome, na.rm = TRUE)
            cat("I: Calculated colSums, number of matches, and mean for chromosome ", chromosome, "\n", sep = "")
        }
    }
    col_sums_total = colSums(col_sums_by_chromosome)
    col_sums_binary_total = colSums(col_sums_by_chromosome_binary)
    return(list(
        # Number of matches within the context data for each binding site
        context_data = context_data,
        # Indication if more than 0 reads were found in the context data for each binding site
        context_data_binary = context_data_binary,
        # Shifts of the cofactors' binding sites in the context data (e.g. for TP73)
        context_shifts = context_shifts,
        # FIXME: Position in gene list passed to filter?
        context_matches = context_matches,
        # FIXME: broken
        downstream_genes = downstream_genes,
        # Columns with cut&run data
        cutandrun_data = cutandrun_data,
        # Genetic coordinates of the binding sites (Chr, From, To, Gene (TP73), Score, Strand)
        coordinates = coordinates,
        # Sum of binding sites within the context data for each factor (multiple occurences in same context counted) per chromosome
        col_sums_by_chromosome = col_sums_by_chromosome,
        # Total sum of binding sites within the context data for each factor (multiple occurences in same context counted)
        col_sums_total = col_sums_total,
        # Sum of binding sites within the context data for each factor (no redundant bindings per binding site counted)
        col_sums_by_chromosome_binary = col_sums_by_chromosome_binary,
        # Total sum of binding sites within the context data for each factor (no redundant bindings per binding site counted)
        col_sums_binary_total = col_sums_binary_total,
        # Number of binding sites per chromosome after filter is applied
        num_matches_by_chromosome = num_matches_by_chromosome,
        # Total number of binding sites across all chromosomes after filter is applied
        num_matches_total = sum(num_matches_by_chromosome),
        # Matrix with mean number of occurrences of each factor in the context data for each chromosome, multiple occurrences in the same context counted
        mean_by_chromosome = mean_by_chromosome,
        # Mean number of occurrences of each factor in the context data for each chromosome
        mean_total = colMeans(context_data, na.rm = TRUE),
        # Mean number of occurrences of each factor in the context data across all chromosomes, counting factor only once per binding
        mean_total_binary = colMeans(context_data_binary, na.rm = TRUE)
    ))
}