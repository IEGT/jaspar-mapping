
prettyIdentifierJaspar <- function(X) {
    A<-strsplit(X,"_")
    B<-sapply(A,function(X) {
        if (2>length(X)) return(X)
        paste(X[1]," (",X[2],")",sep="")
    })
    B
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

# Generate summary for "skmel29_2_any"
nicer.row.names <- function(X) {
    X <- strsplit(X, split="_")
    X <- sapply(X, function(x) {
        if (length(x) >= 2) {
            paste0(x[1]," ","(",x[2],")")
        } else {
            x
        }
    })
}