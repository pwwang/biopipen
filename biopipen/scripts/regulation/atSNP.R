snpinfo2atsnp <- function(snpinfo) {
    # c("chrom", "start", "end", "name", "score", "strand", "ref", "alt", "ref_seq", "alt_seq")
    if (any(nchar(snpinfo$ref) != 1) || any(nchar(snpinfo$alt) != 1)) {
        stop("Only SNVs are supported by atSNP. Consider using motifbreakR instead if you have indels.")
    }
    base_encodings <- c(A = 1, C = 2, G = 3, T = 4)
    transition <- matrix(
        c(
            0.3225035, 0.1738422, 0.24915044, 0.2545039,
            0.3451410, 0.2642147, 0.05245011, 0.3381942,
            0.2813089, 0.2136604, 0.26749171, 0.2375390,
            0.2149776, 0.2071733, 0.25309238, 0.3247568
        ),
        nrow = 4,
        byrow = TRUE
    )
    rownames(transition) <- colnames(transition) <- names(base_encodings)
    list(
        sequence_matrix = unname(sapply(
            snpinfo$ref_seq,
            function(s) as.integer(base_encodings[strsplit(s, "")[[1]]])
        )),
        ref_base = as.integer(base_encodings[snpinfo$ref]),
        snp_base = as.integer(base_encodings[snpinfo$alt]),
        snpids = snpinfo$name,
        transition = transition,
        prior = c(A = 0.287, C = 0.211, G = 0.213, T = 0.289),
        rsid.na = NULL,
        rsid.rm = NULL,
        rsid.duplicate = NULL,
        rsid.missing = NULL
    )
}