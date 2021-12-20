
.normUnit = function(unit) {
    if ("count" %in% unit) {
        return("count")
    }
    return(unit)
}

glenFromGFFExons = function(exonfile) {
    gff  = read.table(exonfile, header = F, row.names = NULL)
    # V4: start, V5: end, V10: gene name
    glen = aggregate(V5-V4+1 ~ V10, gff, sum)
    genes = glen[,1]
    glen = glen[,-1,drop=TRUE]
    names(glen) = genes
    return(glen)
}

count2tpm = function(x, args) {
    if (is.null(args$genelen)) {
        stop("Gene lengths are required to convert count to TPM.")
    }
    glengenes = names(args$genelen)
    mygenes = rownames(x)
    missing = setdiff(mygenes, glengenes)
    warning(paste(length(missing), "gene cannot be found in gene length data"))
    warning(paste(missing, sep=", "))

    genes = intersect(mygenes, glengenes)
    x = x[genes, , drop=FALSE]

    # see: https://gist.github.com/slowkow/c6ab0348747f86e2748b
    # and https://support.bioconductor.org/p/91218/
    out = x / unlist(args$genelen[genes])
    out = t(t(out) * 1e6 / colSums(out))
    rownames(out) = genes
    colnames(out) = colnames(x)

    return(out)
}


unit_conversion = function(x, inunit, outunit, args=list()) {
    inunit = .normUnit(inunit)
    outunit = .normUnit(outunit)
    func = get(paste0(inunit, "2", outunit))
    func(x, args)
}
