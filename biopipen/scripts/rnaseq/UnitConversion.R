library(rlang)
library(glue)
library(biopipen.utils)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
inunit <- {{envs.inunit | r}}
outunit <- {{envs.outunit | r}}
refexon <- {{envs.refexon | r}}
meanfl <- {{envs.meanfl | r}}
nreads <- {{envs.nreads | r}}

log <- get_logger()

log$info("Reading input data ...")
indata = read.table(infile, header = TRUE, sep = "\t", row.names = 1, check.names = F)
samples = colnames(indata)

# parse the inunit to see if there is any transformation
parsable <- function(arg) {	is.call(arg) || is_symbol(arg) }

check_call_args <- function(arg1, arg2) {
    if (parsable(arg1) && parsable(arg2)) {
        stop(glue("Can't parse the call. Multiple names or calls detected: {arg1}, {arg2}\n"))
    }
    if (!parsable(arg1) && !parsable(arg2)) {
        stop(glue("Can't parse the call. Both arguments are constants: {arg1}, {arg2}. Use the result directly\n"))
    }
}

parse_call <- function(call, expr = "indata") {
    if (!is.call(call)) {
        call <- match.arg(
            as_string(call),
            c(
                "count", "counts", "rawcount", "rawcounts",
                "cpm",
                "fpkm", "rpkm",
                "fpkmuq", "rpkmuq",
                "tpm",
                "tmm"
            )
        )
        return(glue("{as_string(call)} = {expr}"))
    }
    cn <- as_string(call_name(call))
    args <- call_args(call)
    if (length(args) == 1) {
        # This should be those supported functions
        cn <- match.arg(cn, c("log", "log2", "log10", "exp", "sqrt"))
        if (cn == "log") return(parse_call(args[[1]], glue("e ^ ({expr})")))
        if (cn == "log2") return(parse_call(args[[1]], glue("2 ^ ({expr})")))
        if (cn == "log10") return(parse_call(args[[1]], glue("10 ^ ({expr})")))
        if (cn == "exp") return(parse_call(args[[1]], glue("log({expr})")))
        if (cn == "sqrt") return(parse_call(args[[1]], glue("({expr}) ^ 2")))
    } else {
        check_call_args(args[[1]], args[[2]])
        if (cn == "+") {
            if (parsable(args[[1]])) return(parse_call(args[[1]], glue("{expr} - {args[[2]]}")))
            return(parse_call(args[[2]], glue("{expr} - {args[[1]]}")))
        }
        if (cn == "-") {
            if (parsable(args[[1]])) return(parse_call(args[[1]], glue("{expr} + {args[[2]]}")))
            return(parse_call(args[[2]], glue("{args[[1]]} - {expr}")))
        }
        if (cn == "*") {
            if (parsable(args[[1]])) return(parse_call(args[[1]], glue("({expr}) / ({args[[2]]})")))
            return(parse_call(args[[2]], glue("({expr}) / ({args[[1]]})")))
        }
        if (cn == "/") {
            if (parsable(args[[1]])) return(parse_call(args[[1]], glue("({expr}) * ({args[[2]]})")))
            return(parse_call(args[[2]], glue("({args[[1]]}) / ({expr})")))
        }
        if (cn == "^") {
            if (parsable(args[[1]])) return(parse_call(args[[1]], glue("{expr} * (1 / ({args[[2]]}))")))
            return(parse_call(args[[2]], glue("log({expr}, {args[[1]]})")))
        }
        stop(paste0("Unknown function to parse: {cn}\n"))
    }
}

glenFromExon = function(exonfile, data) {
    gff  = read.table(exonfile, header = F, row.names = NULL)
    # V4: start, V5: end, V10: gene name
    glen = aggregate(V5-V4+1 ~ V10, gff, sum)
    genes = glen[,1]
    glen = glen[,-1,drop=F]
    rownames(glen) = genes

    mygenes  = rownames(data)
    outgenes = intersect(genes, mygenes)
    if (length(outgenes) < length(mygenes))
        logger('Genes not found in refexon: ', paste(setdiff(mygenes, outgenes), collapse = ','), level = 'WARNING')

    glen[outgenes, , drop = FALSE]
}

meanflFromFile = function(samples, mflfile) {
    if (is.numeric(mflfile)) {
        ret = matrix(mflfile, nrow = length(samples), ncol = 1)
        rownames(ret) = samples
    } else {
        ret = read.table(mflfile, header = F, row.names = 1, check.names = F, sep = "\t")
        ret = ret[samples,,drop = F]
    }
    ret
}

nreadsFromFile = function(samples, nreads) {
    if (is.numeric(nreads)) {
        ret = matrix(nreads, nrow = length(samples), ncol = 1)
        rownames(ret) = samples
    } else {
        ret = read.table(nreads, header = F, row.names = 1, check.names = F, sep = "\t")
        ret = ret[samples,,drop = F]
    }
    ret
}

count2cpm <- function(data) {
    edgeR::cpm(data)
}

count2fpkm = function(data) {
    # may lose some genes
    glen = glenFromExon(refexon, data)
    data = data[rownames(glen), , drop = F]
    dge  = edgeR::DGEList(counts=data)

    dge$genes$Length = glen
    edgeR::rpkm(dge)
}

count2fpkmuq = function(data) {
    # may lose some genes
    glen = glenFromExon(refexon, data)
    data = data[rownames(glen), , drop = FALSE]

    fld  = meanflFromFile(samples, meanfl)
    expr = sapply(samples, function(s){
        RC75 = quantile(data[, s], .75)
        exp( log(data[, s]) + log(1e9) - log(glen - fld[s, ] + 1) - log(RC75) )
    })
    rownames(expr) = rownames(data)
    expr
}

count2tpm = function(data) {
    glen = glenFromExon(refexon, data)
    data = data[rownames(glen), , drop = F]
    fld  = meanflFromFile(samples, meanfl)

    # see: https://gist.github.com/slowkow/c6ab0348747f86e2748b
    expr = as.data.frame(sapply(samples, function(s){
        rate  = log(data[, s]) - log(glen - fld[s, ] + 1)
        denom = log(sum(exp(rate)))
        exp(rate - denom + log(1e6))
    }))
    colnames(expr) = colnames(data)
    rownames(expr) = rownames(data)
    expr
}

count2tmm = function(data) {
    dge = edgeR::DGEList(counts=data)
    dge = edgeR::calcNormFactors(dge, method = "TMM")
    edgeR::cpm(dge)
}

fpkm2count = function(data) {
    glen    = glenFromExon(refexon, data)
    data    = data[rownames(glen), , drop = F]
    fld     = meanflFromFile(samples, meanfl)
    totalnr = nreadsFromFile(samples, nreads)

    expr    = sapply(samples, function(s){
        N = totalnr[s, ]
        exp( log(data[, s]) + log(N) + log(glen - fld[s, ] + 1) - log(1e9) )
    })
    rownames(expr) = rownames(data)
    expr
}

fpkm2tpm = function(data) {
    expr = sapply(samples, function(s) {
        exp( log(data[, s]) - log(sum(data[, s])) + log(1e6) )
    })
    rownames(expr) = rownames(data)
    expr
}

fpkm2cpm = function(data) {
    glen = glenFromExon(refexon, data)
    data = data[rownames(glen), , drop = F]
    expr = sapply(samples, function(s) {
        exp( log(data[, s]) - log(1e3) - log(glen - fld[s, ] + 1) )
    })
    rownames(expr) = rownames(data)
    expr
}

tpm2count = function(data) {
    totalnr = nreadsFromFile(samples, nreads)
    ngenes  = nrow(data)

    expr     = sapply(samples, function(s){
        # counts to tpm:
        # rate <- log(counts) - log(effLen)
        # denom <- log(sum(exp(rate)))
        # tpm = exp(rate - denom + log(1e6))
        # so:
        # log(tpm) = rate - denom + log(1e6)
        # rate = log(tpm) + denom - log(1e6)
        # log(counts) - log(effLen) = log(tpm) + log(sum(exp(rate))) - log(1e6)
        # log(counts) - log(effLen) = log(tpm) + log(sum(exp(log(counts) - log(effLen)))) - log(1e6)
        # log(counts) - log(effLen) = log(tpm) + log(sum(exp(log(counts))/exp(log(effLen)))) - log(1e6)
        # log(counts) - log(effLen) = log(tpm) + log(sum(counts/effLen)) - log(1e6)
        #                                                ?????????????
        # ??? estimated by sum(counts)/sum(effLen) * length(effLen)
        # log(counts) = log(effLen) + log(tpm) + log(sum(counts)) - log(effLen) + log(length(effLen))) - log(1e6)
        # counts = expr( log(tpm) + log(nreads) + log(length(effLen)) - log(1e6) )
        exp( log(data[, s]) + log(totalnr[s, ]) + log(ngenes) - log(1e6) )
    })
    rownames(expr) = rownames(data)
    expr
}

tpm2fpkm = function(data) {
    totalnr = nreadsFromFile(samples, nreads)
    expr = sapply(samples, function(s) {
        exp( log(data[, s]) - log(1e6) + log(totalnr[s, ]) )
    })
    rownames(expr) = rownames(data)
    expr
}

tpm2cpm = function(data) {
    glen   = glenFromExon(refexon, data)
    data   = data[rownames(glen), , drop = F]
    fld    = meanflFromFile(samples, meanfl)
    ngenes = length(outgenes)

    expr = sapply(samples, function(s) {
        exp( log(data[, s]) + log(glen - fld[s, ] + 1) - log(sum(glen - fld[s, ] + 1)) + log(ngenes) )
    })
    rownames(expr) = rownames(data)
    expr
}

cpm2count = function(data) {
    totalnr = nreadsFromFile(samples, nreads)

    expr = sapply(samples, function(s) {
        exp( log(data[, s]) + log(totalnr[s, ]) - log(1e6) )
    })
    rownames(expr) = rownames(data)
    expr
}

cpm2fpkm = function(data) {
    glen = glenFromExon(refexon, data)
    data = data[rownames(glen), , drop = F]
    expr = sapply(samples, function(s) {
        exp( log(data[, s]) + log(1e3) - log(glen - fld[s, ] + 1) )
    })
    rownames(expr) = rownames(data)
    expr
}

cpm2tpm = function(data) {
    glen   = glenFromExon(refexon, data)
    data   = data[rownames(glen), , drop = F]
    ngenes = nrow(glen)
    expr   = sapply(samples, function(s) {
        exp( log(data[, s]) - log(glen - fld[s, ] + 1) - log(sum(glen - fld[s, ] + 1)) + log(ngenes) )
    })
    rownames(expr) = rownames(data)
    expr
}

is.count  = function(unit) {unit %in% c('count', 'counts', 'rawcount', 'rawcounts')}
is.cpm    = function(unit) {unit == 'cpm'}
is.fpkm   = function(unit) {unit %in% c('fpkm', 'rpkm')}
is.fpkmuq = function(unit) {unit %in% c('fpkmuq', 'rpkmuq')}
is.tpm    = function(unit) {unit == 'tpm'}
is.tmm    = function(unit) {unit == 'tmm'}

# log2(count + 1) -> count = 2 ^ indata - 1
parsed_transformation <- parse_call(parse_expr(inunit))
splits <- strsplit(parsed_transformation, " = ")[[1]]
if (is.count(splits[[1]])) {
    intype <- "count"
} else if (is.cpm(splits[[1]])) {
    intype <- "cpm"
} else if (is.fpkm(splits[[1]])) {
    intype <- "fpkm"
} else if (is.fpkmuq(splits[[1]])) {
    intype <- "fpkmuq"
} else if (is.tpm(splits[[1]])) {
    intype <- "tpm"
} else if (is.tmm(splits[[1]])) {
    intype <- "tmm"
} else {
    stop(glue("Can't find a supported unit in the inunit: {inunit}\n"))
}
splits[1] <- intype
eval(parse_expr(paste(splits, collapse = " = ")))
indata <- get(intype)

# find out the outtype
if (grepl('rawcounts|rawcount|counts|count', outunit)) {
    outtype <- 'count'
    outunit <- gsub('rawcounts|rawcount|counts|count', 'count', outunit)
} else if (grepl('fpkmuq|rpkmuq', outunit)) {
    outtype <- 'fpkmuq'
    outunit <- gsub('fpkmuq|rpkmuq', 'fpkmuq', outunit)
} else if (grepl('fpkm|rpkm', outunit)) {
    outtype <- 'fpkm'
    outunit <- gsub('fpkm|rpkm', 'fpkm', outunit)
} else if (grepl('tpm', outunit)) {
    outtype <- 'tpm'
} else if (grepl('cpm', outunit)) {
    outtype <- 'cpm'
} else if (grepl('tmm', outunit)) {
    outtype <- 'tmm'
} else {
    stop(glue("Can't find a supported unit in the outunit: {outunit}\n"))
}

log$info("Transforming data by resolving {inunit} ...")
if (intype == outtype) {
    fun <- identity
} else {
    fun <- glue("{intype}2{outtype}")
    fun <- tryCatch(
        { get(fun) },
        error = function(e) { stop(glue("Unsupported conversion from {intype} to {outunit}\n")) }
    )
}
assign(outtype, fun(indata))
out <- eval(parse_expr(outunit))

log$info("Saving output data ...")
write.table(out, outfile, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
