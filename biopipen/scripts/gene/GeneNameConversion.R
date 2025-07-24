library(biopipen.utils)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
notfound <- {{envs.notfound | r}}
genecol <- {{envs.genecol | r}}
output <- {{envs.output | r}}
dup <- {{envs.dup | r}}
infmt <- {{envs.infmt | r}}
outfmt <- {{envs.outfmt | r}}
species <- {{envs.species | r}}

log <- get_logger()

if (is.na(notfound)) {
    notfound = "na"
}

df <- read.table(infile, header=TRUE, sep="\t", check.names=FALSE)

if (genecol == 0) {
    log$warn("envs.genecol should be 1-based, but 0 was given. Using 1 instead.")
    genecol <- 1
}

if (is.numeric(genecol)) { genecol <- colnames(df)[genecol] }
if (dup == "combine") { dup <- ";" }

genes <- df[[genecol]]
converted <- gene_name_conversion(
    genes = genes,
    species = species,
    infmt = infmt,
    outfmt = outfmt,
    notfound = notfound,
    dup = dup,
    suppress_messages = FALSE
)
#    <genecol> <outfmt>
# 1  1255_g_at   GUCA1A
# 2    1316_at     THRA
# 3    1320_at   PTPN21
# 4    1294_at  MIR5193

# order the converted dataframe by the original gene column
converted <- converted[order(match(converted$query, genes)), , drop=FALSE]
outcol <- outfmt

if (notfound == "skip" || notfound == "ignore") {
    df <- df[df[[genecol]] %in% converted$query, , drop=FALSE]
}

if (output == "append") {
    if (outfmt %in% colnames(df)) {
        log$warn("The output column name already exists in the input dataframe. Appending with a suffix `_1`.")
        outcol <- paste(outfmt, "_1", sep="")
    }
    df[[outcol]] <- converted[[outfmt]]
} else if (output == "replace") {
    df[[genecol]] <- converted[[outfmt]]
} else if (output == "with-query") {
    df <- converted
} else {
    df <- converted[, outfmt, drop=FALSE]
}

write.table(df, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
