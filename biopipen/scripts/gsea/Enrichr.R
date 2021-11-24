source("{{biopipen_dir}}/utils/io.R")
source("{{biopipen_dir}}/utils/gene.R")
source("{{biopipen_dir}}/utils/gsea.R")

infile = {{in.infile | quote}}
outdir = {{out.outdir | quote}}
genecol = {{envs.genecol | r}}
genename = {{envs.genename | r}}
dbs = {{envs.dbs | r}}
inopts = {{envs.inopts | r}}

if (is.integer(genecol)) {
    genecol = genecol + 1
}

indata = read.table.opts(infile, inopts)
genes = indata[, genecol]

if (genename != "symbol") {
    genedf = gene_name_conversion(
        genes,
        NULL,
        infmt = genename,
        outfmt = "symbol",
        notfound = "skip"
    )

    genes = genedf$symbol
}

runEnrichr(genes, dbs, outdir)
