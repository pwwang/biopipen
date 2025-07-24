{{ biopipen_dir | joinpaths: "utils", "io.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "gene.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "gsea.R" | source_r }}

infile = {{in.infile | r}}
outdir = {{out.outdir | r}}
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
