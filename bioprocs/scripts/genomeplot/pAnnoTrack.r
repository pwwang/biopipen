#library (Gviz)
region = unlist(strsplit({{i.chrom | quote}}, ':', fixed = T))
chrom  = region[1]

infile = {{i.infile | quote}}
ext    = {{i.infile | ext | [1:] | quote}}
if (ext == 'bedx') {
	bname   = basename(infile)
	bedfile = file.path({{job.outdir | quote}}, substr(bname, 1, nchar(bname) - 1))
	file.symlink(infile, bedfile)
	infile = bedfile
}

params = list(
	range      = infile,
	genome     = {{args.genome | quote}},
	chromosome = chrom,
	name       = {{i.name | quote}}
)
params          = c(params, {{args.params | Rlist}})
ret             = list()
ret$trackType   = 'AnnotationTrack'
ret$trackParams = params

saveRDS (ret, {{o.outfile | quote}})
