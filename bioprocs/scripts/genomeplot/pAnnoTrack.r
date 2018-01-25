library (Gviz)
region = unlist(strsplit({{in.chrom | quote}}, ':', fixed = T))
chrom  = region[1]

infile = {{in.infile | quote}}
ext    = {{in.infile | ext | [1:] | quote}}
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
	name       = {{in.name | quote}}
)
params = c(params, {{args.params | Rlist}})
track  = do.call(AnnotationTrack, params)
saveRDS (track, {{out.outfile | quote}})
