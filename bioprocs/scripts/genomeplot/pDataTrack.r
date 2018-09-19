#library (Gviz)
region = unlist(strsplit({{i.chrom | quote}}, ':', fixed = T))
chrom  = region[1]
params = list(
	range      = {{i.infile | quote}},
	genome     = {{args.genome | quote}},
	chromosome = chrom,
	name       = {{i.name | quote}}
)
params          = c(params, {{args.params | Rlist}})

ret             = list()
ret$trackType   = 'DataTrack'
ret$trackParams = params

saveRDS (ret, {{o.outfile | quote}})
