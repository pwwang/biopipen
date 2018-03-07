#library (Gviz)
region = unlist(strsplit({{in.region | quote}}, ':', fixed = T))
chrom  = region[1]
region = unlist(strsplit(region[2], '-', fixed = T))
start  = as.numeric(region[1])
end    = as.numeric(region[2])
params = list(
	genome     = {{args.genome | quote}},
	chromosome = chrom,
	from       = start,
	to         = end,
	name       = {{in.name | quote}},
	trackType  = {{in.trackType | quote}},
	track      = {{in.track | quote}}
)
params = c(params, {{args.params | Rlist}})

ret = list()
ret$trackType = 'UcscTrack'
ret$trackParams = params

saveRDS (ret, {{out.outfile | quote}})
