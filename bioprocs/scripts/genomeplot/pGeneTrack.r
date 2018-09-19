#library (Gviz)
region = unlist(strsplit({{i.region | quote}}, ':', fixed = T))
chrom  = region[1]
region = unlist(strsplit(region[2], '-', fixed = T))
start  = as.numeric(region[1])
end    = as.numeric(region[2])
geneParams = list(
	genome     = {{args.genome | quote}},
	chromosome = chrom,
	from       = start,
	to         = end,
	name       = {{i.name | quote}},
	trackType  = 'GeneRegionTrack',
	table      = 'knownGene',
	track      = 'UCSC Genes',
	rstarts    = 'exonStarts',
	rends      = 'exonEnds',
	gene       = "name",
	symbol     = "name2",
	transcript = "name",
	strand     = "strand"
)
geneParams = c(geneParams, {{args.params | Rlist}})

ret             = list()
ret$trackType   = 'UcscTrack'
ret$trackParams = geneParams

saveRDS (ret, {{o.outfile | quote}})
