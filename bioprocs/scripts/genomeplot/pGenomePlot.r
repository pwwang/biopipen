library (Gviz)
tracks = NULL
region = unlist(strsplit({{in.region | quote}}, ':', fixed = T))
chrom  = region[1]
region = unlist(strsplit(region[2], '-', fixed = T))
start  = as.numeric(region[1])
end    = as.numeric(region[2])
if ({{ args.showIdeo | R }}) {
    tracks = c (tracks, IdeogramTrack(genome={{args.genome | quote}}, chromosome=chrom))
}
if ({{ args.showAxis | R }}) {
    tracks = c (tracks, GenomeAxisTrack())
}
if ({{args.showGenes | R}}) {
	geneParams = list(
		genome     = {{args.genome | quote}},
		chromosome = chrom,
		from       = start,
		to         = end,
		name       = 'refGene',
		trackType  = 'GeneRegionTrack',
		table      = 'refGene',
		track      = 'refGene',
		rstarts    = 'exonStarts',
		rends      = 'exonEnds',
		gene       = "name",
		symbol     = "name2",
		transcript = "name",
		strand     = "strand"
	)
	geneParams = c(geneParams, {{args.params.geneTrack | Rlist}})
	geneTrack  = do.call(UcscTrack, geneParams)
	tracks     = c(tracks, geneTrack)
}
for (t in c({{in.trkfiles | acquote}})) {
	if (dir.exists(t)) next  # allow empty track file list
    tracks = c (tracks, readRDS(t))
}

do.call(png, c(list(file = {{out.outfile | quote}}), {{args.devpars | Rlist}}))
params = list(
	trackList = as.list(tracks),
	from      = start,
	to        = end
)
params = c(params, {{args.params.general | Rlist}})
do.call(plotTracks, params)
dev.off()