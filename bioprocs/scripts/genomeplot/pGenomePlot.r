library (Gviz)
tracks = NULL
region = unlist(strsplit({{in.region | quote}}, ':', fixed = T))
chrom  = region[1]
region = unlist(strsplit(region[2], '-', fixed = T))
start  = as.numeric(region[1])
end    = as.numeric(region[2])
if ({{ args.ideoTrack | R }} != F) {
{% if args.ideoTrack | R | lambda x: x == 'TRUE' %}
	tracks = c (tracks, IdeogramTrack(genome={{args.genome | quote}}, chromosome=chrom))
{% else %}
	tracks = c (tracks, IdeogramTrack(genome={{args.genome | quote}}, chromosome=chrom, bands = {{args.ideoTrack | R}}))
{% endif %}
}
if ({{ args.axisTrack | R }}) {
    tracks = c (tracks, GenomeAxisTrack())
}


extracks = NULL
if ({{args.geneTrack | R}} != F) {
	# should be a gtf file
	data(geneModels)
	geneParams = list(
		range = {{args.geneTrack | quote}},
		name  = 'Gene',
		chromosome = chrom
	)
	geneParams = c(geneParams, {{args.params.geneTrack | Rlist}})
	geneTrack  = do.call(GeneRegionTrack, geneParams)
	extracks   = c(extracks, geneTrack)
}
for (t in c({{in.trkfiles | acquote}})) {
	if (dir.exists(t)) next  # allow empty track file list
    extracks = c (extracks, readRDS(t))
}

tracks2plot = as.list(tracks)
highlights  = {{in.highlight | R}}
if (nchar(highlights) > 0) {
	hregs   = unlist(strsplit(highlights, ';', fixed = T))
	hstarts = NULL
	hends   = NULL
	for (hreg in hregs) {
		reg = unlist(strsplit(hreg, '-', fixed = T))
		hstarts = c(hstarts, as.numeric(reg[1]))
		hends   = c(hends  , as.numeric(reg[2]))
	}
	hltrack = HighlightTrack(trackList = as.list(extracks), start = hstarts, end = hends, chromosome = chrom)
	tracks2plot = c(tracks2plot, hltrack)
} else {
	tracks2plot = c(tracks2plot, as.list(extracks))
}

do.call(png, c(list(filename = {{out.outfile | quote}}), {{args.devpars | Rlist}}))
params = list(
	trackList = tracks2plot,
	from      = start,
	to        = end
)
params = c(params, {{args.params.general | Rlist}})
do.call(plotTracks, params)
dev.off()