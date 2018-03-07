library (Gviz)
tracks = NULL
region = unlist(strsplit({{in.region | quote}}, ':', fixed = T))
chrom  = region[1]
region = unlist(strsplit(region[2], '-', fixed = T))
start  = as.numeric(region[1])
end    = as.numeric(region[2])
if ({{ args.ideoTrack | R }} != F) {
	iparams = {{args.params.ideoTrack | Rlist}}
	iparams$chromosome = chrom
	iparams$genome     = {{args.genome | quote}}
	{% if args.ideoTrack | R | lambda x: x != 'TRUE' %}
	iparams$bands      = read.table({{args.ideoTrack | lambda x: 'r:gzfile("%s")' % x if x.endswith('.gz') else x | R}})
	colnames(iparams$bands) = c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')
	{% endif %}
	tracks = c(tracks, do.call(IdeogramTrack, iparams))
}
if ({{ args.axisTrack | R }}) {
	aparams = {{args.params.axisTrack | Rlist}}
    tracks  = c(tracks, do.call(GenomeAxisTrack, aparams))
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
	#if (dir.exists(t)) next  # allow empty track file list
	tdata = readRDS(t)
	if (tdata$trackType == 'InteractionTrack') {
		library(GenomicInteractions)
		irdata              = do.call(makeGenomicInteractionsFromFile, tdata$interactionParams)
		tdata$trackParams$x = irdata
		irtrackParams       = list(
			x          = '',
			chromosome = "",
			name       = NULL,
			start      = NULL,
  			end        = NULL
		)
		for (name in names(irtrackParams)) {
			irtrackParams    [[name]] = tdata$trackParams[[name]]
			tdata$trackParams[[name]] = NULL
		}
		track               = do.call(InteractionTrack, irtrackParams)
		displayPars(track)  = tdata$trackParams
	} else {
		track = do.call(tdata$trackType, tdata$trackParams)
	} 
    extracks = c (extracks, track)
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
	hstarts = sapply(hstarts, function(x) max(x, start))
	hends   = sapply(hends, function(x) min(x, end))
	hltrack = HighlightTrack(trackList = as.list(extracks), start = hstarts, end = hends, chromosome = chrom)
	tracks2plot = c(tracks2plot, hltrack)
} else {
	tracks2plot = c(tracks2plot, as.list(extracks))
}

devpars = {{args.devpars | Rlist}}
devpars$height = length(extracks) * devpars$height + 500 * as.integer(as.logical(length(tracks)))
do.call(png, c(list(filename = {{out.outfile | quote}}), devpars))
params = list(
	trackList = tracks2plot,
	from      = start,
	to        = end
)
params = c(params, {{args.params.general | Rlist}})
do.call(plotTracks, params)
dev.off()