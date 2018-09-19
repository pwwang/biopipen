{{rimport}}('__init__.r')
library (Gviz)

tracks = NULL
region = unlist(strsplit({{i.region | quote}}, ':', fixed = T))
chrom  = region[1]
region = unlist(strsplit(region[2], '-', fixed = T))
start  = as.numeric(region[1])
end    = as.numeric(region[2])

trackFiles = {{i.trkfiles | R}}
highlights = {{i.highlight | R}}
genome     = {{args.genome | R}}
ideoTrack  = {{args.ideoTrack | R}}
axisTrack  = {{args.axisTrack | R}}
geneTrack  = {{args.geneTrack | R}}
params     = {{args.params | R}}
scheme     = {{args.scheme | R}}
devpars    = {{args.devpars | R}}
outfile    = {{o.outfile | R}}

scheme.default = getScheme()
scheme         = update.list(scheme.default, scheme, recursive = T)

addScheme(scheme, "bioprocs")
options(Gviz.scheme = "bioprocs")

ntracks = 0
tracks  = list()

# ideoTrack
if (is.logical(ideoTrack) && ideoTrack) {
	ideofile = textConnection(readLines(gzcon(url('http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz'))))
} else if (!is.logical(ideoTrack)) {
	ideofile = ifelse(endsWith(ideoTrack, '.gz'), gzfile(ideoTrack), ideoTrack)
}
if (!is.logical(ideoTrack) || ideoTrack) {
	ideoParams = list(
		chromosome = chrom,
		genome     = genome,
		bands      = read.table(ideofile)
	)
	colnames(ideoParams$bands) = c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')
	tracks = c(tracks, do.call(IdeogramTrack, ideoParams))
}

# axisTrack
if (axisTrack) {
	tracks = c(tracks, GenomeAxisTrack())
}

# highlightable
hiableTracks = NULL
# geneTrack
if (!is.logical(geneTrack) || geneTrack) {
	ntracks = ntracks + 1
	data(geneModels)
	geneParams = list(
		range      = geneTrack,
		name       = 'Gene',
		chromosome = chrom
	)
	hiableTracks = c(hiableTracks, do.call(GeneRegionTrack, geneParams))
}

# other tracks
for (trackfile in trackFiles) {
	ntracks   = ntracks + 1
	trackData = readRDS(trackfile)

	# interaction track
	if (trackData$trackType == 'InteractionTrack') {
		library(GenomicInteractions)
		trackData$trackParams$x = do.call(makeGenomicInteractionsFromFile, trackData$interactionParams)
		irtrackParams = list()
		for (name in c('x', 'chromosome', 'name', 'start', 'end')) {
			irtrackParams[[name]] = trackData$trackParams[[name]]
			trackData$trackParams[[name]] = NULL
		}
		
		track = do.call(InteractionTrack, irtrackParams)
		displayPars(track)  = trackData$trackParams
	} else {
		track = do.call(trackData$trackType, trackData$trackParams)
	}
	hiableTracks = c(hiableTracks, track)
}

# highlight tracks
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
	hltrack = HighlightTrack(trackList = hiableTracks, start = hstarts, end = hends, chromosome = chrom)
	tracks  = c(tracks, hltrack)
} else {
	tracks  = c(tracks, hiableTracks)
}

devpars$height = (1 + ntracks) * devpars$height
do.call(png, c(list(filename = outfile), devpars))
params.default = list(
	trackList = tracks,
	from      = start,
	to        = end
)
params = c(params.default, params)
do.call(plotTracks, params)
dev.off()