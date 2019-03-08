{{rimport}}('sampleinfo.r')
library(superFreq)
options(stringsAsFactors = FALSE)

indir    = {{i.indir | quote}}
gfile    = {{i.gfile | quote}}
outdir   = {{o.outdir | quote}}
params   = {{args.params | R}}
nthread  = {{args.nthread | R}}
baits    = {{args.baits | quote}}
resdir   = {{args.resdir | quote}}
ref      = {{args.ref | quote}}
genome   = {{args.genome | quote}}
cachedir = {{job.cachedir | quote}}
if (!file.exists(cachedir)) {
	dir.create(cachedir, showWarnings = FALSE, recursive = FALSE, mode = "0777")
}

# construct meta file
metafile = file.path(outdir, 'superFreq.meta.txt')
saminfo  = SampleInfo2$new(gfile, checkPaired = TRUE)
metadata = saminfo$as.superfreq.meta(datadir = indir)
write.table(metadata, metafile, row.names = FALSE, sep = "\t", quote = FALSE)

# try to find index of all bam files
for (bam in unique(unlist(metadata$BAM))) {
	reticulate::import('bioprocs')$utils$reference$bamIndex(bam, samtools = NULL, nthread = nthread)
}

# convert baits to bed
if (baits == "") {
	# try to find it in resource directory
	if (file.exists(file.path(resdir, genome))) {
		allbeds   = Sys.glob(file.path(resdir, genome, '*.bed'))
		namedbeds = Sys.glob(file.path(resdir, genome, '*.named*.bed'))
		baits     = setdiff(allbeds, namedbeds)
	}
} else if (endsWith(baits, ".gff") || endsWith(baits, ".gtf")) {
	baitfile = file.path(outdir, 'bait.bed')
	reticulate::import('gff')$Gff$toBedfile(baits, baitfile, name = '{attributes[gene_id]}')
	baits = baitfile
} 

# normal directory
normaldir_cache = file.path(cachedir, 'normalDirectory')
if (!file.exists(normaldir_cache)) {
	dir.create(normaldir_cache, showWarnings = FALSE, recursive = FALSE, mode = "0777")
}
normaldir = file.path(outdir, 'normalDirectory')
file.symlink(normaldir_cache, normaldir)
nbamdir   = file.path(normaldir, 'bam')
if (!file.exists(nbamdir)) {
	dir.create(nbamdir, showWarnings = FALSE, recursive = FALSE, mode = "0777")
}
for (normal in metadata[which(metadata$NORMAL == "YES"), "BAM", drop = TRUE]) {
	file.symlink(normal, file.path(nbamdir, basename(normal)))
	file.symlink(paste0(normal, '.bai'), file.path(nbamdir, basename(paste0(normal, '.bai'))))
}

Rdir_cache = file.path(cachedir, 'R')
Rdir = file.path(outdir, 'R')
if (!file.exists(Rdir_cache)) {
	dir.create(Rdir_cache, showWarnings = FALSE, recursive = FALSE, mode = "0777")
}
file.symlink(Rdir_cache, Rdir)

plotsdir_cache = file.path(cachedir, 'plots')
plotsdir = file.path(outdir, 'plots')
if (!file.exists(plotsdir_cache)) {
	dir.create(plotsdir_cache, showWarnings = FALSE, recursive = FALSE, mode = "0777")
}
file.symlink(plotsdir_cache, plotsdir)

# resource dir
resourcedir = file.path(outdir, 'resource')
file.symlink(resdir, resourcedir)

params$metaDataFile      = metafile
params$captureRegions    = baits
params$normalDirectory   = normaldir
params$Rdirectory        = Rdir
params$plotDirectory     = plotsdir
params$reference         = ref
params$genome            = genome
params$cpus              = nthread
params$resourceDirectory = resourcedir
tryCatch({
	data = do.call(superFreq, params)
}, error = function(e) {
	cat("Error encountered:\n", file = stderr())
	cat(paste0(e, "\n"), file = stderr())
})
printHTML(metaDataFile = metafile, outputFile = file.path(outdir, 'plots', 'superFreq.html'))
