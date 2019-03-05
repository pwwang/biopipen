{{rimport}}('sampleinfo.r')
library(superFreq)
options(stringsAsFactors = FALSE)

indir   = {{i.indir | quote}}
gfile   = {{i.gfile | quote}}
outdir  = {{o.outdir | quote}}
params  = {{args.params | R}}
nthread = {{args.nthread | R}}
baits   = {{args.baits | quote}}
ref     = {{args.ref | quote}}
genome  = {{args.genome | quote}}

# construct meta file
metafile = file.path(outdir, 'superFreq.meta.txt')
saminfo  = SampleInfo2(gfile, checkPaired = TRUE)
metadata = saminfo$as.superfreq.meta()
write.table(metadata, metafile, row.names = FALSE, sep = "\t", quote = FALSE)

# convert baits to bed
if (endsWith(baits, ".gff") || endsWith(baits, ".gtf")) {
	library(reticulate)
	baitfile = file.path(outdir, 'bait.bed')
	import('gff')$Gff$toBedfile(baits, baitfile, name = '{attributes[exon_id]}')
	baits = baitfile
}

# normal directory
normaldir = file.path(outdir, 'normalDirectory')
nbamdir   = file.path(normaldir, 'bam')
dir.create(nbamdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
for (normal in metadata[which(meta$NORMAL == "YES"), "BAM", drop = TRUE]) {
	file.symlink(normal, file.path(nbamdir, basename(normal)))
}

# resource dir
resourcedir = file.path(outdir, 'resource')

params$metaDataFile      = metafile
params$captureRegions    = baits
params$normalDirectory   = normaldir
params$Rdirectory        = file.path(outdir, 'R')
params$plotDirectory     = file.path(outdir, 'plots')
params$reference         = ref
params$genome            = genome
params$cpus              = nthread
params$resourceDirectory = resourcedir

data = do.call(superFreq, params)
printHTML(metaDataFile = metafile, outputFile = file.path(outdir, 'plots', 'superFreq.html'))
