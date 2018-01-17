require('maftools')

maffiles = c(Sys.glob(file.path({{in.indir | quote}}, "*.maf")), Sys.glob(file.path({{in.indir | quote}}, "*.maf.gz")))
if (length(maffiles) == 0) {
	stop('No maf files found in input directory!')
}
maffile  = maffiles[1]
if (length(maffiles) > 1) {
	cat(paste0('pyppl.log.warning: You multiple MAF files in input directory, using the first one: ', maffile, '\n'), file = stderr())
}

annofiles = c(Sys.glob(file.path({{in.indir | quote}}, "*annot.tsv")), Sys.glob(file.path({{in.indir | quote}}, "*annotion.txt")))
if (length(annofiles) == 0) {
	annofile = NULL
} else {
	annofile  = annofiles[1]
	if (length(annofiles) > 1) {
		cat(paste0('pyppl.log.warning: You multiple annotation files in input directory, using the first one: ', annofile, '\n'), file = stderr())
	}
}

cntables = c(Sys.glob(file.path({{in.indir | quote}}, "*cnv.tsv")), Sys.glob(file.path({{in.indir | quote}}, "*cnv.txt")))
if (length(cntables) == 0) {
	cntable = NULL
} else {
	cntable  = cntables[1]
	if (length(cntables) > 1) {
		cat(paste0('pyppl.log.warning: You multiple cnTable files in input directory, using the first one: ', cntable, '\n'), file = stderr())
	}
}

alllesions = Sys.glob(file.path({{in.indir | quote}}, 'all_lesions.conf_*.txt'))
if (length(alllesions) == 0) {
	alllesion = NULL
} else {
	alllesion  = alllesions[1]
	if (length(alllesions) > 1) {
		cat(paste0('pyppl.log.warning: You multiple all-lesion files in input directory, using the first one: ', alllesion, '\n'), file = stderr())
	}
}

ampgenes = Sys.glob(file.path({{in.indir | quote}}, 'amp_genes.conf_*.txt'))
if (length(ampgenes) == 0) {
	ampgene = NULL
} else {
	ampgene  = ampgenes[1]
	if (length(ampgenes) > 1) {
		cat(paste0('pyppl.log.warning: You multiple amp-gene files in input directory, using the first one: ', ampgene, '\n'), file = stderr())
	}
}

delgenes = Sys.glob(file.path({{in.indir | quote}}, 'del_genes.conf_*.txt'))
if (length(delgenes) == 0) {
	delgene = NULL
} else {
	delgene  = delgenes[1]
	if (length(delgenes) > 1) {
		cat(paste0('pyppl.log.warning: You multiple del-gene files in input directory, using the first one: ', delgene, '\n'), file = stderr())
	}
}

gisscores = Sys.glob(file.path({{in.indir | quote}}, 'scores.gistic'))
if (length(gisscores) == 0) {
	gisscore = NULL
} else {
	gisscore  = gisscores[1]
}

hasGistic = FALSE
if (!is.null(alllesion) && !is.null(ampgene) && !is.null(delgene) && !is.null(gisscore)) {
	hasGistic = TRUE
	lamlGistic = readGistic(gisticAllLesionsFile = alllesion, gisticAmpGenesFile = ampgene, gisticDelGenesFile = delgene, gisticScoresFile = gisscore, isTCGA = {{args.isTCGA | R}})
}

laml    = read.maf(
	maf                  = maffile,
	clinicalData         = annofile,
	gisticAllLesionsFile = alllesion,
	gisticAmpGenesFile   = ampgene,
	gisticDelGenesFile   = delgene,
	gisticScoresFile     = gisscore,
	cnTable              = cntable,
	isTCGA               = {{args.isTCGA | R}},
	vc_nonSyn            = {% if args.mutypes %}{{args.mutypes | Rvec}}{% else %}NULL{% endif %}
)

genesum = getGeneSummary(laml)
samsum  = getSampleSummary(laml)
genes   = genesum$Hugo_Symbol[1:{{args.ngenes}}]
nsample = nrow(samsum)
samples = samsum$Tumor_Sample_Barcode

devpars = {{args.devpars | Rlist}}
devpars2 = devpars
devpars2$width = devpars2$width * 2
#### summary
{% if args.plots.summary %}
cat('## Plotting summary ...\n', file = stderr())
summaryplot = file.path({{out.outdir | quote}}, 'summary.png')
do.call(png, c(list(filename = summaryplot), devpars2))
tryCatch({
	do.call(plotmafSummary, c(list(maf = laml), {{args.params.summary | Rlist}}))
}, error = function(e) {
	cat('pyppl.log.warning: Failed to plot summary:', e, '\n', file = stderr())	
})
dev.off()
{% endif %}

#### oncoplot
{% if args.plots.oncoplot %}
cat('## Plotting oncoplot ...\n', file = stderr())
oncoplotfile = file.path({{out.outdir | quote}}, 'oncoplot.png')
do.call(png, c(list(filename = oncoplotfile), devpars2))
params = c(list(maf = laml), {{args.params.oncoplot | Rlist}})
if (!'top' %in% names(params) && !'genes' %in% names(params)){
	params$top = {{args.ngenes}}
}

if (!is.null(annofile)) {
	params$sortByAnnotation = TRUE
}
tryCatch({
	do.call(oncoplot, params)
}, error = function(e) {
	cat('pyppl.log.warning: Failed to plot oncoplot:', e, '\n', file = stderr())	
})
dev.off()
{% endif %}

#### oncostrip
{% if args.plots.oncostrip %}
cat('## Plotting oncostrip ...\n', file = stderr())
oncostripfile = file.path({{out.outdir | quote}}, 'oncostrip.png')
params        = {{args.params.oncostrip | Rlist}}
params$maf    = laml
if (!'genes' %in% names(params)) {
	params$genes = genes
}
do.call(png, c(list(filename = oncostripfile), devpars2))
tryCatch({
	do.call(oncostrip, params)
}, error = function(e) {
	cat('pyppl.log.warning: Failed to plot oncostrip:', e, '\n', file = stderr())	
})
dev.off()
{% endif %}

#### titv
{% if args.plots.titv %}
cat('## Plotting titv ...\n', file = stderr())
titvplot = file.path({{out.outdir | quote}}, 'titv.png')
titvobj  = titv(maf = laml, plot = FALSE, useSyn = TRUE)
params   = {{args.params.titv | Rlist}}
params$res = titvobj
do.call(png, c(list(filename = titvplot), devpars))
tryCatch({
	do.call(plotTiTv, params)
}, error = function(e) {
	cat('pyppl.log.warning: Failed to plot titv:', e, '\n', file = stderr())	
})
dev.off()
{% endif %}

#### lollipop
{% if args.plots.lollipop %}
cat('## Plotting lollipops ...\n', file = stderr())

params     = {{args.params.lollipop | Rlist}}
params$maf = laml
if ('genes' %in% names(params)) {
	if ('top' %in% names(params)) params$top = NULL
	lpgenes = params$genes
	params$genes = NULL
} else if ('top' %in% names(params)) {
	lpgenes = genes[1:params$top]
	params$top = NULL
} else {
	lpgenes = genes
}

lollipopPlotSingle = function(gene) {
	ps = params
	ps$gene = gene
	lpdir  = file.path({{out.outdir | quote}}, 'lollipops')
	dir.create(lpdir, showWarnings = F)
	lpplot = file.path(lpdir, paste0(ps$gene, '.lollipop.png'))
	do.call(png, c(list(filename = lpplot), devpars2))
	tryCatch({
		do.call(lollipopPlot, ps)
	}, error = function(e) {
		cat(paste('pyppl.log.warning: Failed to plot lollipop for', gene ,':', e, '\n'), file = stderr())	
	})
	dev.off()
}
{% 	if args.nthread | lambda x: x == 1 %}
for (gene in lpgenes) {
	lollipopPlotSingle(gene)
}
{% 	else %}
library(doParallel)
cl = makeCluster({{args.nthread}})
registerDoParallel(cl)
foreach(i=1:length(lpgenes), .verbose = T, .packages=c('maftools')) %dopar% {
	lollipopPlotSingle(lpgenes[i])
}
stopCluster(cl)
{% 	endif %}
{% endif %}

#### cbsseg
{% if args.plots.cbsseg %}
cat('## Plotting cbsseg ... \n', file = stderr())
params  = c(list(maf = laml), {{args.params.cbsseg | Rlist}})
sams    = if (!'tsb' %in% names(params)) samples else params$tsb

cbssegSingle = function(sam) {
	segfiles = c(
		Sys.glob(file.path({{in.indir | quote}}, paste0(sam, '*.seg.txt'))), 
		Sys.glob(file.path({{in.indir | quote}}, paste0(gsub('-', '.', sam), '*.seg.txt'))),
		Sys.glob(file.path({{in.indir | quote}}, paste0(gsub('_', '.', sam), '*.seg.txt')))
	)
	if (length(segfiles) == 0) {
		cat(paste('No seg file found for sample:', sam ,', skip.\n'), file = stderr())
	} else {
		segdir  = file.path({{out.outdir | quote}}, 'cbssegs')
		dir.create(segdir, showWarnings = F)
		segplot = file.path(segdir, paste0(sam, '.cbsseg.png'))
		do.call(png, c(list(filename = segplot), devpars2))
		tryCatch({
			do.call(plotCBSsegments, c(list(cbsFile = segfiles[1]), params))
		}, error = function(e) {
			cat(paste('pyppl.log.warning: Failed to plot cbsseg:', e, '\n'), file = stderr())	
		})
		dev.off()
	}
}
{% 	if args.nthread | lambda x: x == 1 %}
for (sam in sams) {
	cbssegSingle(sam)
}
{% 	else %}
library(doParallel)
cl = makeCluster({{args.nthread}})
registerDoParallel(cl)
foreach(i=1:length(sams), .verbose = T, .packages=c('maftools')) %dopar% {
	cbssegSingle(sams[i])
}
stopCluster(cl)
{% 	endif %}
{% endif %}

#### rainfall
{% if args.plots.rainfall %}
cat('## Plotting rainfall ... \n', file = stderr())
params = c(list(maf = laml), {{args.params.rainfall | Rlist}})
if ('tsb' %in% names(params)) {
	sams = params$tsb
	prams$tsb = NULL
} else {
	sams = samples
}

rainfallSingle = function(sam) {
	rfdir  = file.path({{out.outdir | quote}}, 'rainfalls')
	dir.create(rfdir, showWarnings = F)
	rfplot = file.path(rfdir, paste0(sam, '.rainfall.png'))
	do.call(png, c(list(filename = rfplot), devpars2))
	tryCatch({
		do.call(rainfallPlot, c(list(tsb = sam), params))
	}, error = function(e){
		if (!params$detectChangePoints) {
			cat(paste('pyppl.log.warning: Failed to plot rainfall without detecting change points for sample:', sam, '\n'), file = stderr())
		} else {
			tryCatch({
				do.call(rainfallPlot, c(list(tsb = sam, detectChangePoints = F), params))
			}, error = function(e){
				cat(paste('pyppl.log.warning: Failed to plot rainfall without detecting change points for sample', sam, ':', e, '\n'), file = stderr())
			})
		}
	})
	dev.off()
}
{% 	if args.nthread | lambda x: x == 1 %}
for (sam in sams) {
	rainfallSingle(sam)
}
{% 	else %}
library(doParallel)
cl = makeCluster({{args.nthread}})
registerDoParallel(cl)
foreach(i=1:length(sams), .verbose = T, .packages=c('maftools')) %dopar% {
	rainfallSingle(sams[i])
}
stopCluster(cl)
{% 	endif %}
{% endif %}

#### tcgacomp
{% if args.plots.tcgacomp %}
cat('## Plotting tcgacompare ... \n', file = stderr())
tcgacompplot = file.path({{out.outdir | quote}}, 'tcgacomp.png')
do.call(png, c(list(filename = tcgacompplot), devpars2))
params = c(list(maf = laml), {{args.params.tcgacomp | Rlist}})
if (!'cohortName' %in% names(params)) {
	params$cohortName = unlist(strsplit(basename(maffile), '.', fixed = T))[1]
}
tryCatch({
	do.call(tcgaCompare, params)
}, error = function(e){
	cat(paste('pyppl.log.warning: Failed to plot tcgacomp:', e, '\n'), file = stderr())
})
dev.off()
{% endif %}

#### vaf
{% if args.plots.vaf %}
cat('## Plotting vaf ... \n', file = stderr())
params = c(list(maf = laml), {{args.params.vaf | Rlist}})
if (!'vafCol' %in% names(params)) {
	cat('pyppl.log.warning: No vafCol(args.params.vaf.vafCol) provided, skip plotting VAF.\n', file = stderr())
} else {
	vafplot = file.path({{out.outdir | quote}}, 'vaf.png')
	do.call(png, c(list(filename = vafplot), devpars))
	tryCatch({
		do.call(plotVaf, params)
	}, error = function(e){
		cat(paste('pyppl.log.warning: Failed to plot vaf:', e, '\n'), file = stderr())
	})
	dev.off()
}
{% endif %}

#### genecloud
{% if args.plots.genecloud %}
cat('## Plotting genecloud ... \n', file = stderr())
gcplot = file.path({{out.outdir | quote}}, 'genecloud.png')
do.call(png, c(list(filename = gcplot), devpars))
params = c(list(input = laml), {{args.params.genecloud | Rlist}})
tryCatch({
	do.call(geneCloud, params)
}, error = function(e){
	cat(paste('pyppl.log.warning: Failed to plot genecloud:', e, '\n'), file = stderr())
})
dev.off()
{% endif %}


#### gisticGenome
{% if args.plots.gisticGenome %}
cat('## Plotting gisticGenome ... ', file = stderr())
if (!hasGistic) {
	cat('No gistic files found, skip\n', file = stderr())
} else {
	gisgplot = file.path({{out.outdir | quote}}, 'gisticGenome.png')
	do.call(png, c(list(filename = gisgplot), devpars2))
	params = c(list(gistic = lamlGistic), {{args.params.gisticGenome | Rlist}})
	tryCatch({
		do.call(gisticChromPlot, params)
	}, error = function(e){
		cat(paste('pyppl.log.warning: Failed to plot gisticGenome:', e, '\n'), file = stderr())
	})
	dev.off()
}
{% endif %}

#### gisticBubble
{% if args.plots.gisticBubble %}
cat('## Plotting gisticBubble ...\n', file = stderr())
if (!hasGistic) {
	cat('No gistic files found, skip\n', file = stderr())
} else {
	gisbplot = file.path({{out.outdir | quote}}, 'gisticBubble.png')
	do.call(png, c(list(filename = gisbplot), devpars))
	params = c(list(gistic = lamlGistic), {{args.params.gisticBubble | Rlist}})
	tryCatch({
		do.call(gisticBubblePlot, params)
	}, error = function(e){
		cat(paste('pyppl.log.warning: Failed to plot gisticBubble:', e, '\n'), file = stderr())
	})
	dev.off()
}
{% endif %}

#### gisticOncoplot
{% if args.plots.gisticOncoplot %}
cat('## Plotting gisticOncoplot ...\n', file = stderr())
if (!hasGistic) {
	cat('No gistic files found, skip\n', file = stderr())
} else {
	gisoncoplotfile = file.path({{out.outdir | quote}}, 'gisticOncoplot.png')
	params = c(list(gistic = lamlGistic), {{args.params.gisticOncoplot | Rlist}})
	if (!is.null(annofile)) {
		params$clinicalData     = getClinicalData(x = laml)
		params$sortByAnnotation = TRUE
		if (!'top' %in% names(params)) 
			params$top = {{args.ngenes}}
	}
	do.call(png, c(list(filename = gisoncoplotfile), devpars2))
	tryCatch({
		do.call(gisticOncoPlot, params)
	}, error = function(e) {
		if (is.null(annofile)) {
			cat('pyppl.log.warning: Cannot generate gisticOncoplot without clinic features, skip.\n', file = stderr())
		} else {
			tryCatch({
				params$clinicalFeatures = NULL
				do.call(gisticOncoPlot, params)
			}, error = function(e) {
				cat('pyppl.log.warning: Cannot generate gisticOncoplot even without clinic features, skip.\n', file = stderr())
			})
		}
	})
	dev.off()
}
{% endif %}

#### somInteraction
{% if args.plots.somInteraction %}
cat('## Plotting somInteraction ...\n', file = stderr())
somInteractionfile = file.path({{out.outdir | quote}}, 'somInteraction.png')
params = c(list(maf = laml), {{args.params.somInteraction | Rlist}})
if (!'top' %in% names(params)) 
	params$top = {{args.ngenes}}
do.call(png, c(list(filename = somInteractionfile), devpars))
tryCatch(
	{
		do.call(somaticInteractions, params)
	}, error = function(e) {
		cat(paste('-  Unable to plot somaticInteractions, with error:', e,'\n'), file = stderr())
	}
)
dev.off()
{% endif %}

#### oncodrive
{% if args.plots.oncodrive %}
cat('## Plotting oncodrive ...\n', file = stderr())
oncodrivefile = file.path({{out.outdir | quote}}, 'oncodrive.png')
OcParams = list(maf = laml, AACol = NULL, minMut = 5, pvalMethod = "zscore",
  				nBgGenes = 100, bgEstimate = TRUE, ignoreGenes = NULL)		
params = c({{args.params.oncodrive | Rlist}})
for (name in names(params)) {
	if (name %in% names(OcParams))
		OcParams[[name]] = params[[name]]
}
tryCatch(
	{
		sig = do.call(oncodrive, OcParams)

		PocParams = list(res = sig, fdrCutOff = 0.05, useFraction = FALSE,
						colCode = NULL, labelSize = 2)
		for (name in names(params)) {
			if (name %in% names(PocParams))
				PocParams[[name]] = params[[name]]
		}
		do.call(png, c(list(filename = oncodrivefile), devpars))
		do.call(plotOncodrive, PocParams)
		dev.off()
	}, error = function(e) {
		cat(paste('-  Unable to plot oncodrive, with error:', e,'\n'), file = stderr())	
	}
)
{% endif %}

#### pfam
{% if args.plots.pfam %}
cat('## Plotting pfam ...\n', file = stderr())
pfamfile = file.path({{out.outdir | quote}}, 'pfam.png')
params = c(list(maf = laml), {{args.params.pfam | Rlist}})
if (!'top' %in% names(params)) 
	params$top = {{args.ngenes}}
do.call(png, c(list(filename = pfamfile), devpars))
tryCatch(
	{
		do.call(pfamDomains, params)
	}, error = function(e) {
		cat(paste('-  Unable to plot pfam, with error:', e,'\n'), file = stderr())	
	}
)
dev.off()
{% endif %}

#### pancan
{% if args.plots.pancan %}
siggenes = c(Sys.glob(file.path({{in.indir | quote}}, '*sig_genes.txt.gz')), Sys.glob(file.path({{in.indir | quote}}, '*sig_genes.txt')))
if (length(siggenes) > 0) {
	siggene = siggenes[1]
	if (length(siggenes) > 1) {
		cat(paste0('pyppl.log.warning: You multiple mutsig files in input directory, using the first one: ', siggene, '\n'), file = stderr())
	}
	params  = {{args.params.pancan | Rlist}}
	params$mutsigResults = siggene
	if (!'cohortName' %in% names(params)) {
		params$cohortName = unlist(strsplit(basename(maffile), '.', fixed = T))[1]
	}
	if (!'inputSampleSize ' %in% names(params)) {
		params$inputSampleSize  = nsample
	}
	pancanfile = file.path({{out.outdir | quote}}, 'pancan.png')
	do.call(png, c(list(filename = pancanfile), devpars))
	tryCatch({
		do.call(pancanComparison, params)
	}, error = function(e){
		cat(paste('pyppl.log.warning: Failed to plot pancan:', e, '\n'), file = stderr())
	})
	dev.off()
}
{% endif %}

#### survival
{% if args.plots.survival %}
cat('## Plotting survivals ...\n', file = stderr())

params        = {{args.params.survival | Rlist}}
params$maf    = laml
params$isTCGA = {{args.isTCGA | R}}
if (!'time' %in% names(params) || !'Status' %in% names(params)) {
	cat('-  No time/Status column specified for survival analysis!', file = stderr())
} else {
	survivalSingle = function(gene) {
		ps = params
		ps$genes = gene
		svdir  = file.path({{out.outdir | quote}}, 'survivals')
		dir.create(svdir, showWarnings = F)
		svplot = file.path(svdir, paste0(ps$genes, '.survival.png'))
		do.call(png, c(list(filename = svplot), devpars))
		tryCatch({
			do.call(mafSurvival, ps)
		}, error = function(e){
			cat(paste('pyppl.log.warning: Failed to plot survival for gene', gene ,':', e, '\n'), file = stderr())
		})
		dev.off()
	}
	if ('genes' %in% names(params)) {
		svgenes = params$genes
		params$genes = NULL
	} else {
		svgenes = genes
	}
	{% 	if args.nthread | lambda x: x == 1 %}
	for (gene in svgenes) {
		survivalSingle(gene)
	}
	{% 	else %}
	library(doParallel)
	cl = makeCluster({{args.nthread}})
	registerDoParallel(cl)
	foreach(i=1:length(svgenes), .verbose = T, .packages=c('maftools')) %dopar% {
		survivalSingle(svgenes[i])
	}
	stopCluster(cl)
	{% 	endif %}
}
{% endif %}

#### heterogeneity
{% if args.plots.heterogeneity %}
cat('## Plotting heterogeneity ...\n', file = stderr())
if ('tsb' %in% names(params)) {
	sams = params$tsb
	params$tsb = NULL
 } else {
	sams = samples
 }

ihParams = list(maf = laml, top = 5, vafCol = NULL,
  ignChr = NULL, minVaf = 0, maxVaf = 1,
  dirichlet = FALSE)

params = c({{args.params.heterogeneity | Rlist}})
for (name in names(params)) {
	if (name %in% names(ihParams))
		ihParams[[name]] = params[[name]]
}

PocParams = list(genes = NULL, showCNvars = FALSE,
  savePlot = FALSE, width = 6, height = 5, colors = NULL)
for (name in names(params)) {
	if (name %in% names(PocParams))
		PocParams[[name]] = params[[name]]
}

if (is.null(ihParams$vafCol)) {
	cat('pyppl.log.warning: No vafCol(args.params.heterogeneity.vafCol) provided, skip plotting heterogeneity.\n', file = stderr())
} else {
	heteroSingle = function(sam) {
		segfiles = c(
			Sys.glob(file.path({{in.indir | quote}}, paste0(sam, '*.seg.txt'))), 
			Sys.glob(file.path({{in.indir | quote}}, paste0(gsub('-', '.', sam), '*.seg.txt'))),
			Sys.glob(file.path({{in.indir | quote}}, paste0(gsub('_', '.', sam), '*.seg.txt')))
		)
		tryCatch({
			if (length(segfiles) == 0) {
				hetero = do.call(inferHeterogeneity, c(list(tsb = sam), ihParams))
			} else {
				hetero = do.call(inferHeterogeneity, c(list(tsb = sam, segFile = segfiles[1]), ihParams))
			}
			hgdir  = file.path({{out.outdir | quote}}, 'heterogeneities')
			dir.create(hgdir, showWarnings = F)
			heterogeneityfile = file.path(hgdir, paste0(sam, '.heterogeneity.png'))
			do.call(png, c(list(filename = heterogeneityfile), devpars))
			do.call(plotClusters, c(list(clusters = hetero, tsb = sam), PocParams))
			dev.off()
		}, error = function(e) {
			cat('pyppl.log.warning: Failed to plot heterogeneity for sample:', sam, ', skip.\n', file = stderr())
		})
	}
	{% 	if args.nthread | lambda x: x == 1 %}
	for (sam in sams) {
		heteroSingle(sam)
	}
	{% 	else %}
	library(doParallel)
	cl = makeCluster({{args.nthread}})
	registerDoParallel(cl)
	foreach(i=1:length(sams), .verbose = T, .packages=c('maftools')) %dopar% {
		heteroSingle(sams[i])
	}
	stopCluster(cl)
	{% 	endif %}
}
{% endif %}

#### signature
{% if args.plots.signature %}
cat('## Plotting signature ...\n', file = stderr())
require('NMF')
params = {{args.params.signature | Rlist}}
tmParams = list(prefix = NULL, add = TRUE, ignoreChr = NULL, useSyn = TRUE, fn = NULL)
for (name in names(params)) {
	if (name %in% names(tmParams)) {
		tmParams[[name]] = params[[name]]
		params  [[name]] = NULL
	}
}
tmParams$ref_genome = {{args.ref | quote}}
tmParams$maf = laml
tm = do.call(trinucleotideMatrix, tmParams)
# apobec
apobecfile = file.path({{out.outdir | quote}}, 'apobec.png')
do.call(png, c(list(filename = apobecfile), devpars))
tryCatch(
	{
		plotApobecDiff(tnm = tm, maf = laml)
	}, error = function(e) {
		cat(paste('-  Unable to plot apobec signature:', e), file = stderr())
	}
)
dev.off()

# sigs
sgParams = list(n = NULL, nTry = 6, plotBestFitRes = FALSE, parallel = NULL)
for (name in names(params)) {
	if (name %in% names(sgParams)) {
		sgParams[[name]] = params[[name]]
		params  [[name]] = NULL
	}
}
sgParams$mat = tm
tryCatch({
	sigs  = do.call(extractSignatures, sgParams)
	sigfile = file.path({{out.outdir | quote}}, 'signature.png')
	plotSignatures(sigs)
	dev.off()
	# cosine similarity against validated signatures
	heatmap(sigs$coSineSimMat, file.path({{out.outdir | quote}}, 'signature-sim.png'), devpars = devpars, dendro = 'col')
	# signature contrib
	sigcfile = file.path({{out.outdir | quote}}, 'signature-contrib.png')
	plotSignatures(nmfRes = sigs, contributions = TRUE)
}, error = function(e) {
	cat('pyppl.log.warning: Failed to extract signatures', ', skip.\n', file = stderr())
	
})


{% endif %}





