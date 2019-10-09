{{rimport}}('__init__.r', 'plot.r')
library(maftools)

indir    = {{i.indir | quote}}
msdir    = {{i.msdir | quote}}
outdir   = {{o.outdir | quote}}
isTCGA   = {{args.isTCGA | R}}
extypes  = {{args.extypes | R}}
devpars  = {{args.devpars | R}}
plots    = {{args.plot | R}}
mtparams = {{args.params | R}}
ngenes   = {{args.ngenes | R}}
nthread  = {{args.nthread | R}}
genome   = {{args.genome | quote}}

setwd(outdir)
if (file.exists(indir) && dir.exists(indir)) {
	# you probably have multiple maf files
	maffiles = c(Sys.glob(file.path(indir, "*.maf")), Sys.glob(file.path(indir, "*.maf.gz")))

	if (length(maffiles) == 0) {
		stop('No maf files found in input directory!')
	}
	maffile  = maffiles[1]
	if (length(maffiles) > 1) {
		log2pyppl('You multiple MAF files in input directory, using the first one:', maffile, level = 'warning')
	}
} else if (file.exists(indir)) {
	maffile = indir
	indir = {{job.indir | quote}}
}

# get annotation file
annofiles = c(Sys.glob(file.path(indir, "*annot.tsv")), Sys.glob(file.path(indir, "*annotion.txt")))
if (length(annofiles) == 0) {
	annofile = NULL
} else {
	annofile  = annofiles[1]
	if (length(annofiles) > 1) {
		log2pyppl('You multiple annotation files in input directory, using the first one:', annofile, level = 'warning')
	}
}

# custom cnv data, see https://support.bioconductor.org/p/94113/
cntables = c(Sys.glob(file.path(indir, "*cnv.tsv")), Sys.glob(file.path(indir, "*cnv.txt")))
if (length(cntables) == 0) {
	cntable = NULL
} else {
	cntable  = cntables[1]
	if (length(cntables) > 1) {
		log2pyppl('You multiple cnTable files in input directory, using the first one:', cntable, level = 'warning')
	}
}

# Gistic cnv data
alllesions = Sys.glob(file.path(indir, 'all_lesions.conf_*.txt'))
if (length(alllesions) == 0) {
	alllesion = NULL
} else {
	alllesion  = alllesions[1]
	if (length(alllesions) > 1) {
		log2pyppl('You multiple all-lesion files in input directory, using the first one:', alllesion, level = 'warning')
	}
}

ampgenes = Sys.glob(file.path(indir, 'amp_genes.conf_*.txt'))
if (length(ampgenes) == 0) {
	ampgene = NULL
} else {
	ampgene  = ampgenes[1]
	if (length(ampgenes) > 1) {
		log2pyppl('You multiple amp-gene files in input directory, using the first one:', ampgene, level = 'warning')
	}
}

delgenes = Sys.glob(file.path(indir, 'del_genes.conf_*.txt'))
if (length(delgenes) == 0) {
	delgene = NULL
} else {
	delgene  = delgenes[1]
	if (length(delgenes) > 1) {
		log2pyppl('You multiple del-gene files in input directory, using the first one:', delgene, level = 'warning')
	}
}

gisscores = Sys.glob(file.path(indir, 'scores.gistic'))
if (length(gisscores) == 0) {
	gisscore = NULL
} else {
	gisscore  = gisscores[1]
}

# is gistic data complete?
hasGistic = FALSE
if (!is.null(alllesion) && !is.null(ampgene) && !is.null(delgene) && !is.null(gisscore)) {
	hasGistic = TRUE
	lamlGistic = readGistic(gisticAllLesionsFile = alllesion, gisticAmpGenesFile = ampgene, gisticDelGenesFile = delgene, gisticScoresFile = gisscore, isTCGA = isTCGA)
}

mafdata = read.table.inopts(maffile, list(cnames = TRUE, rnames = FALSE, fill = TRUE))
alltypes = levels(as.factor(mafdata$Variant_Classification))

laml = read.maf(
	maf                  = maffile,
	clinicalData         = annofile,
	gisticAllLesionsFile = alllesion,
	gisticAmpGenesFile   = ampgene,
	gisticDelGenesFile   = delgene,
	gisticScoresFile     = gisscore,
	cnTable              = cntable,
	isTCGA               = isTCGA,
	vc_nonSyn            = setdiff(alltypes, extypes)
)

genesum = getGeneSummary(laml)
samsum  = getSampleSummary(laml)
genes   = genesum$Hugo_Symbol[1:ngenes]
nsample = nrow(samsum)

if (is.true(msdir) && dir.exists(msdir)) {
	siggene = c(
		Sys.glob(file.path(msdir, '*sig_genes.txt.gz')), 
		Sys.glob(file.path(msdir, '*sig_genes.txt')))[1]
} else if (is.true(msdir) && file.exists(msdir)) {
	siggene = msdir
} else {
	siggene = NULL
}
if (!is.null(siggene)) {
	msdata = read.table.inopts(siggene, list(rnames = TRUE, cnames = TRUE))
	genes = rownames(msdata)[1:ngenes]
}

samples = samsum$Tumor_Sample_Barcode

devpars2       = devpars
devpars2$width = devpars2$width * 2
#### summary
if (plots$summary) {
	logger('## Plotting summary ...')
	summaryplot = file.path(outdir, 'summary.png')
	do.call(png, c(list(filename = summaryplot), devpars2))
	tryCatch({
		do.call(plotmafSummary, c(list(maf = laml), mtparams$summary))
	}, error = function(e) {
		log2pyppl('Failed to plot summary:', e, level = 'warning')
	})
	dev.off()
}

#### oncoplot
if (plots$oncoplot) {
	cat('## Plotting oncoplot ...\n')
	oncoplotfile = file.path(outdir, 'oncoplot.png')
	do.call(png, c(list(filename = oncoplotfile), devpars2))
	params = c(list(maf = laml), mtparams$oncoplot)
	if (!'top' %in% names(params) && !'genes' %in% names(params)){
		#params$top = ngenes
		params$genes = genes
	} else if ('top' %in% names(params)) {
		params$genes = genes[1:params$top]
	}

	if (!is.null(annofile)) {
		params$sortByAnnotation = TRUE
	}
	tryCatch({
		do.call(oncoplot, params)
	}, error = function(e) {
		log2pyppl('Failed to plot oncoplot:', e, level = 'warning')
	})
	dev.off()
}

#### oncostrip
if (plots$oncostrip) {
	logger('## Plotting oncostrip ...')
	oncostripfile = file.path(outdir, 'oncostrip.png')
	params        = mtparams$oncostrip
	params$maf    = laml
	if (!'genes' %in% names(params)) {
		params$genes = genes
	}
	do.call(png, c(list(filename = oncostripfile), devpars2))
	tryCatch({
		do.call(oncostrip, params)
	}, error = function(e) {
		log2pyppl('Failed to plot oncostrip:', e, level = 'warning')
	})
	dev.off()
}

#### titv
if (plots$titv) {
	logger('## Plotting titv ...')
	titvplot = file.path(outdir, 'titv.png')
	titvobj  = titv(maf = laml, plot = FALSE, useSyn = TRUE)
	params   = mtparams$titv
	params$res = titvobj
	do.call(png, c(list(filename = titvplot), devpars))
	tryCatch({
		do.call(plotTiTv, params)
	}, error = function(e) {
		log2pyppl('Failed to plot titv:', e, level = 'warning')
	})
	dev.off()
}

#### lollipop
if (plots$lollipop) {
	logger('## Plotting lollipops ...')

	params     = mtparams$lollipop
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
		lpdir  = file.path(outdir, 'lollipops')
		dir.create(lpdir, showWarnings = F)
		lpplot = file.path(lpdir, paste0(ps$gene, '.lollipop.png'))
		do.call(png, c(list(filename = lpplot), devpars2))
		tryCatch({
			print(do.call(lollipopPlot, ps))
		}, error = function(e) {
			log2pyppl('Failed to plot lollipop for', gene, ':', e, level = 'warning')
		})
		dev.off()
	}

	if (nthread == 1) {
		for (gene in lpgenes) {
			lollipopPlotSingle(gene)
		}
	} else {
		library(doParallel)
		cl = makeCluster(nthread)
		registerDoParallel(cl)
		foreach(i=1:length(lpgenes), .verbose = T, .packages=c('maftools')) %dopar% {
			lollipopPlotSingle(lpgenes[i])
		}
		stopCluster(cl)
	}
}

#### cbsseg
if (plots$cbsseg) {
	logger('## Plotting cbsseg ...')
	params  = c(list(maf = laml), mtparams$cbsseg)
	sams    = ifelse (!'tsb' %in% names(params), samples, params$tsb)

	cbssegSingle = function(sam) {
		segfiles = c(
			Sys.glob(file.path(indir, paste0(sam, '*.seg.txt'))),
			Sys.glob(file.path(indir, paste0(gsub('-', '.', sam), '*.seg.txt'))),
			Sys.glob(file.path(indir, paste0(gsub('_', '.', sam), '*.seg.txt')))
		)
		if (length(segfiles) == 0) {
			logger('   No seg file found for sample:', sam ,', skip.')
		} else {
			segdir  = file.path(outdir, 'cbssegs')
			dir.create(segdir, showWarnings = F)
			segplot = file.path(segdir, paste0(sam, '.cbsseg.png'))
			do.call(png, c(list(filename = segplot), devpars2))
			tryCatch({
				do.call(plotCBSsegments, c(list(cbsFile = segfiles[1]), params))
			}, error = function(e) {
				log2pyppl('Failed to plot cbsseg:', e, level = 'warning')
			})
			dev.off()
		}
	}

	if (nthread == 1) {
		for (sam in sams) {
			cbssegSingle(sam)
		}
	} else {
		library(doParallel)
		cl = makeCluster(nthread)
		registerDoParallel(cl)
		foreach(i=1:length(sams), .verbose = T, .packages=c('maftools')) %dopar% {
			cbssegSingle(sams[i])
		}
		stopCluster(cl)
	}
}

#### rainfall
if (plots$rainfall) {
	logger('## Plotting rainfall ...')
	params = c(list(maf = laml), mtparams$rainfall)
	if ('tsb' %in% names(params)) {
		sams = params$tsb
		prams$tsb = NULL
	} else {
		sams = samples
	}

	rainfallSingle = function(sam) {
		rfdir  = file.path(outdir, 'rainfalls')
		dir.create(rfdir, showWarnings = F)
		rfplot = file.path(rfdir, paste0(sam, '.rainfall.png'))
		do.call(png, c(list(filename = rfplot), devpars2))
		tryCatch({
			print(do.call(rainfallPlot, c(list(tsb = sam), params)))
		}, error = function(e){
			if (!params$detectChangePoints) {
				log2pyppl('Failed to plot rainfall without detecting change points for sample:', sam, level = 'warning')
			} else {
				tryCatch({
					do.call(rainfallPlot, c(list(tsb = sam, detectChangePoints = F), params))
				}, error = function(e){
					log2pyppl('Failed to plot rainfall without detecting change points for sample', sam, e, level = 'warning')
				})
			}
		})
		dev.off()
	}
	if (nthread == 1) {
		for (sam in sams) {
			rainfallSingle(sam)
		}
	} else {
		library(doParallel)
		cl = makeCluster(nthread)
		registerDoParallel(cl)
		foreach(i=1:length(sams), .verbose = T, .packages=c('maftools')) %dopar% {
			rainfallSingle(sams[i])
		}
		stopCluster(cl)
	}
}

#### tcgacomp
if (plots$tcgacomp) {
	logger('## Plotting tcgacompare ...')
	tcgacompplot = file.path(outdir, 'tcgacomp.png')
	do.call(png, c(list(filename = tcgacompplot), devpars2))
	params = c(list(maf = laml), mtparams$tcgacomp)
	if (!'cohortName' %in% names(params)) {
		params$cohortName = unlist(strsplit(basename(maffile), '.', fixed = T))[1]
	}
	tryCatch({
		do.call(tcgaCompare, params)
	}, error = function(e){
		log2pyppl('Failed to plot tcgacomp:', e, level = 'warning')
	})
	dev.off()
}

#### vaf
if (plots$vaf) {
	logger('## Plotting vaf ...')
	params = c(list(maf = laml), mtparams$vaf)
	if (!'vafCol' %in% names(params)) {
		log2pyppl('No vafCol(args.params.vaf.vafCol) provided, skip plotting VAF.', level = 'warning')
	} else {
		vafplot = file.path(outdir, 'vaf.png')
		do.call(png, c(list(filename = vafplot), devpars))
		tryCatch({
			do.call(plotVaf, params)
		}, error = function(e){
			log2pyppl('Failed to plot vaf:', e, level = 'warning')
		})
		dev.off()
	}
}

#### genecloud
if (plots$genecloud) {
	logger('## Plotting genecloud ...')
	gcplot = file.path(outdir, 'genecloud.png')
	do.call(png, c(list(filename = gcplot), devpars))
	params = c(list(input = laml), mtparams$genecloud)
	tryCatch({
		do.call(geneCloud, params)
	}, error = function(e){
		log2pyppl('Failed to plot genecloud:', e, level = 'warning')
	})
	dev.off()
}


#### gisticGenome
if (plots$gisticGenome) {
	logger('## Plotting gisticGenome ... ')
	if (!hasGistic) {
		logger('   No gistic files found, skip')
	} else {
		gisgplot = file.path(outdir, 'gisticGenome.png')
		do.call(png, c(list(filename = gisgplot), devpars2))
		params = c(list(gistic = lamlGistic), mtparams$gisticGenome)
		tryCatch({
			do.call(gisticChromPlot, params)
		}, error = function(e){
			log2pyppl('Failed to plot gisticGenome:', e, level = 'warning')
		})
		dev.off()
	}
}

#### gisticBubble
if (plots$gisticBubble) {
	logger('## Plotting gisticBubble ...')
	if (!hasGistic) {
		logger('   No gistic files found, skip\n')
	} else {
		gisbplot = file.path(outdir, 'gisticBubble.png')
		do.call(png, c(list(filename = gisbplot), devpars))
		params = c(list(gistic = lamlGistic), mtparams$gisticBubble)
		tryCatch({
			do.call(gisticBubblePlot, params)
		}, error = function(e){
			log2pyppl('Failed to plot gisticBubble:', e, level = 'warning')
		})
		dev.off()
	}
}

#### gisticOncoplot
if (plots$gisticOncoplot) {
	logger('## Plotting gisticOncoplot ...')
	if (!hasGistic) {
		logger('   No gistic files found, skip')
	} else {
		gisoncoplotfile = file.path(outdir, 'gisticOncoplot.png')
		params = c(list(gistic = lamlGistic), mtparams$gisticOncoplot)
		if (!is.null(annofile)) {
			params$clinicalData     = getClinicalData(x = laml)
			params$sortByAnnotation = TRUE
			if (!'top' %in% names(params))
				params$top = ngenes
		}
		do.call(png, c(list(filename = gisoncoplotfile), devpars2))
		tryCatch({
			do.call(gisticOncoPlot, params)
		}, error = function(e) {
			if (is.null(annofile)) {
				log2pyppl('Cannot generate gisticOncoplot without clinic features, skip.', level = 'warning')
			} else {
				tryCatch({
					params$clinicalFeatures = NULL
					do.call(gisticOncoPlot, params)
				}, error = function(e) {
					log2pyppl('Failed to plot gisticOncoplot:', e, level = 'warning')
				})
			}
		})
		dev.off()
	}
}

#### somInteraction
if (plots$somInteraction) {
	logger('## Plotting somInteraction ...')
	somInteractionfile = file.path(outdir, 'somInteraction.png')
	params = c(list(maf = laml), mtparams$somInteraction)
	if (!'top' %in% names(params) && !'genes' %in% names(params)){
		#params$top = ngenes
		params$genes = genes
	} else if ('top' %in% names(params)) {
		params$genes = genes[1:params$top]
	}
	do.call(png, c(list(filename = somInteractionfile), devpars))
	tryCatch(
		{
			do.call(somaticInteractions, params)
		}, error = function(e) {
			log2pyppl('Failed to plot somaticInteractions:', e, level = 'warning')
		}
	)
	dev.off()
}

#### oncodrive
if (plots$oncodrive) {
	logger('## Plotting oncodrive ...')
	oncodrivefile = file.path(outdir, 'oncodrive.png')
	OcParams = list(maf = laml, AACol = NULL, minMut = 5, pvalMethod = "zscore",
					nBgGenes = 100, bgEstimate = TRUE, ignoreGenes = NULL)
	params = mtparams$oncodrive
	for (name in names(params)) {
		if (name %in% names(OcParams))
			OcParams[[name]] = params[[name]]
	}
	tryCatch(
		{
			sig = do.call(oncodrive, OcParams)

			PocParams = list(res = sig, fdrCutOff = 0.05, useFraction = FALSE,
							colCode = NULL) #, labelSize = 2) # plotOncodrive hangs with labelSize argument
			for (name in names(params)) {
				if (name %in% names(PocParams))
					PocParams[[name]] = params[[name]]
			}
			do.call(png, c(list(filename = oncodrivefile), devpars))
			do.call(plotOncodrive, PocParams)
			dev.off()
		}, error = function(e) {
			log2pyppl('Failed to plot oncodrive:', e, level = 'warning')
		}
	)
}

#### pfam
if (plots$pfam) {
	logger('## Plotting pfam ...')
	pfamfile = file.path(outdir, 'pfam.png')
	params = c(list(maf = laml), mtparams$pfam)
	if (!'top' %in% names(params))
		params$top = ngenes
	do.call(png, c(list(filename = pfamfile), devpars))
	tryCatch(
		{
			pfam = do.call(pfamDomains, params)
			pfsum = pfam$proteinSummary
			# add link to domains 
			pfsum$DomainLabel = sapply(
				pfsum$DomainLabel, 
				function(x) {sprintf("[%s](https://www.ncbi.nlm.nih.gov/cdd/?term=%s)", x, x)})
			write.table(
				pfsum[!is.na(pfsum$DomainLabel), 1:7, with=FALSE], 
				file.path(outdir, 'pfam.csv'),
				row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

		}, error = function(e) {
			log2pyppl('Failed to plot pfam:', e, level = 'warning')
		}
	)
	dev.off()
}

#### pancan
if (plots$pancan) {
	if (!is.null(siggene)) {
		params  = mtparams$pancan
		params$mutsigResults = siggene
		if (!'cohortName' %in% names(params)) {
			params$cohortName = unlist(strsplit(basename(maffile), '.', fixed = T))[1]
		}
		if (!'inputSampleSize ' %in% names(params)) {
			params$inputSampleSize  = nsample
		}
		pancanfile = file.path(outdir, 'pancan.png')
		do.call(png, c(list(filename = pancanfile), devpars))
		tryCatch({
			do.call(pancanComparison, params)
		}, error = function(e){
			log2pyppl('Failed to plot pancan:', e, level = 'warning')
		})
		dev.off()
	}
}

#### survival
if (plots$survival) {
cat('## Plotting survivals ...\n')

	params        = mtparams$survival
	params$maf    = laml
	params$isTCGA = isTCGA
	if (!'time' %in% names(params) || !'Status' %in% names(params)) {
		cat('-  No time/Status column specified for survival analysis!')
	} else {
		survivalSingle = function(gene) {
			ps = params
			ps$genes = gene
			svdir  = file.path(outdir, 'survivals')
			dir.create(svdir, showWarnings = F)
			svplot = file.path(svdir, paste0(ps$genes, '.survival.png'))
			do.call(png, c(list(filename = svplot), devpars))
			tryCatch({
				do.call(mafSurvival, ps)
			}, error = function(e){
				log2pyppl('Failed to plot survival for gene', gene, ':', e, level = 'warning')
			})
			dev.off()
		}
		if ('genes' %in% names(params)) {
			svgenes = params$genes
			params$genes = NULL
		} else {
			svgenes = genes
		}

		if (nthread == 1) {
			for (gene in svgenes) {
				survivalSingle(gene)
			}
		} else {
			library(doParallel)
			cl = makeCluster(nthread)
			registerDoParallel(cl)
			foreach(i=1:length(svgenes), .verbose = T, .packages=c('maftools')) %dopar% {
				survivalSingle(svgenes[i])
			}
			stopCluster(cl)
		}
	}
}

#### heterogeneity
if (plots$heterogeneity) {
	logger('## Plotting heterogeneity ...')
	if ('tsb' %in% names(params)) {
		sams = params$tsb
		params$tsb = NULL
	} else {
		sams = samples
	}

	ihParams = list(maf = laml, top = 5, vafCol = NULL,
	ignChr = NULL, minVaf = 0, maxVaf = 1,
	dirichlet = FALSE)

	params = mtparams$heterogeneity
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
		log2pyppl('No vafCol(args.params.heterogeneity.vafCol) provided, skip plotting heterogeneity.', level = 'warning')
	} else {
		heteroSingle = function(sam) {
			segfiles = c(
				Sys.glob(file.path(indir, paste0(sam, '*.seg.txt'))),
				Sys.glob(file.path(indir, paste0(gsub('-', '.', sam), '*.seg.txt'))),
				Sys.glob(file.path(indir, paste0(gsub('_', '.', sam), '*.seg.txt')))
			)
			tryCatch({
				if (length(segfiles) == 0) {
					hetero = do.call(inferHeterogeneity, c(list(tsb = sam), ihParams))
				} else {
					hetero = do.call(inferHeterogeneity, c(list(tsb = sam, segFile = segfiles[1]), ihParams))
				}
				hgdir  = file.path(outdir, 'heterogeneities')
				dir.create(hgdir, showWarnings = F)
				heterogeneityfile = file.path(hgdir, paste0(sam, '.heterogeneity.png'))
				do.call(png, c(list(filename = heterogeneityfile), devpars))
				do.call(plotClusters, c(list(clusters = hetero, tsb = sam), PocParams))
				dev.off()
			}, error = function(e) {
				log2pyppl('Failed to plot heterogeneity for sample', sam, ':', e, level = 'warning')
			})
		}

		if (nthread == 1) {
			for (sam in sams) {
				heteroSingle(sam)
			}
		} else {
			library(doParallel)
			cl = makeCluster(nthread)
			registerDoParallel(cl)
			foreach(i=1:length(sams), .verbose = T, .packages=c('maftools')) %dopar% {
				heteroSingle(sams[i])
			}
			stopCluster(cl)
		}
	}
}

#### signature
if (plots$signature) {
	logger('## Plotting signature ...')
	require('NMF')
	params = mtparams$signature
	tmParams = list(prefix = '', add = TRUE, ignoreChr = NULL, useSyn = TRUE, fn = NULL)
	for (name in names(params)) {
		if (name %in% names(tmParams)) {
			tmParams[[name]] = params[[name]]
			params  [[name]] = NULL
		}
	}
	tmParams$ref_genome = paste0('BSgenome.Hsapiens.UCSC.', genome)
	tmParams$maf = laml
	tryCatch({
		tm = do.call(trinucleotideMatrix, tmParams)
	}, error = function(e){
		tmParams$add <<- FALSE
		tryCatch({
			tm <<- do.call(trinucleotideMatrix, tmParams)
		}, error = function(e){
			tmParams$prefix = NULL
			tm <<- do.call(trinucleotideMatrix, tmParams)
		})
	})
	# apobec
	apobecfile = file.path(outdir, 'apobec.png')
	do.call(png, c(list(filename = apobecfile), devpars))
	tryCatch(
		{
			plotApobecDiff(tnm = tm, maf = laml)
		}, error = function(e) {
			log2pyppl('Failed to plot apobec signature:', e, level = 'warning')
		}
	)
	dev.off()

	# sigs
	sgParams = list(n = NULL, nTry = 6, pConstant = 1, plotBestFitRes = FALSE, parallel = NULL)
	for (name in names(params)) {
		if (name %in% names(sgParams)) {
			sgParams[[name]] = params[[name]]
			params  [[name]] = NULL
		}
	}
	sgParams$mat = tm
	tryCatch({
		cophenfile = file.path(outdir, 'signature-cophen.png')
		do.call(png, c(list(filename = cophenfile), devpars))
		dev.off()
		sigs  = do.call(extractSignatures, sgParams)
		sigfile = file.path(outdir, 'signature.png')
		do.call(png, c(list(filename = sigfile), devpars))
		plotSignatures(sigs)
		dev.off()
		# cosine similarity against validated signatures
		plot.heatmap(sigs$coSineSimMat, file.path(outdir, 'signature-sim.png'), devpars = devpars, params = list(dendro = 'col'))
		# signature contrib
		sigcfile = file.path(outdir, 'signature-contrib.png')
		do.call(png, c(list(filename = sigcfile), devpars))
		plotSignatures(nmfRes = sigs, contributions = TRUE)
		dev.off()
	}, error = function(e) {
		log2pyppl('Failed to plot extract signatures:', e, level = 'warning')
	})


}





