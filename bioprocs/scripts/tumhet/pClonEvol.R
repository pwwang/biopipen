{{rimport}}('__init__.r', 'sampleinfo.r')
library(clonevol)

# use to save converted pyclone results
indir     = {{job.indir | quote}}
mutfile   = {{i.mutfile | quote}}
indrivers = {{i.drivers | quote}}
samfile   = {{i.samfile | quote}}
outdir    = {{o.outdir | quote}}
inopts    = {{args.inopts | R}}
devpars   = {{args.devpars | R}}
drivers   = {{args.drivers | R}}
refgene   = {{args.refgene | R}}
bedtools  = {{args.bedtools | R}}

if (is.true(indrivers)) {
	if (file.exists(indrivers)) {
		drvfile = indrivers
		drvtype = 'plain'
		if (dir.exists(indrivers)) {
			drvfile = Sys.glob(file.path(indrivers, "*.sig_genes.txt"))[1]
			drvtype = 'mutsig'
		} else if (endsWith(indrivers, '.sig_genes.txt')) {
			drvtype = 'mutsig'
		}

		if (drvtype == 'plain') {
			drivers = read.table.inopts(drvfile, list(cnames = FALSE, rnames = FALSE))
			drivers = unlist(drivers[, 1])
		} else {
			drivers = read.table.inopts(drvfile, list(cnames = TRUE, rnames = FALSE))
			drivers = drivers[1:10, , drop = FALSE]
			drivers = unlist(drivers[drivers$p < 0.05, 1, drop = TRUE])
		}
	} else {
		drivers = unlist(strsplit(indrivers, ",", fixed = TRUE))
	}
}

params.plot.variant.clusters =
	{{args.params['plot.variant.clusters'] | R}}
params.plot.cluster.flow =
	{{args.params['plot.cluster.flow'] | R}}
params.infer.clonal.models =
	{{args.params['infer.clonal.models'] | R}}
params.transfer.events.to.consensus.trees =
	{{args.params['transfer.events.to.consensus.trees'] | R}}
params.convert.consensus.tree.clone.to.branch =
	{{args.params['convert.consensus.tree.clone.to.branch'] | R}}
params.plot.clonal.models =
	{{args.params['plot.clonal.models'] | R}}

saminfo = SampleInfo2$new(samfile)

samples = saminfo$get.samples()
sample.groups = c()
for (sam in samples) {
	sample.groups[[sam]] = saminfo$sample.info(sam, info = 'Group')[1]
}

# expect mutfile:
# cluster gene is.driver P.vaf R.vaf P.ccf R.ccf P.ref.count P.var.count P.depth R.ref.count R.var.count R.depth P R
# 1 SMC3 TRUE 50.28 49.13 100.56 98.26 3027 2611 5639 2035 1990 4024 50.28 49.13
# 1 PTPRT TRUE 45.24 45.31 90.48 90.62 2159 1799 3958 1715 1579 3294 45.24 45.31
# 1 - FALSE 43.83 41.37 87.66 82.74 2931 2364 5296 2224 1501 3725 43.83 41.37

# in case mutfile is pylone output
if (dir.exists(mutfile) && endsWith(mutfile, '.pyclone')) {
	# tables/loci.tsv
	#mutation_id	sample_id	cluster_id	cellular_prevalence	cellular_prevalence_std	variant_allele_frequency	gene
	#chr10:17875816	NA12878	2	0.828865827591	0.0355432431918	0.904883227176	Gene1
	#chr10:17875816	NA12156	2	0.741264908504	0.11050033704	0.560311284047	Gene2

	# <sample>.mut.beds
	##CHR   START   END     NAME            GT  REF VAR CASE    VAF
	#chr2	3482633	3482633	chr2:3482633	AB	103	2	sample1	0.019
	#chr2	9661972	9661972	chr2:9661972	AB	60	6	sample1	0.091
	logger('Extracting information from pyclone results ...')
	mutations = read.table.inopts(file.path(mutfile, "tables/loci.tsv"), list(rnames = FALSE))
	mutations = mutations[, c(
		"mutation_id", "sample_id", "cluster_id", "variant_allele_frequency", "cellular_prevalence"),
		drop = FALSE]
	rownames(mutations) = make.unique(do.call(paste, c(mutations[c("sample_id", "mutation_id")], sep = "-")))

	samples = levels(as.factor(mutations$sample_id))
	mutbeds = NULL
	for (sample in samples) {
		mutbed = file.path(mutfile, paste0(sample, '.muts.bed'))
		mutbed_withgene = file.path(outdir, paste0(sample, '.muts.bed'))
		runcmd(sprintf("%s intersect -a '%s' -b '%s' -loj > '%s'", bedtools, mutbed, refgene, mutbed_withgene))

		mbtmp = read.table.inopts(mutbed_withgene, list(cnames = FALSE, rnames = FALSE))
		mbtmp = mbtmp[, c(1:9, 18), drop = FALSE]
		colnames(mbtmp) = c("CHR", "START", "END", "NAME", "GENOTYPE",
							"REFCOUNTS", "VARCOUNTS", "CASE", "VARFREQ", "GENE")
		mbtmp$GENE = sapply(mbtmp$GENE, function(g) {
			if (g == '.') {
				'-'
			} else {
				#gene_id "DCDC2C"
				geneid = unlist(strsplit(g, "; ", fixed = TRUE))[1]
				substring(geneid, 10, nchar(geneid)-1)
			}
		})
		mbtmp$SAMPLE = sample
		mbtmp$DEPTH  = mbtmp$VARCOUNTS / mbtmp$VARFREQ
		if (is.null(mutbeds)) {
			mutbeds = mbtmp
		} else {
			mutbeds = rbind(mutbeds, mbtmp)
		}
	}
	rownames(mutbeds) = make.unique(do.call(paste, c(mutbeds[c("SAMPLE", "NAME")], sep = "-")))
	mutbeds = mutbeds[rownames(mutations),,drop = FALSE]

	mutations = cbind(mutations, mutbeds)
	mutations = mutations[, c(
		"cluster_id", "variant_allele_frequency", "cellular_prevalence",
		"REFCOUNTS", "VARCOUNTS", "DEPTH", "SAMPLE", "mutation_id", "GENE"
	), drop = FALSE]

	sample0 = samples[1]
	mutinfo = data.frame(
		#Cluster identities must be contiguous integer starting at 1
		cluster   = mutations[mutations$SAMPLE == sample0, "cluster_id"] + 1,
		gene      = mutations[mutations$SAMPLE == sample0, "GENE"],
		is.driver = mutations[mutations$SAMPLE == sample0, "GENE"] %in% drivers)
	rownames(mutinfo) = make.unique(mutations[mutations$SAMPLE == sample0, "mutation_id"])
	for (sample in samples) {
		samdata = mutations[mutations$SAMPLE == sample,,drop = FALSE]
		rownames(samdata) = make.unique(samdata$mutation_id)
		samdata = samdata[rownames(mutinfo),,drop = FALSE]
		mutinfo[[sample]] = samdata$variant_allele_frequency * 100
		mutinfo[[paste(sample, 'vaf', sep = '.')]] = samdata$variant_allele_frequency * 100
		mutinfo[[paste(sample, 'ccf', sep = '.')]] = samdata$cellular_prevalence * 100
		mutinfo[[paste(sample, 'ref.count', sep = '.')]] = samdata$REFCOUNTS
		mutinfo[[paste(sample, 'var.count', sep = '.')]] = samdata$VARCOUNTS
		mutinfo[[paste(sample, 'depth', sep = '.')]] = samdata$DEPTH
	}
	write.table(mutinfo, file.path(indir, paste0(basename(mutfile), '.txt')), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
} else {
	mutinfo = read.table.inopts(mutfile, inopts)
}

# plot.variant.clusters
logger('Running plot.variant.clusters ...')
variant.cluster.outfile = file.path(outdir, 'variant_cluster.png')
do.call(png, c(list(filename = variant.cluster.outfile), devpars))
do.call(plot.variant.clusters,
		c(list(df = mutinfo, vaf.col.names = samples), params.plot.variant.clusters))
dev.off()

# plot.cluster.flow
logger('Running plot.cluster.flow ...')
cluster.flow.outfile = file.path(outdir, 'cluster_flow.png')
do.call(png, c(list(filename = cluster.flow.outfile), devpars))
do.call(plot.cluster.flow,
		c(list(variants = mutinfo, vaf.col.names = samples), params.plot.cluster.flow))
dev.off()

# infer.clonal.models
logger('Running infer.clonal.models ...')
# determin founding cluster if it is NULL as default
# if (is.null(params.infer.clonal.models$founding.cluster)) {
# 	allVafs = aggregate(as.formula(paste(paste(samples, collapse = '+'), 'cluster', sep = '~')), mutinfo, sum)
# 	params.infer.clonal.models$founding.cluster = allVafs[which.max(allVafs[,2]), 'cluster'][1]
# }

passed = FALSE
# try to remove invalid samples, instead of stop the program
while (!passed) {
	tryCatch({
		capture = capture.output(y <<- do.call(infer.clonal.models,
			c(list(variants = mutinfo,
				   vaf.col.names = samples,
				   sample.groups = sample.groups),
			params.infer.clonal.models)), type = 'message')
		if (!grepl("ERROR: No clonal models for sample:", capture[1], fixed = TRUE)) {
			passed <<- TRUE
		} else {
			exclude.sample = sub(" $", "", substring(capture[1], 37))
			logger(paste('Sample', exclude.sample, 'not valid for modeling', level = 'warning'))
			samples = setdiff(samples, exclude.sample)
			if (len(samples) < 2) {
				stop('Not enough valid samples.')
			}
		}
	}, error = function(e) { passsed <<- FALSE })
}

if (any(mutinfo$is.driver)) {
	# transfer.events.to.consensus.trees
	logger('Running transfer.events.to.consensus.trees ...')
	y = do.call(transfer.events.to.consensus.trees,
				c(list(x = y,
					events = mutinfo[mutinfo$is.driver,]),
				params.transfer.events.to.consensus.trees))
}

# convert.consensus.tree.clone.to.branch
logger('Running convert.consensus.tree.clone.to.branch ...')
y = do.call(convert.consensus.tree.clone.to.branch,
			c(list(x=y), params.convert.consensus.tree.clone.to.branch))

logger('Running plot.clonal.models ...')
do.call(plot.clonal.models,
		c(list(y,
			   out.dir = outdir,
			   out.format = 'png',
			   resolution = devpars$res,
			   # let it choose
			   #width = devpars$width / 96 * 1.5,
			   #height = devpars$height / 96,
			   overwrite.output = TRUE), params.plot.clonal.models))
