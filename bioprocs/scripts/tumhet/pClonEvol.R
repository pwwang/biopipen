{{rimport}}('__init__.r', 'sampleinfo.r')
library(clonevol)

mutfile = {{i.infile | quote}}
saminfo = {{i.saminfo | quote}}
outdir  = {{o.outdir | quote}}
inopts  = {{args.inopts | R}}
devpars = {{args.devpars | R}}

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

mutinfo = read.table.inopts(mutfile, inopts)
saminfo = SampleInfo2(saminfo)

samples = saminfo$get.samples()
sample.groups = c()
for (sam in samples) {
	sample.groups[[sam]] = saminfo$sample.info(sam, info = 'Group')
}

vaf.cols = paste(samples, 'vaf', sep = '.')

# plot.variant.clusters
variant.cluster.outfile = file.path(outdir, 'variant_cluster.png')
do.call(png, c(list(filename = variant.cluster.outfile), devpars))
do.call(plot.variant.clusters,
		c(list(df = mutinfo, vaf.col.names = vaf.cols), params.plot.variant.clusters))
dev.off()

# plot.cluster.flow
cluster.flow.outfile = file.path(outdir, 'cluster_flow.png')
do.call(png, c(list(filename = cluster.flow.outfile), devpars))
do.call(plot.cluster.flow,
		c(list(variants = mutinfo, vaf.col.names = vaf.cols), params.plot.cluster.flow))
dev.off()

# infer.clonal.models
y = do.call(infer.clonal.models,
			c(list(variants = mutinfo,
				   vaf.col.names = vaf.cols,
				   sample.groups = sample.groups),
			params.infer.clonal.models))

# convert.consensus.tree.clone.to.branch
y = do.call(convert.consensus.tree.clone.to.branch,
			c(list(x=y), params.consensus.tree.clone.to.branch))

do.call(plot.clonal.models,
		c(list(x=y, out.dir = outdir, out.format = 'png'), params.plot.clonal.models))






