{{rimport}}('__init__.r', 'sampleinfo.r')
library(clonevol)

mutfile = {{i.mutfile | quote}}
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
saminfo = SampleInfo2$new(saminfo)

samples = saminfo$get.samples()
sample.groups = c()
for (sam in samples) {
	sample.groups[[sam]] = saminfo$sample.info(sam, info = 'Group')
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
head(mutinfo)
print(samples)
print(sample.groups)
y = do.call(infer.clonal.models,
			c(list(variants = mutinfo,
				   vaf.col.names = samples,
				   sample.groups = sample.groups),
			params.infer.clonal.models))

# transfer.events.to.consensus.trees
logger('Running transfer.events.to.consensus.trees ...')
y = do.call(transfer.events.to.consensus.trees,
			c(list(x = y,
				   events = mutinfo[mutinfo$is.driver,]),
			params.transfer.events.to.consensus.trees))

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
