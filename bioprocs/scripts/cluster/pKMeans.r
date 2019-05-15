{{rimport}}('__init__.r')

infile  = {{ i.infile | quote}}
outfile = {{ o.outfile | quote}}
outdir  = {{ o.outdir | quote}}
inopts  = {{ args.inopts | R}}
k       = {{ args.k | R}}
doplot  = {{args.plot | R}}
devpars = {{args.devpars | R}}
ggs     = {{args.ggs | R}}

if (!is.numeric(k)) {
	stop('Number of clusters `k` must be an integer.')
}

if (doplot) {
	{{rimport}}('plot.r')
}

indata  = read.table.inopts(infile, inopts)
samples = rownames(indata)
if (is.null(samples)) {
	samples = 1:nrow(indata)
}

cluster = kmeans(indata, k)
# Sample,Cluster,Center,Withinss,Mean_Feature1,Mean_Feature2,...
ret = NULL
for (i in 1:k) {
	indexes = which(cluster$cluster == i)
	cSample = samples[indexes]
	tmp = data.frame(
		Cluster  = as.factor(rep(i, length(cSample))),
		Withinss = cluster$withinss[i]
	)
    cmeans = t(colMeans(indata[cSample, ]))
    colnames(cmeans) = sapply(colnames(cmeans), function(s) paste0("Mean_", s))
	tmp = cbind(tmp, cmeans)
    rownames(tmp) = cSample
	if (is.null(ret)) {
		ret = tmp
	} else {
		ret = rbind(ret, tmp)
	}
}
ret = ret[samples,]
write.table(ret, outfile, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

if (doplot) {
	plotdata = cbind(indata, ret)
	plotfile = file.path(
		outdir,
		paste0(
			tools::file_path_sans_ext(tools::file_path_sans_ext(basename(outfile))),
			'.png'
		)
	)
	ggs = c(list(
		geom_point = list(
			mapping = aes(x = ret[,3], y = ret[,4], color = 'Center', shape = 'Center'),
			size = 4
		)
	), ggs)
	if (ncol(indata) == 2) {
		plot.scatter(
			plotdata,
			plotfile,
			x       = 1,
			y       = 2,
			params  = list(mapping = aes_string(color = 'Cluster', shape = 'Cluster')),
			ggs     = ggs,
			devpars = devpars
		)
	} else {
		stop('Not implemented yet.')
	}
}
