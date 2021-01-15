{{'plot.r', '__init__.r' | rimport}}

infiles  = {{i.infiles | R}}
metafile = {{i.metafile | R}}
outfile  = {{o.outfile | R}}
inopts   = {{args.inopts | R}}
cats     = {{args.cats | ?!: "NULL"}}

allelements = c()
categories  = c()
datalist    = list()
for (infile in infiles) {
	indata = read.table.inopts(infile, inopts)
	if (is.function(cats)) {
		if (length(formals(cats)) == 1) {
			if (is.true(inopts$cnames)) {
				category = cats(colnames(indata)[1])
			} else {
				category = cats(tools::file_path_sans_ext(basename(infile)))
			}
		} else if (!is.true(inopts$cnames)) {
			stop('`args.cats` has 2 arguments but `inopts.cnames` is `False`.')
		} else {
			category = cats(tools::file_path_sans_ext(basename(infile)), colnames(indata)[1])
		}
	} else if (is.true(inopts$cnames)) {
		category = colnames(indata)[1]
	} else {
		category = tools::file_path_sans_ext(basename(infile))
	}
	categories  = make.unique(c(categories, category))
	catdata     = as.vector(unlist(indata[,1]))
	datalist[[category]] = catdata
	allelements = unique(c(allelements, catdata))
}
data = NULL
for (name in names(datalist)) {
	tmpdata = data.frame(tmp = as.integer(allelements %in% datalist[[name]]))
	colnames(tmpdata) = name
	if (!is.null(data)) {
		data = cbind(data, tmpdata)
	} else {
		data = tmpdata
	}
}
rownames(data) = allelements

tool    = {{args.tool | R}}
ncols   = ncol(data)
cnames  = colnames(data)
devpars = {{args.devpars | R}}
# use real venn plot
if ((ncols <= 3 && tool == 'auto') || tool == 'venn') {
	params = {{args.params | R}}
	plot.venn(data, outfile, params, devpars)
} else {
	metadata = data.frame(counts = colSums(data), ratio = colSums(data)/nrow(data))
	if (metafile != "") {
		exmeta   = read.table(metafile, sep = "\t", header = TRUE, row.names = 1, check.names = F)
		metadata = cbind.fill(metadata, exmeta)
	}
	metadata = cbind(metadata, sets = rownames(metadata))
	missingcols = rownames(metadata[metadata$counts == 0,,drop=F])
	if (length(missingcols) > 0) {
		log2pyppl('No data in columns:', missingcols)
	}
	metadata = metadata[metadata$counts>0, , drop = F]
	data = data[, rownames(metadata), drop = F]
	params = {{args.params | R}}
	plot.upset(data, outfile, params, devpars)
}

write.table(data, file.path(dirname(outfile), "data.txt"), sep = "\t", quote = FALSE)
