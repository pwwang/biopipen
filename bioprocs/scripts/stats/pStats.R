{{rimport}}('__init__.r', 'plot.r')
library(gridExtra)
# usage: https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html

options(stringsAsFactors=FALSE)
infile = {{i.infile | R}}
outdir = {{o.outdir | R}}
inopts = {{args.inopts | R}}
types  = {{args.types | R}}
groups = {{args.groups | R}}
ignore = {{args.ignore | R}}

indata = read.table.inopts(infile, inopts)
ninst  = nrow(indata)

# top 10 records of the data
write.table(
	indata[1:10],
	file.path(outdir, {{i.infile | stem | @append: ".top10.txt" | quote}}),
	row.names = inopts$rnames,
	col.names = inopts$cnames,
	sep = inopts$delimit,
	quote = FALSE)

do_continuous = function (feat, groups, outdir) {
	#histp = plot.histo(data, 'return')
	statfile = file.path(outdir, paste0('feature-stat.', feat, '.txt'))
	testfile = file.path(outdir, paste0('feature-test.', feat, '.txt'))
	statdata = data.frame()
	fdata = as.numeric(indata[, feat])
	statdata = rbind(statdata, list(
		Group  = '__ALL__',
		N      = ninst,
		Max    = max(fdata, na.rm = TRUE),
		Min    = min(fdata, na.rm = TRUE),
		Mean   = mean(fdata, na.rm = TRUE),
		Median = median(fdata, na.rm = TRUE),
		Sd     = sd(fdata, na.rm = TRUE)
	))
	testdata = data.frame()
	for (grup in groups) {
		lvls = levels(factor(indata[, grup]))
		for (lvl in lvls) {
			lvldata = as.numeric(indata[indata[, grup] == lvl, feat])
			statdata = rbind(statdata, list(
				Group  = paste0(feat, ' [', lvl, ']'),
				N      = length(lvldata),
				Max    = max(lvldata, na.rm = TRUE),
				Min    = min(lvldata, na.rm = TRUE),
				Mean   = mean(lvldata, na.rm = TRUE),
				Median = median(lvldata, na.rm = TRUE),
				Sd     = sd(lvldata, na.rm = TRUE)
			))
		}
	}
	write.table(statdata, statfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	if (nrow(testdata)) {
		write.table(testdata, testfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	}
}

do_categorical = function(feat, groups, outdir) {}

features = colnames(indata)
features = setdiff(features, ignore)
# get the number of cases for each group
# if (is.true(groups)) {
# 	grpfile = file.path(outdir, 'groups.txt')
# 	grpdata = data.frame(`Group` = character(), `Number of cases` = character())
# 	for (grp in groups) {
# 		lvls = levels(factor(indata[, grp]))
# 		grprow = sapply(lvls, function(lvl) {
# 			paste(lvl, ':', nrow(indata[indata[, grp] == lvl, ]))
# 		})
# 		grprow = data.frame(`Group` = grp, `Number of cases` = paste(grprow, collapse = ' | '))
# 		grpdata = rbind(grpdata, grprow)
# 	}
# 	write.table(grpdata, grpfile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# }
for (feat in features) {
	ftype = list.get(types, feat, 'auto')
	if (ftype == 'auto') {
		lvls  = levels(factor(indata[, feat]))
		if (length(lvls) < 10 && length(lvls) < ninst * .8) {
			ftype = 'categorical'
		} else {
			ftype = 'continuous'
		}
	}
	ftype = ifelse(startsWith(ftype, 'cont'), 'continuous', 'categorical')
	do.call(paste0('do_', ftype), list(feat = feat, groups = groups, outdir = outdir))
}


