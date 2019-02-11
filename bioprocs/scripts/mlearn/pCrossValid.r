{{rimport}}('__init__.r', 'plot.r')
library(caret)
options(stringsAsFactors = FALSE)

{% python from bioprocs.utils import alwaysList %}
{% python from os import path %}
infile   = {{ i.infile | quote}}
prefix   = {{ o.outdir, fn2(i.infile) | *path.join | quote }}
outmodel = {{ o.outmodel | quote}}
outdir   = {{ o.outdir | quote}}
inopts   = {{ args.inopts | R}}
seed     = {{ args.seed | R}}
plots    = {{ args.plots | alwaysList | R}}
ctrl     = {{ args.ctrl | R}}
nthread  = {{ args.nthread | R}}
trainps  = {{ args.train | R}}
devpars  = {{ args.devpars | R}}

set.seed(seed)

if (is.null(ctrl$method) || ctrl$method == '') {
	stop('Method has not been specified for trainControl (args.ctrl.method).')
}

if (is.null(trainps$method) || trainps$method == '') {
	stop('Method has not been specified for training (args.train.method).')
}

if (is.null(trainps$form) || trainps$form == '') {
	stop('Formula has not been specified for training (args.train.form).')
}
trainps$form = as.formula(trainps$form)

# method specific packages
# make it a complete list in the future
if (trainps$method %in% c('glm')) {
	library(e1071)
}

if (nthread > 1) {
	library(doParallel)
	cl = makePSOCKcluster(nthread)
	registerDoParallel(cl)
}

indata            = read.table.inopts(infile, inopts)
ctrl$verboseIter  = TRUE
train.control     = do.call(trainControl, ctrl)
trainps$data      = indata
trainps$trControl = train.control
m                 = do.call(train, trainps)

if (nthread > 1) {
	stopCluster(cl)
}
saveRDS(m, outmodel)

if ('model' %in% plots) {
	modelplot = paste0(prefix, '.model.png')
	do.call(png, c(list(modelplot), devpars))
	tryCatch({
		print(ggplot(m))
	}, error = function(e) {
		logger(e, level = 'WARNNING')
	}, finally = {
		dev.off()
	})
}

if ('roc' %in% plots) {
	klass = m$pred$obs
	if (is.null(m$pred$Resample)) {
		m$pred$Resample = 'M'
	}
	if (length(levels(factor(klass))) > 2) {
		klass = (klass - min(klass)) / (max(klass) - min(klass))
		# use different cutoff to bin
		for (tfcut in c(.5, .6, .7, .8, .9)) {
			rocplot = sprintf('%s.roc%d.png', prefix, 10*tfcut)
			rocdata = data.frame(D = klass > tfcut, M = m$pred[,4], name = m$pred$Resample)
			aucs = plot.roc(rocdata, rocplot, stacked = TRUE, devpars = devpars)
		}
	} else {
		rocplot = sprintf('%s.roc.png', prefix)
		rocdata = data.frame(D = klass, M = m$pred[,4], name = m$pred$Resample)
		aucs = plot.roc(rocdata, rocplot, stacked = TRUE, devpars = devpars)
	}
	aucfile = sprintf('%s.aucs.txt', prefix)
	aucs    = t(as.data.frame(aucs))
	colnames(aucs) = 'AUC'
	write.table(aucs, aucfile, quote = FALSE, sep = "\t")
}

if ('varimp' %in% plots) {
	vi = varImp(m)
	viplot = sprintf('%s.varimp.png', prefix)
	vifile = sprintf('%s.varimp.txt', prefix)
	write.table(vi$importance, vifile, quote = FALSE, sep = "\t")
	do.call(png, c(list(viplot), devpars))
	print(ggplot(vi))
	dev.off()
}



