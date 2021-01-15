{{rimport}}('__init__.r', 'plot.r')
library(caret)
options(stringsAsFactors = FALSE)

{% python from pyppl.utils import alwaysList %}
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
# I tried na.action for train function, but I got an error:
# Error in as.character(call_obj$na.action) :
#	cannot coerce type 'closure' to vector of type 'character'
# haven't figured out a way to fixed, just remove all NA records
indata = indata[complete.cases(indata),,drop = FALSE]
ycol = all.vars(trainps$form)[1]
# for glm or binary output
if (!is.numeric(indata[, ycol])) {
	indata[, ycol] = as.factor(indata[, ycol])
}
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
	tryCatch({
		save.plot(ggplot(m), modelplot, devpars = devpars)
	}, error = function(e) {
		logger(e, level = 'WARNNING')
	})
}

klass = m$pred$obs
if (is.null(m$pred$Resample)) {
	m$pred$Resample = 'M'
}
if (length(levels(factor(klass))) > 2) {
	rocdata = data.frame(D = klass >= .5, M = m$pred[,4], name = m$pred$Resample)
} else {
	rocdata = data.frame(D = klass, M = m$pred[,4], name = m$pred$Resample)
}
aucs = plot.roc(rocdata, 'return',
	params = list(returnTable = TRUE),
	ggs    = list(scale_color_discrete = list(guide = FALSE)),
	stacked = TRUE, devpars = devpars)
aucs$table$mean_auc = mean(aucs$table$auc)
aucfile = sprintf('%s.aucs.txt', prefix)
write.table(
	pretty.numbers2(aucs$table[aucs$table$is.best,, drop=FALSE], threshold..sensitivity..specificity..auc..auc_95ci1..auc_95ci2..mean_auc = '%.3f'),
	aucfile, row.names = FALSE, quote = FALSE, sep = "\t")

if ('roc' %in% plots) {
	rocplot = sprintf('%s.roc.png', prefix)
	save.plot(aucs$plot, rocplot, devpars)
}

vi = varImp(m)
vifile = sprintf('%s.varimp.txt', prefix)
write.table(vi$importance, vifile, quote = FALSE, sep = "\t")
if ('varimp' %in% plots) {
	viplot = sprintf('%s.varimp.png', prefix)
	save.plot(ggplot(vi), viplot, devpars = devpars)
}
