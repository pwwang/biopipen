{{rimport}}('__init__.r')
library(methods)
library(glmnet)
options(stringsAsFactors = FALSE)

inD     = {{i.D | quote}}
inPt    = {{i.Pt | quote}}
inY     = {{i.Y | quote}}
outW    = {{o.W | quote}}
outdir  = {{o.outdir | quote}}
seed    = {{args.seed | R}}
inopts  = {{args.inopts | R}}
svdP    = {{args.svdP | R}}
predY   = {{args.predY | R}}
WPt     = {{args.WPt | R}}
WtDtY   = {{args.WtDtY | R}}
nfold   = {{args.nfold | R}}
tfrac   = {{args.tfrac | R}}
nthread = {{args.nthread | R}}
method  = {{args.method | quote}}

set.seed(seed)

logger('Reading input files ...')
D  = as.matrix(read.table.inopts(inD, inopts))
Pt = as.matrix(read.table.inopts(inPt, inopts))
Y  = as.matrix(read.table.inopts(inY, inopts))

logger('Deciding training and testing sets ...')
samples = colnames(Y)
nsams   = length(samples)
tnsams  = sample(samples, as.integer(nsams * tfrac))
ttsams  = setdiff(samples, tnsams)

trainY  = Y [, tnsams]
trainPt = Pt[, tnsams]
testY   = Y [, ttsams]
testPt  = Pt[, ttsams]

# Normlization funcs
L2 = function(mat) {
	t(t(mat)/sqrt(colSums(mat^2)))
}

logger('Normalizing data ...')
D            = L2(D)
trainYrmean  = rowMeans(trainY)
trainY       = t(scale(t(trainY), scale = FALSE))
trainY       = L2(trainY)
trainPtrmean = rowMeans(trainPt)
trainP       = scale(t(trainPt), scale = FALSE)

testY = testY - trainYrmean
testY = L2(testY)
testP = t(testPt) - trainPtrmean

if (svdP > 0) {
	logger('Using SVD to approximate X ...')
	idx   = 1:svdP
	if (is.installed('corpcor')) {
		library(corpcor)
		P.svd = fast.svd(t(trainP))
	} else {
		P.svd = svd(t(trainP))
	}
	diagP = diag(P.svd$d)[idx, idx]
	U     = as.matrix(P.svd$u[, idx])
	V     = as.matrix(P.svd$v[, idx])
	XR    = t(diagP %*% t(V))
	XL    = t(trainY) %*% D
	X     = kronecker(XR, XL)
} else {
	logger('Calculating X ...')
	X = kronecker(trainP, t(trainY) %*% D)
}
vecYtY = as.matrix(as.vector(t(trainY) %*% trainY))

# train 
alpha = 0
if (nthread > 1) {
	library(doMC)
	registerDoMC(cores = nthread)
}
logger('Doing cv.glmnet ...')
cvfit = cv.glmnet(X, vecYtY, alpha=alpha, nfolds=nfold, lower = -1, upper = 1, parallel = nthread > 1)
if (method == 'glmnet') {
	logger('Calculating fit using glmnet ...')
	fit  = glmnet(X, vecYtY, alpha=alpha, lower = -1, upper = 1)
	vecW = predict(fit, X, s=cvfit$lambda.min, type="coefficients")
} else {
	logger('Calculating fit using admm ...')
	if (nthread > 1) {
		vecW = admm_enet(x, y)$penalty(exp(-2), alpha = alpha)$parallel$fit()
	} else {
		vecW = admm_enet(x, y)$penalty(exp(-2), alpha = alpha)$fit()
	}
}

logger('Calculating W ...')
if (svdP > 0) {
	hatW = matrix(vecW[-1], ncol = svdP)
	W    = hatW %*% t(U)
} else {
	W = matrix(as.vector(vecW[-1]), nrow = ncol(D))
}

rownames(W) = colnames(D)
colnames(W) = colnames(trainP)
write.table(round(W, 4), outW, quote = F, sep = "\t")

if (predY) {
	logger('Predicting Y ...')
	predY = D %*% W %*% t(testP)
	predyfile = file.path(outdir, 'predY.txt')
	write.table(predY, predyfile, sep="\t", quote=F)

	spcors  = NULL
	for (tts in ttsams) {
		spcor  = cor.test(predY[, tts], testY[, tts], method = 'spearman')
		#spcors = rbind(spcors, sample = tts, cor = spcor$rho, p = spcor$p.value)
		spcors = rbind(spcors, c(tts, spcor$estimate, spcor$p.value))
	}
	colnames(spcors) = c('sample', 'cor', 'p')
	predycorfile = file.path(outdir, 'predYcor.txt')
	write.table(spcors, predycorfile, sep="\t", quote=F, row.names = F)

	predycmsfile = file.path(outdir, 'predYcorMeanSd.txt')
	cormeansd = matrix(c(mean(as.numeric(spcors[, 'cor'])), sd(as.numeric(spcors[, 'cor']))), ncol = 1)
	rownames(cormeansd) = c('mean', 'sd')
	write.table(cormeansd, predycmsfile, sep="\t", quote=F, col.names = F)
}

if (WPt) {
	logger('Calculating WPt ...')
	outfile = file.path(outdir, 'WPt.txt')
	write.table(W %*% as.matrix(Pt), outfile, sep="\t", quote=F)
}

if (WtDtY) {
	logger('Calculating WtDtY ...')
	outfile = file.path(outdir, 'WtDtY.txt')
	write.table(t(W) %*% t(as.matrix(D)) %*% as.matrix(Y), outfile, sep="\t", quote=F)
}

