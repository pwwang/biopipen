options(stringsAsFactors = FALSE)

inD  = {% if i.D.endswith('.gz') %}gzfile({{i.D  | quote}}){% else %}{{i.D  | quote}}{% endif %}
inPt = {% if i.Pt.endswith('.gz') %}gzfile({{i.Pt | quote}}){% else %}{{i.Pt | quote}}{% endif %}
inY  = {% if i.Y.endswith('.gz') %}gzfile({{i.Y  | quote}}){% else %}{{i.Y  | quote}}{% endif %}
D  = read.table(inD,  sep = "\t", check.names = F, row.names = 1, header = T)
Pt = read.table(inPt, sep = "\t", check.names = F, row.names = 1, header = T)
Y  = read.table(inY,  sep = "\t", check.names = F, row.names = 1, header = T)

set.seed({{args.seed | lambda x: 'NULL' if x<0 else x}})
samples = colnames(Y)
nsams   = length(samples)
tnsams  = sample(samples, as.integer(nsams * {{args.tfrac}}))
ttsams  = setdiff(samples, tnsams)

trainY  = Y [, tnsams]
trainPt = Pt[, tnsams]
testY   = Y [, ttsams]
testPt  = Pt[, ttsams]

# Normlization funcs
L2 = function(mat) {
	t(t(mat)/sqrt(colSums(mat^2)))
}

meanCenter = function(mat) {
	rmean = rowMeans(mat)
	mcmat = mat - rmean
	L2(mcmat)
}

{% if args.normD %}
D = {{args.normD}}(D)
{% endif %}

{% if args.normPt %}
trainPt = {{args.normPt}}(trainPt)
testPt  = {{args.normPt}}(testPt)
{% endif %}

{% if args.normY %}
trainY = {{args.normY}}(trainY)
testY  = {{args.normY}}(testY)
{% endif %}

{% if args.svdP %}
idx    = 1:{{args.svdP}}
svdP   = svd(trainPt)
diagP  = diag(svdP$d)[idx, idx]
U      = as.matrix(svdP$u[, idx])
V      = as.matrix(svdP$v[, idx])
XR     = t(diagP %*% t(V))
XL     = t(trainY) %*% D
X      = kronecker(XR, XL)
{% else %}
X      = kronecker(t(trainPt), t(trainY) %*% D)
{% endif %}
vecYtY = as.matrix(as.vector(t(trainY) %*% trainY))

# train 
require(methods)

require(glmnet)
alpha = 0
{% if args.parallel %}
require(doMC)
registerDoMC(cores={{args.nfolds}})
{% endif %}
cvfit = cv.glmnet(X, vecYtY, alpha=alpha, nfolds={{args.nfolds}}, lower = -1, upper = 1, parallel = {{args.parallel | R}})

{% if args.method == 'glmnet' %}
fit   = glmnet(X, vecYtY, alpha=alpha, lower = -1, upper = 1)
# learn
vecW  = predict(fit, X, s=cvfit$lambda.min, type="coefficients")
{% elif args.method == 'admm' %}
{%    if args.parallel %}
vecW  = admm_enet(x, y)$penalty(exp(-2), alpha = alpha)$parallel$fit()
{% 	  else %}
vecW  = admm_enet(x, y)$penalty(exp(-2), alpha = alpha)$fit()
{% 	  endif %}
{% endif %}

{% if args.svdP %}
hatW  = matrix(vecW[-1], ncol = {{args.svdP}})
W     = hatW %*% t(U) 
{% else %}
W     = matrix(as.vector(vecW[-1]), nrow = ncol(D))
{% endif %}
rownames(W) = colnames(D)
colnames(W) = rownames(trainPt)

write.table(round(W, 4), {{o.W | quote}}, quote = F, sep = "\t")

{% if args.predY %}
predY = D %*% W %*% testPt
predyfile = file.path({{o.outdir | quote}}, 'predY.txt')
write.table(predY, predyfile, sep="\t", quote=F)

spcors  = NULL
for (tts in ttsams) {
	spcor  = cor.test(predY[, tts], testY[, tts], method = 'spearman')
	#spcors = rbind(spcors, sample = tts, cor = spcor$rho, p = spcor$p.value)
	spcors = rbind(spcors, c(tts, spcor$estimate, spcor$p.value))
}
colnames(spcors) = c('sample', 'cor', 'p')
predycorfile = file.path({{o.outdir | quote}}, 'predYcor.txt')
write.table(spcors, predycorfile, sep="\t", quote=F, row.names = F)

predycmsfile = file.path({{o.outdir | quote}}, 'predYcorMeanSd.txt')
cormeansd = matrix(c(mean(as.numeric(spcors[, 'cor'])), sd(as.numeric(spcors[, 'cor']))), ncol = 1)
rownames(cormeansd) = c('mean', 'sd')
write.table(cormeansd, predycmsfile, sep="\t", quote=F, col.names = F)
{% endif %}

{% if args.WPt %}
outfile = file.path({{o.outdir | quote}}, 'WPt.txt')
write.table(W %*% as.matrix(Pt), outfile, sep="\t", quote=F)
{% endif %}

{% if args.WtDtY %}
outfile = file.path({{o.outdir | quote}}, 'WtDtY.txt')
write.table(t(W) %*% t(as.matrix(D)) %*% as.matrix(Y), outfile, sep="\t", quote=F)
{% endif %}

