library(RUVcorr)

seed     = {{i.seed | R}}
outfile  = {{o.outfile | quote}}
nsamples = {{args.nsamples | R}}
ngenes   = {{args.ngenes | R}}
slabel   = {{args.slabel | quote}}
glabel   = {{args.glabel | quote}}
params   = {{args.params | R}}

set.seed(seed)
update.list = function (list1, list2, recursive = F) {
	names2 = names(list2)
	for (name in names2) {
		if (is.list(list1[[name]]) && is.list(list2[[name]]) && recursive) {
			list1[[name]] = update.list(list1[[name]], list2[[name]], recursive = recursive)
		} else {
			list1[[name]] = list2[[name]]
		}
	}
	return (list1)
}
params   = update.list(params, list(
	k             = 10,
	size.alpha    = 2,
	corr.strength = 3,
	g             = NULL,
	Sigma.eps     = 1,
	# An integer setting the number of negative controls.
	nc            = ngenes/2,
	# An integer setting the number of strongly expressed genes.
	ne            = ngenes/2,
	intercept     = TRUE,
	check         = TRUE
))
params$n = ngenes
params$m = nsamples
exprs    = do.call(simulateGEdata, params)$Y
colnames(exprs) = paste0(glabel, 1:ngenes)
rownames(exprs) = paste0(slabel, 1:nsamples)

write.table(t(exprs), outfile, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


