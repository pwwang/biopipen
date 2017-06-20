from pyppl import proc

"""
@name:
	pRWR
@description:
	Do random walk with restart (RWR)
@input:
	`Wfile:file`: The adjecent matrix
	`Efile:file`: The start vector
@output:
	`outfile:file`: The output of final probabilities
@args:
	`c`:       The restart probability
	`eps`:     The convergent cutoff || R(i+1) - R(i) ||
	`tmax`:    Max iterations to stop
	`Wformat`: The format of Wfile, rds or mat/txt, default: rds
	`Eformat`: The format of Efile, rds or mat/txt, default: rds
	`Rformat`: The format of the output file, rds or mat/txt, default: rds
	`normW`:   Weather to normalize W or not, default False. Laplacian normalization is used (more to add).
	`normE`:   Weather to normalize E or not, default False. E will be normalized as: E = E/sum(E)
@requires:
	if normW = True, R package `NetPreProc` is required.
"""
pRWR = proc ()
pRWR.input     = "Wfile:file, Efile:file"
pRWR.output    = "outfile:file:R-{{Efile | fn}}"
pRWR.args      = {'c': 0.1, 'eps': 1e-5, 'tmax': 10000, 'Wformat': 'rds', 'Eformat': 'rds', 'Rformat': 'rds', 'normW': False, 'normE': False}
pRWR.defaultSh = "Rscript"
pRWR.script    = """
normW   = {{proc.args.normW | Rbool}}
normE   = {{proc.args.normE | Rbool}}
Wformat = "{{proc.args.Wformat}}"
Eformat = "{{proc.args.Eformat}}"
Rformat = "{{proc.args.Rformat}}"

if (Wformat == "rds") {
	W = readRDS ("{{Wfile}}")
} else {
	W = read.table ("{{Wfile}}", sep="\\t", header=T, row.names=1, check.names=F, strip.white=T)
	W = as.matrix(W)
}

if (normW) {
	library(NetPreProc)
	W = abs(W)
	W = Laplacian.norm(W)
}

if (Eformat == "rds") {
	E = readRDS ("{{Efile}}")
} else {
	E = read.table ("{{Efile}}", header=F, row.names=1, check.names=F, strip.white=T)
}
E = as.matrix(E)
E = E[colnames(W), ]

RWR = function (W, e, c = {{proc.args.c}}, eps = {{proc.args.eps}}, tmax={{proc.args.tmax}}) {
	r0 = e
	for (i in 1:tmax) {
		r1 = ((1-c) * W) %*% r0 + c * e
		diff = norm (r1 - r0, 'f')
		if (diff < eps) break
		r0 = r1
	}

	return (list(r=r1, eps=diff, tmax=i))
}

r = RWR (W, E)
print (paste("eps: ", r$eps))
print (paste("tmax:", r$tmax))

if (Rformat == 'rds') {
	saveRDS (r$r, "{{outfile}}")
} else {
	write.table (format(r$r, digits=3), "{{outfile}}", quote=F, col.names=F, row.names = T, sep="\\t")
}
"""




