
normW   = {{args.normW | R}}
normE   = {{args.normE | R}}

W = read.table ({{i.Wfile | quote}}, sep="\t", header=T, row.names=1, check.names=F, strip.white=T)
W = as.matrix(W)

if (normW) {
	library(NetPreProc)
	W = abs(W)
	W = Laplacian.norm(W)
}

E = read.table ({{i.Efile | quote}}, header=F, row.names=1, check.names=F, strip.white=T)
E = as.matrix(E)
E = E[colnames(W), ]

RWR = function (W, e, c = {{args.c}}, eps = {{args.eps}}, tmax={{args.niter}}) {
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

write.table (round(r$r, 3), {{o.outfile | quote}}, quote=F, col.names=F, row.names = T, sep="\t")
