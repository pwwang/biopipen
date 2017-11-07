library("metap")

indata = read.table({{in.infile | quote}}, header = {{args.header | R}}, row.names = NULL, check.names = F)
rnames = unique(as.vector(indata[, 1]))
pcol   = {{args.pcol}}
if (pcol < 0) pcol = ncol(indata) + pcol + 1

ret = NULL
for (rname in rnames) {
	pvals  = as.vector(indata[which(indata[,1] == rname), pcol])
	if(length(pvals) == 1) {
		{% if args.method | lambda x: x == 'logitp' %}
		mp1 = list(t = 0, df = 0, p = pvals[1], validp = pvals[1])
		{% elif args.method | lambda x: x == 'sumz' %}
		mp1 = list(z = 0, p = pvals[1], validp = pvals[1])
		{% elif args.method | lambda x: x == 'votep' %}
		mp1 = list(p = pvals[1], pos = 0, neg = 0, alpha = 0, validp = pvals[1])
		{% elif args.method | lambda x: x == 'sump' %}
		mp1 = list(p = pvals[1], conservativep = pvals[1], validp = pvals[1])
		{% elif args.method | lambda x: x == 'meanp' %}
		mp1 = list(z = 0, p = pvals[1], validp = pvals[1])
		{% elif args.method | lambda x: x == 'wilkinsonp' %}
		mp1 = list(p = pvals[1], pr = 0, r = 0, critp = pvals[1], alpha = 0, validp = pvals[1])
		{% endif %}
	} else {
		mp1 = {{args.method}}(pvals)
	}
	mp2    = unlist(mp1)
	mnames = names(mp1)
	mnames = if ("weights" %in% mnames) mnames[1:length(mnames)-2] else mnames[1:length(mnames)-1]
	mmp    = matrix(mp2[1:length(mnames)], nrow = 1)
	rownames(mmp) = rname
	colnames(mmp) = mnames
	mmp    = cbind(mmp, validp = paste(mp1$validp, collapse=';'))
	if (is.null(ret)) {
		ret = mmp
	} else {

		ret = rbind(ret, mmp)
	}
}
ret = {% if args.poutonly %}ret[order(as.numeric(ret[, 'p'])), 'p', drop=F]{% else %}ret[order(as.numeric(ret[, 'p'])),,drop=F]{% endif %}
write.table(ret, {{out.outfile | quote}}, sep="\t", quote=F, col.names={{args.outheader | R}})
