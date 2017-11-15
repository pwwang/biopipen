library(ggplot2)

mat = read.table({{in.infile | quote}}, header = T, check.names = F, sep = "\t", row.names = as.integer({{args.rownames | R}}))
if (ncol(mat) < 2) {
	stop('Cannot plot scattercompare with less than columns.')
}

cnames = colnames(mat)
maxy   = max(mat[, 1])
miny   = min(mat[, 1])
maxx   = max(mat[, 2])
minx   = min(mat[, 2])
maxa   = max(maxx, maxy) * 1.1
mina   = min(minx, miny) * 1.1

mat$better = apply(mat, 1, function(row) if (row[1] >= row[2]) paste(cnames[1], "better") else paste(cnames[2], "better"))
do.call(png, c(list(filename={{out.outfile | quote}}), {{args.devpars | Rlist}}))
p      = ggplot(mat, aes_string(y = cnames[1], x = cnames[2])) + 
		 coord_fixed() + 
		 xlim(mina, maxa) + 
		 ylim(mina, maxa) 

{% if args.diag %}
p = p + geom_abline(intercept = 0, slope = 1, color = "#AAAAAA") 
{% endif %}


{% if args.corr | lambda x: x in ['pearson', 'spearman', 'kendall'] %}
correlation = cor(mat[,1,drop=T], mat[,2,drop=T], method={{args.corr | quote}})
correlation = format(round(correlation, 3), nsmall=3)
p = p + geom_text(x = maxa * 1.1, y = mina, label = paste("cor", correlation, sep = " = "), size = 3, hjust = 2)
{% endif %}

ggs      = {{args.ggs | Rlist}}
haspoint = F
if (length(ggs) > 0) {
	for (i in 1:length(ggs)) {
		g = ggs[[i]]
		if (!'geom' %in% names(g)) next
		k = class(g$geom)[1]
		if (k == 'GeomPoint') {
			haspoint = T
			break
		}
	}
}
if (!haspoint) p = p + geom_point(aes(color = factor(better))) + theme(legend.title = element_blank())
if (length(ggs) > 0) {
	for (i in 1:length(ggs)) {
		p = p + ggs[[i]]
	}
}

{% if args.regr %}
lm_eqn = function(y, x){
	m  = lm(y ~ x);
	eq = substitute(italic(y) == b ~ italic(x) + a*","~~italic(r)^2~"="~r2, 
		 list(a = format(coef(m)[1], digits = 2), 
			  b = format(coef(m)[2], digits = 2), 
			 r2 = format(summary(m)$r.squared, digits = 3)))
	as.character(as.expression(eq));
}
p = p + geom_text(x = mina * 1.1, y = maxa, label = lm_eqn(mat[,2,drop=T], mat[,1,drop=T]), parse = TRUE, size = 3, hjust = 0)
p = p + geom_smooth(method = 'lm', se=F, fullrange=T)
{% endif %}

print(p)
dev.off()
