
dat = read.table('{{i.infile}}', header = T, row.names = 1, sep="\t", check.names = F)
sim = cor(t(dat), method = '{{args.method}}')
p   = sim
n   = ncol(dat)

{% if args.pval %}
{% if args.method == 'pearson' or args.method == 'spearman' %}
t = sim / sqrt((1-sim*sim)/(n-2))
p = 2*pt(-abs(t), n-2)

{% elif args.method == 'kendall' %}
z = 3*sim * sqrt( n*(n-1)/2/(2*n+5) )
p = 2*pnorm(-abs(z))

{% endif %}
{% endif %}

write.table(format(sim, digits=4, scientific=F), '{{o.outfile}}', sep="\t", row.names=T, col.names=T, quote=F)
write.table(format(p, digits=3, scientific=T)  , '{{o.outpval}}', sep="\t", row.names=T, col.names=T, quote=F)