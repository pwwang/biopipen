
dat = read.table('{{in.infile}}', header = T, row.names = 1, sep="\t", check.names = F)
sim = cor(t(m), method = '{{args.method}}')
p   = sim
n   = ncol(m)

{% if args.pval %}
{% if args.method | lambda x: x == 'pearson' or x == 'spearman' %}
t = sim / sqrt((1-sim*sim)/(n-2))
p = 2*pt(-abs(t), n-2)

{% elif args.method | lambda x: x == 'kendall' %}
z = 3*sim * sqrt( n*(n-1)/2/(2*n+5) )
p = 2*pnorm(-abs(z))

{% endif %}
{% endif %}

write.table(sim, '{{out.outfile}}', sep="\t", row.names=T, col.names=T, quote=F)
write.table(p  , '{{out.outpval}}', sep="\t", row.names=T, col.names=T, quote=F)