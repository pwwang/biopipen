
{% case args.infmt %}
# full
{% when 'full' %}
mat = read.table({{i.infile | quote}}, row.names = 1, header = T, check.names = F, sep = "\t")

# upper
{% when 'upper' %}
mat = read.table({{i.infile | quote}}, row.names = 1, header = T, check.names = F, sep = "\t", fill = T)
mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

# lower
{% when 'lower' %}
mat = read.table({{i.infile | quote}}, row.names = 1, header = T, check.names = F, sep = "\t", fill = T)
mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]

# pair
{% when 'pair' %}
mat0 <- read.table({{i.infile | quote}}, header=F, row.names=NULL, stringsAsFactor=F, check.names=F)
name1 = as.character (mat0[, 1])
name2 = as.character (mat0[, 2])
names = unique(union(name1, name2))
nlen  = length(names)
mat   = matrix(NA, ncol = nlen, nrow = nlen)
rownames(mat) = names
colnames(mat) = names
diag(mat)     = 0
mat[as.matrix(mat0[,1:2,drop=F])] = mat0[,3]

ut = upper.tri(mat)
lt = lower.tri(mat)

utdata = mat[ut]
ltdata = mat[lt]

uttdata = t(mat)[ut]
lttdata = t(mat)[lt]

utdata[is.na(utdata)] = uttdata[!is.na(uttdata)]
ltdata[is.na(ltdata)] = lttdata[!is.na(lttdata)]

mat[ut] = utdata
mat[lt] = ltdata

{% endcase %}

coords = cmdscale(mat, k={{args.k}})
coords = round(coords, 3)
write.table (coords, {{o.outfile | quote}}, col.names=T, row.names=T, quote=F, sep="\t", na="")
