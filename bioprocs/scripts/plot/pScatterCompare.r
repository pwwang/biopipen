{{rimport}}('plot.r', '__init__.r')

library(ggplot2) # make aes available

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
cnames  = {{args.cnames | R }}  # required to be TRUE
rnames  = {{args.rnames | R }}
devpars = {{args.devpars | R}}
x       = {{args.x | repr}}
y       = {{args.y | repr}}

data = read.table(infile, header = cnames, row.names = if (rnames) 1 else NULL, sep = '\t', check.names = F)
if (!is.na(as.numeric(x))) x = as.integer(x)
if (!is.na(as.numeric(y))) y = as.integer(y)
cnames = colnames(data)
m1 = if (is.numeric(x)) cnames[x] else x
m2 = if (is.numeric(y)) cnames[y] else y

data$.group = apply(data[,c(x, y)], 1, function(row) if(row[1]>row[2]) m1 else if (row[1]<row[2]) m2 else 'None')

eval(parse(text = {{args.helper | repr}}))
params  = {{args.params | R}}
ggs     = {{args.ggs | R}}

if ("mapping" %in% names(params)) {
	params$mapping = update.list(aes(color = .group), params$mapping)
} else if ('' %in% names(params)) {
	params[[1]] = update.list(aes(color = .group), params[[1]])
} else {
	params$mapping = aes(color = .group)
}

# plot diagonal line and legend
ggs = c(
	list(
		geom_abline = list(intercept = 0, slop = 1),
		scale_color_discrete = list(name = 'Winner')
	),
	ggs
)

plot.scatter(data, outfile, x, y, params, ggs, devpars = devpars)
