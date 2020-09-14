{{"__init__.R", "plot.R" | rimport}}

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts = {{args.inopts | R}}
devpars = {{args.devpars | R}}
ggs = {{args.ggs | R}}
params = {{args.params | R}}
corr = {{args.corr | R}}
stacked = {{args.stacked | R}}
tsform = {{args.tsform | ?.startswith: "function" | !R}}
col_x = {{args.x | R}}
line = {{args.line | R}}

data = read.table.inopts(infile, inopts)
cnames = colnames(data)

if (stacked) {
	library(reshape2)
	if (is.null(cnames)) {
		cnames = c("Item", "Value", "Group")
		colnames(data) = cnames
	}
	data = dcast(data,
				 as.formula(paste(bQuote(cnames[1]), "~", bQuote(cnames[3]))),
				 value.vars = cnames[2])
	rownames(data) = data[, 1]
	data = data[, -1, drop=FALSE]
}

# transformation if needed
if (is.function(tsform)) {
	data[, 1] = tsform(data[, 1])
	data[, 2] = tsform(data[, 2])
}

# make sure expression is parsable
tryCatch({
	expr = rlang::parse_expr(line)
}, error = function(e) {
	stop(paste("Wrong format of line expression:", str(e)))
})

# find groups
slope = 0
intercept = 0
# get slope and intercept first
x = 0
eval(expr)
tryCatch({
	intercept = y
}, error = function(e) {
	stop(paste("Line expression should be in format: y=2*x+1", str(e)))
})
x = 1
eval(expr)
slope = y - intercept
x = data[, col_x]
eval(expr)
data$Group = data[, 2] >= y
data[data$Group, "Group"] = paste(cnames[2], "outperformed")
data[data$Group == FALSE, "Group"] = paste(cnames[1], "outperformed")

# get correlation
corr_value = NULL
if (corr %in% c("pearson", "spearman", "kendall")) {
	corr_value = cor(data[, 1], data[, 2], method = corr)
}

# plot diagonal line and legend
params = c(list(ggplot2::aes_string(color = "Group")), params)
ggs = c(list(geom_abline = list(intercept = intercept, slope = slope)), ggs)

# add correlation text
ggs = c(list(geom_text = list(x = Inf, y = -Inf,
							  hjust = 1, vjust = 0,
							  label = round(corr_value, 3))), ggs)

plot.scatter(data, outfile,
			 x = col_x,
			 y = ifelse(is.numeric(col_x), 3-col_x, cnames[1:2][cnames[1:2]!=col_x]),
			 params = params,
			 ggs = ggs,
			 devpars = devpars)
