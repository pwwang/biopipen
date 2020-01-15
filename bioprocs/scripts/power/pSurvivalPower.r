{{rimport}}('plot.r', '__init__.r')

set.seed(8525)

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
outdir  = {{o.outdir | R}}
intype  = {{args.intype | R}}
rnames  = {{args.rnames | lambda x: 1 if x else 'NULL' | R}}
ngroup  = as.integer({{args.ngroup | R}})
alphas  = {{args.alphas | lambda x: x if isinstance(x, list) else [float(a.strip()) for a in x.split(',')] if isinstance(x, str) else [x] | R}}
betas   = {{args.betas | lambda x: x if isinstance(x, list) else [float(a.strip()) for a in x.split(',')] if isinstance(x, str) else [x] | R}}
ngroup1 = {{i.ngroup1 | : int(a) if a else 0}}
ngroup2 = {{i.ngroup2 | : int(a) if a else 0}}
ngroup3 = {{i.ngroup3 | : int(a) if a else 0}}
ngroup4 = {{i.ngroup4 | : int(a) if a else 0}}

ncut = function(dat, n, prefix = 'q') {
	diffdata = max(dat) - min(dat)
	if (diffdata == 0) diffdata = 1e-2
	diffdata = diffdata * 1e-8
	fakedata = dat + runif(length(dat), diffdata/10, diffdata)
	q = quantile(fakedata, probs = seq(0, 1, 1/n), include.lowest = T)
	return(cut(fakedata, q, labels = paste0(prefix, 1:n), include.lowest = T))
}

assignGroup = function(dat, ngroups, prefix = 'q') {
	r = rank(dat, ties.method = 'first')
	b = 0
	for (i in 1:length(ngroups)) {
		dat[r > b & r <= b + ngroups[i]] = paste0(prefix, i)
		b = b + ngroups[i]
	}
	return(dat)
}

data = list()
if (intype == 'detail' || intype == 'detailed') {
	detailed  = read.table(infile, header = T, row.names = rnames, sep = '\t', check.names = F)
	censoring = as.integer(levels(factor(detailed[, 2, drop = T])))
	alive = min(censoring)
	cnames = colnames(detailed)
	for (x in 3:length(cnames)) {
		var    = cnames[x]
		groups = levels(factor(detailed[, var, drop = T]))
		lengroup = length(groups)
		if (lengroup < 2) {
			logger('Variable', var, 'skipped, as only one ground found.', level = 'WARN')
			next
		} else if (lengroup > ngroup) {
			if (ngroup1 == 0) {
				detailed[, var] = ncut(detailed[, var], ngroup)
				lengroup = ngroup
				groups   = levels(factor(detailed[, var, drop = T]))
			} else {
				ngroups = c()
				if (ngroup1 > 0) ngroups = c(ngroups, ngroup1)
				if (ngroup2 > 0) ngroups = c(ngroups, ngroup2)
				if (ngroup3 > 0) ngroups = c(ngroups, ngroup3)
				if (ngroup4 > 0) ngroups = c(ngroups, ngroup4)
				detailed[, var] = assignGroup(detailed[, var], ngroups)
				lengroup = length(ngroups)
				groups   = levels(factor(detailed[, var, drop = T]))
			}
		}

		for (i in 1:(lengroup-1)) {
			for (j in (i+1):lengroup) {
				group1 = groups[i]
				group2 = groups[j]
				variable = paste0(var, '(', group1, '-', group2, ')')
				data1  = detailed[which(detailed[[var]] == group1), , drop = F]
				data2  = detailed[which(detailed[[var]] == group2), , drop = F]
				nrow1  = nrow(data1)
				nrow2  = nrow(data2)
				data[[variable]] = list(
					survrate1 = nrow(data1[which(data1[, 2] == alive), var, drop = F]) / nrow1,
					survrate2 = nrow(data2[which(data2[, 2] == alive), var, drop = F]) / nrow2,
					ssratio   = nrow1 / nrow2
				)
			}
		}
	}
} else {
	rdata = read.table(infile, header = T, row.names = rnames, sep = '\t', check.names = F)
	if (!rnames) {
		rownames(rdata) = paste0('Var', 1:nrow(rdata))
	}
	for (rname in rownames(rdata)) {
		data[[rname]] = list(
			survrate1 = rdata[rname, 1],
			survrate2 = rdata[rname, 2],
			ssratio   = rdata[rname, 3]
		)
	}
}

samsize = function(survrate1, survrate2, ssratio, a, b) {
	za     = abs(qnorm(1-a/2))
	zb     = abs(qnorm(1-b))
	sqrtA  = za + zb
	A      = sqrtA * sqrtA
	hr     = (1-survrate1)/(1-survrate2 + 1e-3)
	loghr  = log(hr) # ln
	# B    = pi1*pi2*(logHR)^2
	# B    = (ssratio/(ssratio + 1))*(1/(ssratio + 1)) * loghr * loghr
	B      = ssratio * loghr * loghr / (ssratio + 1) / (ssratio + 1)
	events = A/B

	return(list(
		ssize1 = as.integer(ssratio * events / (ssratio + 1)),
		ssize2 = as.integer(events / (ssratio + 1)),
		total  = as.integer(events),
		# hazard ratio
		hr     = hr
	))
}

outs = matrix(ncol = 10, nrow = 0)
colnames(outs) = c('Variable', 'Alpha', 'Beta', 'Survrate1', 'Survrate2', 'SSRatio', 'HRatio', 'SSize1', 'SSize2', 'Total')

alpha_plots = alphas
beta_plots  = betas
{% if args.plot %}
alpha_plots = c(exp(seq(-7, -2, by = .5)), alphas)
{% endif %}

for (var in names(data)) {
	plot_hr   = NA
	plot_data = matrix(ncol = 3, nrow = 0)
	colnames(plot_data) = c('mlogp', 'ssize', 'power')
	for (a in alpha_plots) {
		for (b in beta_plots) {
			ret = samsize(data[[var]]$survrate1, data[[var]]$survrate2, data[[var]]$ssratio, a, b)
			plot_hr   = ret$hr
			plot_data = rbind(plot_data, c(
				mlogp = -log2(a),
				ssize = ret$total,
				power = 1-b
			))
			if (a %in% alphas && b %in% betas) {
				outs = rbind(outs, c(
					Variable  = var,
					Alpha     = a,
					Beta      = b,
					Survrate1 = round(data[[var]]$survrate1, 3),
					Survrate2 = round(data[[var]]$survrate2, 3),
					SSRatio   = round(data[[var]]$ssratio, 3),
					HRatio    = round(ret$hr, 3),
					SSize1    = ret$ssize1,
					SSize2    = ret$ssize2,
					Total     = ret$total
				))
			}

		}
	}
	{% if args.plot %}
	plotfile = file.path(outdir, paste0(gsub("[[:punct:]]", "_", var), '.png'))
	vldata   = data.frame(x = -log2(alphas), y=-Inf, pval = alphas, label = paste('p = ', alphas))
	ggs = list(
		xlab = list("-log2(p-value)"),
		ylab = list("Sample size"),
		labs = list(color = 'Power'),
		geom_smooth = list(aes(color = factor(power)), method = lm, se = F),
		geom_vline  = list(
			aes(xintercept = x),
			color = '#555555',
			vldata,
			linetype = 'dotdash',
			inherit.aes = FALSE
		),
		geom_text = list(
			aes(x = x, y = y, label = label),
			vldata,
			color = '#555555',
			angle = 90,
			vjust = 1.5,
			hjust = -.2,
			inherit.aes = FALSE
		),
		annotate = list(
			"text",
			label = paste('Hazrd ratio =', round(plot_hr, 2)),
			x = Inf,
			y = -Inf,
			color = '#555555',
			vjust = -.8,
			hjust = 1.1
		)
	)
	plot.points(
		as.data.frame(plot_data),
		plotfile,
		params = list(mapping = aes(color = factor(power))),
		ggs = ggs
	)
	{% endif %}
}

write.table(outs, outfile, sep = '\t', row.names = F, col.names = T, quote = F)
