{{rimport}}('__init__.r')

infile   = {{i.infile | quote}}
survfile = {{i.survfile | quote}}
outdir   = {{o.outdir | quote}}
covfile  = {{args.covfile | R}}
rnames   = {{args.inopts.get('rnames', True) | R}}
methods  = {{args.methods | R}}

survret = read.table(infile, header = T, row.names = NULL, sep = "\t", check.names = F)
survret = unique(survret[, c('var', 'groups')])

survdata = read.table(survfile, header = T, row.names = if (rnames) 1 else NULL, sep = "\t", check.names = F)

if (!is.null(covfile)) {
	covdata = read.table(covfile,  header = T, row.names = 1, sep = "\t", check.names = F)
	covs    = colnames(covdata)
}

for (i in 1:nrow(survret)) {

	var     = as.character(survret[i, 'var'])
	logger('Handling variable', var)
	outfile = file.path(outdir, paste0(gsub("[^0-9a-zA-Z._-]+", "_", var), ".survstats.txt"))
	f = file(outfile)
	lines = c()
	lines = c(lines, paste('# Variable:', var))

	if (!is.null(covfile)) {
		lines = c(lines, paste('  - Covariates:', paste(covs, collapse = ', ')))
	}

	groups = as.integer(unlist(strsplit(as.character(survret[i, 'groups']), ",", fixed = T)))

	vardata = survdata[, var, drop = F]
	if (!is.null(covfile)) {
		vardata = cbind.fill(vardata, covdata)
	}

	vardata = vardata[complete.cases(vardata),, drop = F]
	vardata = vardata[order(vardata[, var]),, drop = F]

	lines = c(lines, '')
	lines = c(lines, paste('# Basic statistics'))
	if (!is.null(covfile)) {
		lines = c(lines, paste('## Main variable:', var))
	}
	lines = c(lines, paste('  - min   :', min(vardata[, var])))
	lines = c(lines, paste('  - max   :', max(vardata[, var])))
	lines = c(lines, paste('  - mean  :', mean(vardata[, var])))
	lines = c(lines, paste('  - median:', median(vardata[, var])))

	if (!is.null(covfile)) {
		for (v in covs) {
			lines = c(lines, '')
			lines = c(lines, paste('## Covariate:', v))
			
			lines = c(lines, paste('  - min   :', min(vardata[, v])))
			lines = c(lines, paste('  - max   :', max(vardata[, v])))
			lines = c(lines, paste('  - mean  :', mean(vardata[, v])))
			lines = c(lines, paste('  - median:', median(vardata[, v])))
		}
	}

	groupdata = list()
	prevgroup = 0
	lines = c(lines, '')
	lines = c(lines, '# Groups:')
	for (j in 1:length(groups)) {
		groupdata[[j]] = vardata[(prevgroup + 1):(prevgroup + groups[j]),,drop = F]
		prevgroup = prevgroup + groups[j]
		lines = c(lines, paste0('  - ', j, '. size: ', groups[j], ', cutting point (min): ', min(groupdata[[j]][, var])))
	}

	if (length(groupdata) >= 2) {
		vars  = var
		if (!is.null(covfile)) {
			vars = c(var, covs)
		}

		for (m in 1:length(groupdata)) {
			for (n in 1:length(groupdata)) {
				if (n <= m) next

				for (v in vars) {
					lines = c(lines, '')
					lines = c(lines, paste0('# Variable "', v, '" (Group', m, ' vs Group', n, ')'))

					for (met in methods) {
						if (met == 't') {
							lines = c(lines, '## T-test')
							tryCatch({
								ret = t.test(groupdata[[m]][, v], groupdata[[n]][, v])
								lines = c(lines, paste('  - P Value:', ret$p.value))
								lines = c(lines, paste('  - ConfInt:', paste(ret$conf.int, collapse = ', ')))
								lines = c(lines, paste0('  - Mean  ', c(m, n), ': ', ret$estimate))
							}, error = function(e) {
								lines = c(lines, paste('  - Not Available:', e))
							})
						} else if (met == 'wilcox') {
							lines = c(lines, '## Wilcox-test')
							ret = wilcox.test(groupdata[[m]][, v], groupdata[[n]][, v])
							lines = c(lines, paste('  - P Value:', ret$p.value))
						} else if (startsWith(met, 'chisq')) {
							if (met == 'chisq') {
								chisqn = 10
							} else {
								chisqn = as.integer(substring(met, 7))
							}

							lvls1 = levels(as.factor(groupdata[[m]][, v]))
							lvls2 = levels(as.factor(groupdata[[n]][, v]))
							if (length(lvls1) > chisqn || length(lvls2) > chisqn) next

							alllvls = union(lvls1, lvls2)
							# contingency table
							t1 = as.data.frame(table(groupdata[[m]][, v]))
							t2 = as.data.frame(table(groupdata[[n]][, v]))
							rownames(t1) = as.character(t1$Var1)
							rownames(t2) = as.character(t2$Var1)
							t1 = t1[, -1, drop=F]
							t2 = t2[, -1, drop=F]
							colnames(t1) = paste0("Group", m)
							colnames(t2) = paste0("Group", n)
							ctable = cbind.fill(t1, t2)
							ctable[is.na(ctable)] = 0
							ctable = ctable[order(rownames(ctable)), , drop = F]
							lines = c(lines, '## Chisquare-test')
							tryCatch({
								ret = chisq.test(ctable)
								lines = c(lines, paste('  - P Value:', ret$p.value))
								lines = c(lines, paste('  - Chi2   :', as.character(ret$statistic)))
								lines = c(lines, '  - Observed (Contingency table):')
								lines = c(lines, paste(c('     ', colnames(ret$observed)), collapse = "\t"))
								retrnames = rownames(ret$observed)
								for (o in 1:nrow(ret$observed)) {
									lines = c(lines, paste(c(paste0('    ', retrnames[o]), ret$observed[o, ]), collapse = "\t"))
								}
								lines = c(lines, '  - Expected:')
								lines = c(lines, paste(c('     ', colnames(ret$expected)), collapse = "\t"))
								retrnames = rownames(ret$expected)
								for (o in 1:nrow(ret$expected)) {
									lines = c(lines, paste(c(paste0('    ', retrnames[o]), ret$expected[o, ]), collapse = "\t"))
								}
							}, function(e) {
								lines = c(lines, paste('  - Not Available:', e))
							})
						}
					}

				}
			}
		}
	}
	writeLines(lines, f)
	close(f)
}



