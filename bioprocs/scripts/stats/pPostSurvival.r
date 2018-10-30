{{rimport}}('__init__.r')
library('openxlsx')
options("openxlsx.maxWidth" = 25)
options("openxlsx.minWidth" = 5)

infile   = {{i.infile | quote}}
survfile = {{i.survfile | quote}}
outfile  = {{o.outfile | quote}}
covfile  = {{args.covfile | R}}
rnames   = {{args.inopts.get('rnames', True) | R}}
chi2n    = {{args.chi2n | R}}

survret = read.table(infile, header = T, row.names = NULL, sep = "\t", check.names = F)
survret = unique(survret[, c('var', 'groups')])

survdata = read.table(survfile, header = T, row.names = if (rnames) 1 else NULL, sep = "\t", check.names = F)

covs = c()
if (!is.null(covfile)) {
	covdata = read.table(covfile,  header = T, row.names = 1, sep = "\t", check.names = F)
	covs    = colnames(covdata)
}

lighten = function(color, factor = 1.5){
	col = col2rgb(color)
	col = col*factor
	col[col>255] = 255
	col = rgb(t(col), maxColorValue=255)
	col
}
# 	A	B
# R1: 1 (6.67%, 1.82%)	6 (15.00%, 10.91%)
# R2: 2 (13.33%, 3.64%)	7 (17.50%, 12.73%)
# R3: 3 (20.00%, 5.45%)	8 (20.00%, 14.55%)
# R4: 4 (26.67%, 7.27%)	9 (22.50%, 16.36%)
# R5: 5 (33.33%, 9.09%)	10 (25.00%, 18.18%)
ctable.summary = function(ctable) {
	lines  = paste(colnames(ctable), collapse = "; ")
	rnames = rownames(ctable)
	for (i in 1:nrow(ctable)) {
        vals = sapply(1:length(ctable[i,]), function(x)  sprintf("%d (%.2f%%, %.2f%%)", ctable[i,x], 100*ctable[i,x]/sum(ctable[,x]), 100*ctable[i,x]/sum(ctable)))
        #vals = sapply(ctable[i,], function(x) sprintf("%d (%.2f%%, %.2f%%)", x, 100*x/sum(ctable[,1]), 100*x/sum(ctable)))
		lines = c(lines, paste0(paste0(rnames[i], ": "), paste(vals, collapse="; ")))
	}
	paste(lines, collapse = "\n")
}

wb = createWorkbook()
for (i in 1:nrow(survret)) {

	var     = as.character(survret[i, 'var'])
	vars    = c(var, covs)
	logger('Handling variable', var)
	groups  = as.integer(unlist(strsplit(as.character(survret[i, 'groups']), ",", fixed = T)))
	ngroups = length(groups)
	# nothing to do if less 2 groups
	if (ngroups < 2) next
	ncomb   = ngroups * (ngroups -1) / 2

	vardata = survdata[, var, drop = F]
	# alive
	alvdata = survdata[which(survdata$status == min(survdata$status)), var, drop = F]
	# dead
	deadata = survdata[setdiff(rownames(survdata), rownames(alvdata)), var, drop = F]
	if (!is.null(covfile)) {
		vardata = cbind.fill(vardata, covdata)
		alvdata = cbind.fill(alvdata, covdata)
		deadata = cbind.fill(deadata, covdata)
	}

	vardata = vardata[complete.cases(vardata),, drop = F]
	vardata = vardata[order(vardata[, var]),, drop = F]

	alvdata = alvdata[complete.cases(alvdata),, drop = F]
	alvdata = alvdata[order(alvdata[, var]),, drop = F]

	deadata = deadata[complete.cases(deadata),, drop = F]
	deadata = deadata[order(deadata[, var]),, drop = F]

	out = matrix(, ncol = 3 + ngroups + 3 + ncomb * 3 + 1 + ncomb + 8, nrow = 5 + length(covs) * 5)
	cnamesout = c('All', 'Alive', 'Deceased', paste0('Group_', 1:ngroups))
	# t-test
	for (tn in c('T_P_', 'T_ConfInt_', 'T_Means_')) {
		cnamesout = c(cnamesout, paste0(tn, 'AvsD'))
	}
	for (tn in c('T_P_', 'T_ConfInt_', 'T_Means_')) {
		for (j in 1:ngroups) {
			for (k in 1:ngroups) {
				if (k >= j) next
				cnamesout = c(cnamesout, paste0(tn, j, 'vs', k))
			}
		}
	}
	# wilcox-test
	cnamesout = c(cnamesout, paste0('Wilcox_P_AvsD'))
	for (j in 1:ngroups) {
		for (k in 1:ngroups) {
			if (k >= j) next
			cnamesout = c(cnamesout, paste0('Wilcox_P_', j, 'vs', k))
		}
	}
	# chi2-test
	cnamesout = c(cnamesout, 'Chi2_P_AvsD', 'Chi2_Chi2_AvsD', 'Chi2_Obs_AvsD', 'Chi2_Exp_AvsD', 'Chi2_P_Groups', 'Chi2_Chi2_Groups', 'Chi2_Obs_Groups', 'Chi2_Exp_Groups')
	colnames(out) = cnamesout
	rnamesout = c()
	for (v in c(var, covs)) {
		rnamesout = c(rnamesout, paste0(v, c('_Min', '_Max', '_Mean', '_Median', '_Stdev')))
	}
	rownames(out) = rnamesout

	for (v in vars) {
		out[paste0(v, '_Min'), 'All']    = round(min(vardata[, v]), 2)
		out[paste0(v, '_Max'), 'All']    = round(max(vardata[, v]), 2)
		out[paste0(v, '_Mean'), 'All']   = round(mean(vardata[, v]), 2)
		out[paste0(v, '_Median'), 'All'] = round(median(vardata[, v]), 2)
		out[paste0(v, '_Stdev'), 'All']  = round(sd(vardata[, v]), 2)

		out[paste0(v, '_Min'), 'Alive']    = round(min(alvdata[, v]), 2)
		out[paste0(v, '_Max'), 'Alive']    = round(max(alvdata[, v]), 2)
		out[paste0(v, '_Mean'), 'Alive']   = round(mean(alvdata[, v]), 2)
		out[paste0(v, '_Median'), 'Alive'] = round(median(alvdata[, v]), 2)
		out[paste0(v, '_Stdev'), 'Alive']  = round(sd(alvdata[, v]), 2)

		out[paste0(v, '_Min'), 'Deceased']    = round(min(deadata[, v]), 2)
		out[paste0(v, '_Max'), 'Deceased']    = round(max(deadata[, v]), 2)
		out[paste0(v, '_Mean'), 'Deceased']   = round(mean(deadata[, v]), 2)
		out[paste0(v, '_Median'), 'Deceased'] = round(median(deadata[, v]), 2)
		out[paste0(v, '_Stdev'), 'Deceased']  = round(sd(deadata[, v]), 2)
	}

	prevgroup = 0
	groupdata = list()
	for (j in 1:ngroups) {
		groupdata[[j]] = vardata[(prevgroup + 1):(prevgroup + groups[j]),,drop = F]
		prevgroup = prevgroup + groups[j]
		for (v in vars) {
			out[paste0(v, '_Min'), paste0('Group_', j)]    = round(min(groupdata[[j]][, v]), 2)
			out[paste0(v, '_Max'), paste0('Group_', j)]    = round(max(groupdata[[j]][, v]), 2)
			out[paste0(v, '_Mean'), paste0('Group_', j)]   = round(mean(groupdata[[j]][, v]), 2)
			out[paste0(v, '_Median'), paste0('Group_', j)] = round(median(groupdata[[j]][, v]), 2)
			out[paste0(v, '_Stdev'), paste0('Group_', j)]  = round(sd(groupdata[[j]][, v]), 2)
		}
	}
	# t-test, wilcox-test
	for (v in vars) {
		ttest = t.test(alvdata[, v], deadata[, v])
		out[paste0(v, '_Min'), 'T_P_AvsD'] = formatC(ttest$p.value, format = "E", digits = 2)
		out[paste0(v, '_Min'), 'T_ConfInt_AvsD'] = paste(round(ttest$conf.int, 2), collapse = ',')
		out[paste0(v, '_Min'), 'T_Means_AvsD'] = paste(round(ttest$estimate, 2), collapse = ',')
		wtest = wilcox.test(alvdata[, v], deadata[, v])
		out[paste0(v, '_Min'), 'Wilcox_P_AvsD'] = formatC(wtest$p.value, format = "E", digits = 2)
	}
	for (j in 1:ngroups) {
		for (k in 1:ngroups) {
			if (k >= j) next
			for (v in vars) {
				tryCatch({
					ttest = t.test(groupdata[[j]][, v], groupdata[[k]][, v])
					out[paste0(v, '_Min'), paste0('T_P_', j, 'vs', k)] = formatC(ttest$p.value, format = "E", digits = 2)
					out[paste0(v, '_Min'), paste0('T_ConfInt_', j, 'vs', k)] = paste(round(ttest$conf.int, 2), collapse = ',')
					out[paste0(v, '_Min'), paste0('T_Means_', j, 'vs', k)] = paste(round(ttest$estimate, 2), collapse = ',')
				}, error = function(e) {
					out[paste0(v, '_Min'), paste0('T_P_', j, 'vs', k)] = 'NA'
				})
				tryCatch({
				wtest = wilcox.test(groupdata[[j]][, v], groupdata[[k]][, v])
					out[paste0(v, '_Min'), paste0('Wilcox_P_', j, 'vs', k)] = formatC(wtest$p.value, format = "E", digits = 2)
				}, error = function(e){
					out[paste0(v, '_Min'), paste0('Wilcox_P_', j, 'vs', k)] = 'NA'
				})
			}
		}
	}
	# chi2-test
	for (v in vars) {

		if (length(levels(as.factor(vardata[, v]))) > chi2n) {
			out[paste0(v, '_Min'), 'Chi2_P_AvsD'] = 'NA'
			out[paste0(v, '_Min'), 'Chi2_P_Groups'] = 'NA'
			next
		}

		ctable = NULL
		tdata  = as.data.frame(table(alvdata[, v]))
		rownames(tdata) = as.character(tdata$Var1)
		tdata  = tdata[, -1, drop = F]
		colnames(tdata) = 'Alive'
		ctable = cbind.fill(ctable, tdata)
		tdata  = as.data.frame(table(deadata[, v]))
		rownames(tdata) = as.character(tdata$Var1)
		tdata  = tdata[, -1, drop = F]
		colnames(tdata) = 'Deceased'
		ctable = cbind.fill(ctable, tdata)
		ctable[is.na(ctable)] = 0
		ctable = ctable[order(rownames(ctable)), , drop = F]
		tryCatch({
			ctest  = chisq.test(ctable)
			out[paste0(v, '_Min'), 'Chi2_P_AvsD']    = formatC(ctest$p.value, format = "E", digits = 2)
			out[paste0(v, '_Min'), 'Chi2_Chi2_AvsD'] = as.character(round(ctest$statistic, 2))
			out[paste0(v, '_Min'), 'Chi2_Obs_AvsD']  = ctable.summary(ctest$observed)
			out[paste0(v, '_Min'), 'Chi2_Exp_AvsD']  = paste(apply(round(ctest$expected, 2), 1, function(r) paste(r, collapse = ',')), collapse = '; ')
		}, error = function(e) {
			out[paste0(v, '_Min'), 'Chi2_P_AvsD']    = 'NA'
		})
		
		ctable   = NULL
		for (j in 1:ngroups) {
			tdata   = as.data.frame(table(groupdata[[j]][, v]))
			rownames(tdata) = as.character(tdata$Var1)
			tdata = tdata[, -1, drop = F]
			colnames(tdata) = paste0("Group", j)
			ctable  = cbind.fill(ctable, tdata)
		}
		ctable[is.na(ctable)] = 0
		ctable = ctable[order(rownames(ctable)), , drop = F]
		tryCatch({
		ctest  = chisq.test(ctable)
			out[paste0(v, '_Min'), 'Chi2_P_Groups']    = formatC(ctest$p.value, format = "E", digits = 2)
			out[paste0(v, '_Min'), 'Chi2_Chi2_Groups'] = as.character(round(ctest$statistic, 2))
			out[paste0(v, '_Min'), 'Chi2_Obs_Groups']  = ctable.summary(ctest$observed)
			out[paste0(v, '_Min'), 'Chi2_Exp_Groups']  = paste(apply(round(ctest$expected, 2), 1, function(r) paste(r, collapse = ',')), collapse = '; ')
		}, error = function(e){
			out[paste0(v, '_Min'), 'Chi2_P_Groups']    = 'NA'
		})
	}
	sheet = paste0(var)
	addWorksheet(wb, sheet)
	writeData(wb, sheet, out, startCol = 2, startRow = 2, rowNames = T)
	colors = lighten(rainbow(length(vars) + ncomb + 14))
	for (j in 1:length(vars)) {
		mergeCells(wb, sheet, cols = 1, rows = (3 + 5 * (j - 1)):(3 + 5 * j - 1)) 
		writeData(wb, sheet, vars[j], startCol = 1, startRow = 3 + 5 * (j - 1))
		style = createStyle(fgFill = substring(colors[j], 1, 7), valign = "center", textDecoration = "bold", border="Bottom", borderColour = "#000000")
		#addStyle(wb, sheet, style, cols = 1, rows = 3 + 5 * (j - 1))
		addStyle(wb, sheet, style, cols = 1:2, rows = (3 + 5 * (j - 1)):(3 + 5 * j - 1), gridExpand = TRUE)

		style = createStyle(fgFill = substring(colors[j], 1, 7), valign = "center", border="LeftBottom", borderColour = "#000000")
		# total columns: 3 + ngroups + 3 + ncomb * 3 + 1 + ncomb + 8
		addStyle(wb, sheet, style, cols = 3:(17+ngroups+ncomb*4), rows = (3 + 5 * (j - 1)):(3 + 5 * j - 1), gridExpand = TRUE)

		for (k in (6+ngroups):(17+ngroups+ncomb*4)) {
			mergeCells(wb, sheet, cols = k, rows = (3 + 5 * (j - 1)):(3 + 5 * j - 1)) 
		}
	}
	mergeCells(wb, sheet, cols = 3:5, rows = 1)
	writeData(wb, sheet, 'Samples', startCol = 3, startRow = 1)
	style = createStyle(fgFill = substring(colors[length(vars)+1], 1, 7), halign = "center", textDecoration = "bold")
	addStyle(wb, sheet, style, cols = 3:5, rows = 1:2, gridExpand = TRUE)

	mergeCells(wb, sheet, cols = 6:(5+ngroups), rows = 1)
	writeData(wb, sheet, 'Groups', startCol = 6, startRow = 1)
	style = createStyle(fgFill = substring(colors[length(vars)+2], 1, 7), halign = "center", textDecoration = "bold")
	addStyle(wb, sheet, style, cols = 6:(5+ngroups), rows = 1:2, gridExpand = TRUE)

	mergeCells(wb, sheet, cols = (6+ngroups):(8+ngroups), rows = 1)
	writeData(wb, sheet, 'T-Test (Alive vs Deceased)', startCol = 6+ngroups, startRow = 1)
	style = createStyle(fgFill = substring(colors[length(vars)+3], 1, 7), halign = "center", textDecoration = "bold")
	addStyle(wb, sheet, style, cols = (6+ngroups):(8+ngroups), rows = 1:2, gridExpand = TRUE)

	mergeCells(wb, sheet, cols = (9+ngroups):(8+ngroups+ncomb*3), rows = 1)
	writeData(wb, sheet, 'T-Test', startCol = 9+ngroups, startRow = 1)
	style = createStyle(fgFill = substring(colors[length(vars)+4], 1, 7), halign = "center", textDecoration = "bold")
	addStyle(wb, sheet, style, cols = (9+ngroups):(8+ngroups+ncomb*3), rows = 1:2, gridExpand = TRUE)
	
	mergeCells(wb, sheet, cols = (9+ngroups+ncomb*3):(9+ngroups+ncomb*4), rows = 1)
	writeData(wb, sheet, 'Wilcox-Test', startCol = 9+ngroups+ncomb*3, startRow = 1)
	style = createStyle(fgFill = substring(colors[length(vars)+5], 1, 7), halign = "center", textDecoration = "bold")
	addStyle(wb, sheet, style, cols = (9+ngroups+ncomb*3):(9+ngroups+ncomb*4), rows = 1:2, gridExpand = TRUE)

	mergeCells(wb, sheet, cols = (10+ngroups+ncomb*4):(13+ngroups+ncomb*4), rows = 1)
	writeData(wb, sheet, 'Chi2-Test (Alive vs Deceased)', startCol = 10+ngroups+ncomb*4, startRow = 1)
	style = createStyle(fgFill = substring(colors[length(vars)+6], 1, 7), halign = "center", textDecoration = "bold")
	addStyle(wb, sheet, style, cols = (10+ngroups+ncomb*4):(13+ngroups+ncomb*4), rows = 1:2, gridExpand = TRUE)

	mergeCells(wb, sheet, cols = (14+ngroups+ncomb*4):(17+ngroups+ncomb*4), rows = 1)
	writeData(wb, sheet, 'Chi2-Test', startCol = 14+ngroups+ncomb*4, startRow = 1)
	style = createStyle(fgFill = substring(colors[length(vars)+7], 1, 7), halign = "center", textDecoration = "bold")
	addStyle(wb, sheet, style, cols = (14+ngroups+ncomb*4):(17+ngroups+ncomb*4), rows = 1:2, gridExpand = TRUE)
	
	setColWidths(wb, sheet, 1:(17+ngroups+ncomb*4), widths = 'auto')

}

saveWorkbook(wb, outfile, overwrite = TRUE)



