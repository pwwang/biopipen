if (!exists('utils')) {
	# don't load them twice
	# python modkit package makes second load of utils$shell2 raise exceptions
	utils  = reticulate::import('bioprocs', delay_load = TRUE)$utils
	shell  = utils$shell2
	runcmd = shell$runcmd
	mem2   = utils$mem2
}
options(stringsAsFactors = FALSE)

cbindfill = function (...) {
	dfs = list(...)
	ret = Reduce(function(...) {
		ret = merge(..., by='row.names', all=T, sort = F)
		rownames(ret) = ret[, "Row.names"]
		ret = ret[, -1, drop=F]
		ret
	}, dfs)
	cnames = unlist(lapply(dfs, colnames))
	if (length(cnames) > 0) {
		colnames(ret) = make.unique(cnames)
	}
	return (ret)
}

rbindfill = function (...) {
	library(data.table)
	dfs = lapply(list(...), as.data.frame)
	ret = as.data.frame(rbindlist(dfs, fill = T, use.names = T))
	rnames = unlist(lapply(dfs, rownames))
	if (length(rnames) > 0) {
		rownames(ret) = make.unique(rnames)
	}
	return (ret)
}

logger = function(..., level = 'INFO') {
	cat(paste0(level, ': ', paste(...), '\n'), file = stderr())
}

log2pyppl = function(..., level = 'LOG') {
	msg = paste(...)
	if (!endsWith(msg, "\n")) msg = paste0(msg, "\n")
	cat(paste0('pyppl.log.', level, ': ', msg), file = stderr())
}

# avoid typos
cbind.fill = cbindfill
rbind.fill = rbindfill

read.table.nodup = function(...) {
	args = list(...)
	if (!'row.names' %in% names(args) || is.null(args$row.names)) {
		return(read.table(...))
	} else {
		# 'row.names' removed
		args$row.names = NULL
		# get it back
		args = c(args, list(row.names = NULL))
		mat = do.call(read.table, args)
		rnames = make.unique(as.character(as.vector(mat[,1,drop = T])))
		mat = mat[,-1,drop=F]
		row.names(mat)  = rnames
		return(mat)
	}
}

update.list = function (list1, list2, recursive = F) {
	names2 = names(list2)
	for (name in names2) {
		if (is.list(list1[[name]]) && is.list(list2[[name]]) && recursive) {
			list1[[name]] = update.list(list1[[name]], list2[[name]], recursive = recursive)
		} else {
			list1[[name]] = list2[[name]]
		}
	}
	return (list1)
}

# allow ifelse to return NULL
ifelse = function(condition, true, false) {
	if(condition) return(true)
	return(false)
}

read.table.inopts = function(infile, inopts) {
	inopts.default = function(key, default) list.get(inopts, key, default, check.names = TRUE)
	optrnames = inopts.default('rnames', TRUE)
	dup = inopts$dup
	try = list.get(inopts, 'try', FALSE)
	inopts$dup = NULL
	params = list(
		infile,
		#header      = ifelse('cnames' %in% opts, inopts$cnames, T),
		header      = inopts.default('cnames', TRUE),
		row.names   = ifelse(!is.null(dup), NULL, ifelse(optrnames, 1, NULL)),
		sep         = inopts.default('delimit', "\t"),
		check.names = F,
		quote       = inopts.default('quote', ""),
		skip        = inopts.default('skip', 0)
	)
	inopts$cnames  = NULL
	inopts$rnames  = NULL
	inopts$delimit = NULL
	inopts$quote   = NULL
	inopts$skip    = NULL
	params = c(params, inopts)
	if (!try) {
		d = do.call(read.table, params)
	} else {
		d = tryCatch({do.call(read.table, params)}, error = function(e) {
			logger('Error encountered while read file:', infile, level = 'WARNING')
			logger(e, level = 'WARNING')
			return (NULL)
		})
	}
	if (is.null(dup) || !optrnames || is.null(d)) {
		#return (d)
	} else {
		rnames = as.character(as.vector(d[,1,drop = T]))
		if (dup == 'drop') {
			rindex = !duplicated(rnames)
			d = d[rindex, -1, drop = FALSE]
			rownames(d) = rnames[rindex]
		} else if (dup == 'mean') {
			d = aggregate(d[,-1,drop=FALSE], by = list(rnames), FUN = mean)
			rownames(d) = as.character(d[,1,drop=TRUE])
			d = d[,-1,drop=FALSE]
		} else  { # keep
			rownames(d) = make.unique(rnames)
			d = d[,-1,drop = FALSE]
		}
	}
	d
}

write.xls = function(x, file, ...) {
	rnames         = rownames(x)
	args           = list(...)
	args$quote     = list.get(args, 'quote', FALSE)
	args$sep       = list.get(args, 'sep', "\t")
	args$row.names = list.get(args, 'row.names', TRUE)
	if (!is.null(rnames) && args$row.names == TRUE) {
		cnames = colnames(x)
		x = cbind(`.` = rnames, x)
		colnames(x) = c('.', cnames)
		args$row.names = FALSE
	}
	do.call(write.table, c(list(x, file), args))
}

# format data.frame to output
pretty.numbers = function(df, formats) {
	# remember set stringsAsFactors as FALSE for the dataframe!!
	if (nrow(df) == 0) {
		return(df)
	}
	allCols      = colnames(df)
	formatedCols = c()
	for (fcols in names(formats)) {
		if (fcols == '.') { # must be last element of formats
			cols = which(!allCols %in% formatedCols)
		} else {
			cols = unlist(strsplit(fcols, '..', fixed = T))
			formatedCols = c(formatedCols, cols)
		}
		cols = intersect(cols, allCols)
		df[, cols] = sprintf(formats[[fcols]], as.numeric(unlist(df[, cols])))
	}
	df
}

# format data.frame to output
pretty.numbers2 = function(df, ...) {
	formats = list(...)
	df = as.data.frame(df)
	if (nrow(df) == 0)
		return(df)

	allCols      = colnames(df)
	if (is.null(allCols))
		allCols = 1:ncol(df)
	formatedCols = c()
	for (fcols in names(formats)) {
		if (fcols == '.') { # must be last element of formats
			if (length(formatedCols) == 0) {
				cols = allCols
			} else {
				cols = which(!allCols %in% formatedCols)
			}
		} else {
			cols = unlist(strsplit(fcols, '..', fixed = T))
			formatedCols = c(formatedCols, cols)
		}
		cols = intersect(cols, allCols)
		df[, cols] = sprintf(formats[[fcols]], as.numeric(unlist(df[, cols])))
	}
	df
}

is.installed = function(pkg) {
	is.element(pkg, installed.packages()[,1])
}

.bQuote = function(s) {
	if (startsWith(s, '`') && endsWith(s, '`')) {
		return (s)
	} else {
		paste0('`', s, '`')
	}
}
bQuote = Vectorize(.bQuote)

.nobQuote = function(s) {
	if (!startsWith(s, '`') || !endsWith(s, '`'))
		return(s)
	return (substring(s, 2, nchar(s) - 1))
}
nobQuote = Vectorize(.nobQuote)

is.true = function(x, collapse = 'all') {
	if (is.null(x)) return (FALSE)
	if (length(x) == 0) return (FALSE)
	if (length(x) == 1) {
		if (is.na(x)) return (FALSE)
		if (is.list(x)) return (TRUE)
		if (is.character(x)) return (nchar(x) > 0)
		tryCatch({
			x = as.logical(x)
		}, error = function(e){
			x <<- TRUE
		})
		if (is.na(x)) return (TRUE)
		return (x)
	} else if (collapse == 'all') {
		for (i in x) {
			if (!is.true(i, 'any')) return (FALSE)
		}
		return (TRUE)
	} else {
		for (i in x) {
			if (is.true(i, 'any')) return (TRUE)
		}
		return (FALSE)
	}
}

is.false = function(x, collapse = 'all') {
	!is.true(x, ifelse(collapse == 'all', 'any', 'all'))
}

list.get = function(l, key, default = NULL, check.names = FALSE) {
	# get the value of a key in list with default
	# @params:
	#	`l`: The list
	#	`key`: The key
	#	`default`: The default value. Default: `NULL`
	#	`check.names`: Check whether the name exists, even with value `NULL`. Default: `FALSE`
	#		- `list.get(list(a = NULL), 'a', default = 1, check.names = TRUE) == NULL`
	#		- `list.get(list(a = NULL), 'a', default = 1, check.names = FALSE) == 1`
	if (!check.names) {
		ifelse(is.null(l[[key]]), default, l[[key]])
	} else {
		ns = names(l)
		if (key %in% ns)
			return (l[[key]])
		return (default)
	}
}
