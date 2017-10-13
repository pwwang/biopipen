if (!exists('params2CmdArgs')) {
	params2CmdArgs = function (params, dash = 'auto', equal = 'auto', noq = c()) {
		ret = c()
		for (key in names(params)) {
			val = params[[key]]
			key = trimws(key)
			if (is.logical(val) && !val) next
			item = if (dash == 'auto' && nchar(key) > 1) '--' else if (dash == 'auto') '-' else dash
			item = paste0(item, key)
			if (!is.logical(val)) {
				item = paste0(item, if (equal == 'auto' && nchar(key)>1) '=' else if (equal == 'auto') ' ' else equal)
				if (length(noq) == 0 || !(key %in% noq)) {
					item = paste0(item, shQuote(val))
				} else {
					item = paste0(item, val)
				}
			}
			ret = c(ret, item)
		}
		return(paste(ret, collapse=' '))
	}

}