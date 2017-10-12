if (!exists('mem2')) {

	autoUnit = function (num) {
		if (num %% (1024 * 1024) == 0) {
			return(list(n = num / 1024 / 1024, u = 'G'))
		} else if (num %% 1024 == 0) {
			return(list(n = num / 1024, u = 'M'))
		} else {
			return(list(n = num, u = 'K'))
		}
	}

	mem2 = function (mem, unit = 'auto') {

		if (grepl('G$', mem) || grepl('g$', mem)) {
			num = as.numeric (substr(mem, 1, nchar(mem)-1)) * 1024 * 1024
		} else if (grepl('M$', mem) || grepl('m$', mem)) {
			num = as.numeric (substr(mem, 1, nchar(mem)-1)) * 1024
		} else if (grepl('K$', mem) || grepl('k$', mem)) {
			num = as.numeric (substr(mem, 1, nchar(mem)-1))
		} else {
			num = as.numeric (mem)
		}

		ret = autoUnit(num)
		if (unit == 'G' || unit == 'g') {
			if (ret$u == 'M') ret$n = ret$n / 1024.0
			if (ret$u == 'K') ret$n = ret$n / 1024.0 / 1024.0
			ret$u = 'G'
		} else if (unit == 'M' || unit == 'm') {
			if (ret$u == 'G') ret$n = ret$n * 1024
			if (ret$u == 'K') ret$n = ret$n / 1024.0 
			ret$u = 'M'
		} else if (unit == 'K' || unit == 'k') {
			if (ret$u == 'G') ret$n = ret$n * 1024.0 * 1024.0
			if (ret$u == 'M') ret$n = ret$n * 1024.0
			ret$u = 'K'
		}

		if (unit == 'java' || unit == 'Java' || unit == 'JAVA') {
			xmx = paste0("-Xmx", ret$n, ret$u)
			r   = autoUnit(num / 8)
			return(paste(paste0('-Xms', r$n, r$u), xmx))
		} else {
			return(paste0(ret$n, ret$u))
		}
	}
}