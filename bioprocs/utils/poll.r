library(methods)
library(filelock)
library(reticulate)
bioprocs = import('bioprocs')
source(file.path(bioprocs$UTILS, '__init__.r'))

getkw = function (kwargs, key, default) {
	kwnames = names(kwargs)
	if (!key %in% kwnames) {
		value = default
	} else {
		value = kwargs[[key]]
		kwargs[[key]] = NULL
	}
	ret = list(kwargs = kwargs)
	ret[[key]] = value
	return (ret)
}

wait = function (func, ...) {
	ret = getkw(as.list(match.call()), 'interval', .1)
	ret$kwargs$func = NULL
	ret$kwargs[[1]] = NULL
	env = parent.frame()
	ret$kwargs = lapply(ret$kwargs, function(row) eval(row, envir = env))
	while (do.call(func, ret$kwargs)) {
		Sys.sleep(ret$interval)
	}
}

Poll = function (workdir, joblen, jobindex) {
	poll = list()
	poll$workdir  = workdir
	poll$joblen   = joblen
	poll$jobindex = jobindex

	poll$first    = function(todo, ...) {
		ret = getkw(as.list(match.call()), 'lockfile', 'poll.first.lock')
		ret$kwargs$todo = NULL
		ret$kwargs[[1]] = NULL
		lockfilename = file.path(workdir, '1', 'output', ret$lockfile)
		if (jobindex == 0) {
			lockfile = lock(lockfilename)
			if (is.function(todo)) {
				do.call(todo, ret$kwargs)
			} else {
				cmd = do.call(sprintf, c(list(todo), ret$kwargs))
				runcmd(cmd)
			}
		}
	}

	poll$non1st    = function(todo, ...) {
		ret = getkw(as.list(match.call()), 'lockfile', 'poll.non1st.lock')
		ret$kwargs$todo = NULL
		ret$kwargs[[1]] = NULL
		lockfilenames = sapply(1:joblen, function (x) file.path(workdir, x, 'output', ret$lockfile))
		if (jobindex == 0) {
			for (i in 1:joblen) {
				if (i == 1) next
				wait(function(x) !file.exists(x), lockfilenames[i])
				unlock(lock(lockfilenames[i]))
			}
		} else {
			lockfile = lock(lockfilenames[jobindex + 1])
			if (is.function(todo)) {
				do.call(todo, ret$kwargs)
			} else {
				cmd = do.call(sprintf, c(list(todo), ret$kwargs))
				runcmd(cmd)
			}
			unlock(lockfile)
		}
	}

	poll$all = function(todo, ...) {
		ret = getkw(as.list(match.call()), 'lockfile', 'poll.all.lock')
		ret$kwargs$todo = NULL
		ret$kwargs[[1]] = NULL
		lockfilenames = sapply(1:joblen, function (x) file.path(workdir, x, 'output', ret$lockfile))
		lockfile = lock(lockfilenames[jobindex + 1])
		if (is.function(todo)) {
			do.call(todo, ret$kwargs)
		} else {
			cmd = do.call(sprintf, c(list(todo), ret$kwargs))
			runcmd(cmd)
		}
		unlock(lockfile)
		for (i in 1:joblen) {
			if (i == jobindex + 1) next
			wait(function(x) !file.exists(x), lockfilenames[i])
			unlock(lock(lockfilenames[i]))
		}
	}

	return(poll)

}
