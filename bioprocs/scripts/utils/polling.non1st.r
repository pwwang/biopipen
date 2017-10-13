if (!exists('pollingNon1st')) {
	pollingNon1st = function (jobid, cmd, flagfile, t = 10) {
		errorfile  = paste0(flagfile, '.error')
		
		if (jobid == 0) {
			tryCatch ({
				runcmd (cmd)
				file.create (flagfile)
			}, error = function(cond) {
				file.create (errorfile)
			})
		} else {

			while (TRUE) {
				if (file.exists(errorfile)) {
					stop('Error happend for job #0!')
				}
				if (file.exists(flagfile)) {
					return ()
				}
				write('Waiting for job #0 ...', stderr())
				Sys.sleep(t)
			}

		}
	}
}