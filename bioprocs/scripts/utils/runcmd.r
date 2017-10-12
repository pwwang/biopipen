if (!exists('runcmd')) {
	runcmd = function (cmd, quit = TRUE, intern=FALSE) {
		write (paste(rep('-', 80), collapse="") , stderr())
		write ("Running command: ", stderr())
		write (cmd, stderr())
		
		if (!intern) {
			rc = try(system (cmd))
			if (rc == 0) {
				write ("Return code: 0", stderr())
				write (paste(rep('-', 80), collapse="") , stderr())
			} else {
				write (paste0("Return code: ", rc), stderr())
				write (paste(rep('-', 80), collapse="") , stderr())
				if (quit) {
					stop ('Quitted.')
				}
			}
			return (rc)
		} else {
			return (system(cmd, intern=T))
		}
	}
}