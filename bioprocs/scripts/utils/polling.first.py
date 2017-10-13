from os import path
from time import sleep
from sys import stderr, exit

if 'pollingFirst' not in vars() or not callable(pollingFirst):
	"""
	Waiting for non-first jobs done.
	@params:
		`workdir`: The workdir of the process
		`length`:  Number of jobs
		`jobid`:   The job index
		`cmd`:     The command for each job
		`flagname`:The flag filename
		`t`:       The interval
	"""
	def pollingFirst (workdir, length, jobid, cmd, flagfname, t = 10):
		errorfname = flagfname + '.error'
		flagfile   = path.join(workdir, str(jobid), 'output', flagfname)
		errorfile  = path.join(workdir, str(jobid), 'output', errorfname)
		
		try:
			runcmd (cmd)
			open (flagfile, 'w').close()
		except SystemExit:
			open (errorfile, 'w').close()
		
		if jobid > 0: return # only job #0 wait

		wait = True
		while wait:
			wait = False
			for i in range(length):
				efile = path.join(workdir, str(i), 'output', errorfname)
				ffile = path.join(workdir, str(i), 'output', flagfname)
				if path.exists (efile):
					stderr.write ('Job '+ str(i) +' failed, I am also exiting ...')
					exit (1)
				if not path.exists (ffile):
					wait = True
					break
			if wait:
				stderr.write('Waiting till other jobs (job # != 0) done ...\n')
				sleep(t)
			else:
				break