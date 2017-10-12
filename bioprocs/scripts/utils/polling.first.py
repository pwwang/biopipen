from os import path
from time import sleep
from sys import stderr

if 'pollingFirst' not in vars() or not callable(pollingFirst):
	def pollingFirst (jobid, cmd, flagfile, t = 10):
		errorfile = flagfile + '.error'
		
		if jobid == 0:
			try:
				runcmd (cmd)
				open (flagfile, 'w').close()
			except SystemExit:
				open (errorfile, 'w').close()
		else:
			while True:
				if path.exists(errorfile):
					stderr.write('Error happened for job #0!\n')
					exit (1)
				if path.exists(flagfile):
					return
				stderr.write('Waiting for job #0 ...\n')
				sleep(t)
				