from subprocess import Popen
from sys import stderr
if 'runcmd' not in vars() or not callable(runcmd):
	def runcmd (cmd, quit = True):
		stderr.write ("%s\n" % ('-'*80))
		stderr.write ("Running command:\n")
		stderr.write (cmd + "\n")
		rc = Popen (cmd, shell=True).wait()
		if rc == 0:
			stderr.write ("Return code: 0\n")
			stderr.write ("%s\n" % ('-'*80))
		else:
			stderr.write ("Return code: %s\n" % str(rc))
			stderr.write ("%s\n" % ('-'*80))
			if quit: exit(1)
		return rc == 0