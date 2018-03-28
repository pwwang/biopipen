import sys
# https://github.com/kennethreitz/delegator.py
import delegator

class RuncmdException(Exception):
	pass

def runcmd(cmd, quit = True):
	c = delegator.run(cmd)
	sys.stderr.write('Running command at PID: %s\n' % c.pid)
	sys.stderr.write(cmd + '\n')
	c.block()
	sys.stderr.write('Return code: %s\n' % c.return_code)
	sys.stderr.write ('-'*80 + '\n')
	if quit and c.return_code != 0:
		raise RuncmdException('Command failed to run.')
	return c.return_code == 0
