
import logging
import sys
import shlex
import re
from pyppl.utils import cmd
from glob import glob
from os import path
from collections import OrderedDict

class RuncmdException(Exception):
	pass

class Mem2Exception(Exception):
	pass

def fs2name(files):
	if not files: return 'nothing.etc'
	bnames  = [path.basename(f).split('.')[0] for f in files]
	compfix = path.commonprefix(bnames)
	compfix = re.sub(r'[^a-zA-Z0-9]$', '', compfix)
	if len(compfix) < 3:
		return bnames[0] + '.etc'
	else:
		return compfix + '.etc'

def dirpat2name(directory, pattern = '*'):
	return fs2name(glob(path.join(directory, pattern)))

def getLogger(name = 'bioprocs', logfmt = "[%(asctime)s][%(levelname)7s] %(message)s", level = logging.INFO):
	logger = logging.getLogger(name)
	for handler in logger.handlers:
		handler.close()
	del logger.handlers[:]

	ch = logging.StreamHandler()
	ch.setFormatter(logging.Formatter(fmt = logfmt))
	logger.addHandler(ch)

	logger.setLevel(level)
	return logger

def log2pyppl(msg, level = None):
	level = '.' + level if level else ''
	sys.stderr.write('pyppl.log%s:%s\n' % (level, msg))

def alwaysList(l):
	ret = []
	if isinstance(l, (list, tuple)):
		for x in l:
			ret += alwaysList(x)
	else:
		ret = [x.strip() for x in l.split(',') if x.strip()]
	return ret

def runcmd(cmd2run, shell = True, quit = True, ret = 'rc'):
	c = cmd.Cmd(cmd2run, shell = shell)
	cmdstr = ' '.join(c.cmd) if isinstance(c.cmd, list) else c.cmd
	logger.info('Running command at PID: %s' % c.pid)
	logger.info(cmdstr)
	#c.run()
	#for line in c.stderr.splitlines():
	
	# while c.p.poll() is None:
	# 	# blocking!!!
	# 	line = c.p.stderr.readline()
	# 	if not line:
	# 		continue
	# 	logger.error('STDERR: {}'.format(line.rstrip()))
	(stdout, stderr) = c.p.communicate()
	if ret == 'rc':
		# don't output if it's returned
		for line in stdout.splitlines():
			sys.stdout.write(line)
	for line in stderr.splitlines():
		logger.error('STDERR: {}'.format(line.rstrip()))

	c.rc = c.p.returncode
	logger.info('-'*80)
	logger.info('Return code: %s' % c.rc)
	if quit and c.rc != 0:
		raise RuncmdException('Command failed to run:\n{}\n'.format(cmdstr))
	return c.rc == 0 if ret == 'rc' else stdout

def call(command, args, quit = True):
	bioprocs = path.join(path.realpath(path.dirname(path.dirname(path.dirname(__file__)))), 'bin', 'bioprocs')
	if not isinstance(args, list):
		args = shlex.split(args)
	args = [bioprocs, command] + args
	return runcmd(args, quit)

def cmdargs(params, dash = 'auto', equal = 'auto', duplistkey = False, ignorefalse = True):
	"""
	Convert a dict of parameters into a command line parameter
	@params:
		`params`: The dict of parameters
		`dash`  : The dash prefix of each parameter key.
			- `auto`: A single dash for single-char key, otherwise double dashes
			- Otherwise use whatever being passed.
		`equal` : The separater between the key and the value.
			- `auto`: Space for single-char key, otherwise equal sign `=`
			- Otherwise use whatever being passed.
		`duplistkey`: Whether duplicate keys for list values. Default: `False`
			- For example for `params = {'a': [1,2,3]}`,
			- if `duplistkey = True`, it will be: `-a 1 -a 2 -a 3`,
			- otherwise it will be `-a 1 2 3`
		`ignorefalse`: Ignore boolean False parameter? Default: `True`
	"""
	if not params: return ''
	try:  # py3
		from shlex import quote
	except ImportError:  # py2
		from pipes import quote

	ret = []
	# decide the order
	if not isinstance(params, OrderedDict):
		params = OrderedDict([(key, params[key]) for key in sorted(params.keys())])
	
	# positional parameters
	positional = None
	if "" in params:
		positional = params[""]
		del params[""]
	
	for key, val in params.items():
		# allow comments for parameters for the same key
		key = key.split('#')[0].strip()
		item = dash if dash != 'auto' else '--' if len(key) > 1 else '-'
		item += key
		if isinstance(val, (tuple, list)):
			# ignore keys with no values
			if not val: continue
			item += equal if equal != 'auto' else '=' if len(key)>1 else ' '
			if duplistkey:
				ret.extend([item + quote(str(v)) for v in val])
			else:
				ret.append(item + quote(str(val[0])))
				ret.extend([quote(str(v)) for v in val[1:]])
		elif isinstance(val, bool):
			if not val and ignorefalse:
				continue
			if val:
				ret.append(item)
			else:
				item += equal if equal != 'auto' else '=' if len(key)>1 else ' '
				item += '0'
				ret.append(item)
		else:
			item += equal if equal != 'auto' else '=' if len(key)>1 else ' '
			item += quote(str(val))
			ret.append(item)

	if positional:
		if not isinstance(positional, (tuple, list)):
			positional = [positional]
		ret.extend([quote(str(pos)) for pos in positional])

	return ' '.join(ret)

def _autoUnit(num):
	if num % (1024 * 1024) == 0:
		return num / (1024*1024), 'G'
	elif num % 1024 == 0:
		return num / 1024, 'M'
	else:
		return num, 'K'

# unit = 'auto'/'G'/'M'/'K'/'java'
def mem2 (mem, unit = 'auto'):
	mem   = str(mem)
	ounit = mem[-1].upper()
	if ounit == 'G':
		num   = int(mem[:-1]) * 1024 * 1024
	elif ounit == 'M':
		num   = int(mem[:-1]) * 1024
	elif ounit == 'K':
		num   = int(mem[:-1])
	elif not ounit.isdigit():
		raise Mem2Exception('Unknown memory unit: ' + ounit)
	else:
		num   = int(mem)

	unit = unit.upper()
	retn, retu = _autoUnit(num)
	if unit == 'G':
		if retu == 'M':   retn /= 1024.
		elif retu == 'K': retn /= (1024. * 1024.)
		retu = 'G'
	elif unit == 'M':
		if retu == 'G':   retn *= 1024
		elif retu == 'K': ret /= 1024.
		retu = 'M'
	elif unit == 'K':
		if retu == 'G':   retn *= 1024 * 1024
		elif retu == 'M': retn *= 1024
		retu = 'K'

	if unit == 'JAVA':
		xmx = "-Xmx" + str(retn) + retu
		n, u = _autoUnit(num / 8)
		return '-Xms' + str(n) + u + ' ' + xmx
	else:
		return str(retn) + retu

def regionOverlap(CHR1, START1, END1, CHR2, START2, END2):
	if CHR1 != CHR2: return False
	if int(END1) < int(START2): return False
	if int(END2) < int(START1): return False
	return True

def funcargs(func):
	if not callable(func):
		raise ValueError('Expect a callable.')
	try:
		from inspect import signature
		return list(signature(func).parameters.keys())
	except ImportError:
		from inspect import getargspec
		return getargspec(func).args

logger = getLogger()
