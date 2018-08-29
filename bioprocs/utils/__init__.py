
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

def runcmd(cmd2run, shell = True, quit = True):
	c = cmd.Cmd(cmd2run, shell = shell)
	cmdstr = ' '.join(c.cmd) if isinstance(c.cmd, list) else c.cmd
	logger.info('Running command at PID: %s' % c.pid)
	logger.info(cmdstr)
	c.run()
	for line in c.stderr.splitlines():
		logger.error('STDERR: {}'.format(line))
	logger.info('Return code: %s' % c.rc)
	logger.info('-'*80)
	if quit and c.rc != 0:
		raise RuncmdException('Command failed to run:\n{}\n'.format(cmdstr))
	return c.rc == 0

def call(command, args, quit = True):
	bioprocs = path.join(path.realpath(path.dirname(path.dirname(path.dirname(__file__)))), 'bin', 'bioprocs')
	if not isinstance(args, list):
		args = shlex.split(args)
	args = [bioprocs, command] + args
	return runcmd(args, quit)

def cmdargs(params, dash = 'auto', equal = 'auto'):
	if not params: return ''
	try:  # py3
		from shlex import quote
	except ImportError:  # py2
		from pipes import quote
	ret = []
	if not isinstance(params, OrderedDict):
		params = OrderedDict([(key, params[key]) for key in sorted(params.keys())])

	for key, val in params.items():
		key = key.split('#')[0].strip()
		if isinstance(val, bool) and not val: continue
		item = '--' if (dash == 'auto' and len(key) > 1) else '-' if dash == 'auto' else dash
		item += key
		if isinstance(val, (tuple, list)):
			item += '=' if (equal == 'auto' and len(key)>1) else ' ' if equal == 'auto' else equal
			for i, v in enumerate(val):
				if i == 0:
					item += quote(str(v))
					ret.append(item)
				else:
					ret.append(quote(str(v)))
		elif not isinstance(val, bool):
			item += '=' if (equal == 'auto' and len(key)>1) else ' ' if equal == 'auto' else equal
			item += quote(str(val))
			ret.append(item)
		else:
			ret.append(item)
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

logger = getLogger()
