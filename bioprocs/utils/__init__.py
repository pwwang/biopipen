
import logging
import sys
import shlex
import re
import gzip
from glob import glob
from os import path
import yaml
from pyppl.utils import Box

class Mem2Exception(Exception):
	pass

def ensureBox(var):
	if isinstance(var, Box):
		return var
	return Box(var)

def fs2name(files):
	if not files: return 'nothing.etc'
	bnames  = [path.basename(f).split('.')[0] for f in files]
	compfix = path.commonprefix(bnames)
	compfix = re.sub(r'[^a-zA-Z0-9]$', '', compfix)
	if len(compfix) < 3:
		return bnames[0] + '.etc'
	else:
		return compfix + '.etc'

def reportdata(outdir, *args, **kwargs):
	data = {}
	for arg in args:
		if not isinstance(arg, dict):
			raise TypeError('Expect a dict as report data.')
		data.update(arg)
	data.update(kwargs)
	datafile = path.join(outdir, 'job.report.data.yaml')
	with open(datafile, 'w') as f:
		yaml.dump(data, f, default_flow_style=False)

def dirpat2name(directory, pattern = '*'):
	return fs2name(glob(path.join(directory, pattern)))

def getLogger(name = 'bioprocs', logfmt = "[%(asctime)s %(levelname)7s] %(message)s", datafmt = "%Y-%m-%d %H:%M:%S", level = logging.INFO):
	logger = logging.getLogger(name)
	for handler in logger.handlers:
		handler.close()
	del logger.handlers[:]

	ch = logging.StreamHandler()
	ch.setFormatter(logging.Formatter(logfmt, datafmt))
	logger.addHandler(ch)

	logger.setLevel(level)
	return logger

def log2pyppl(msg, level = None):
	level = '.' + level if level else ''
	sys.stderr.write('pyppl.log%s:%s\n' % (level, msg))

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
		xmx = "-Xmx" + str(int(retn)) + retu
		n, u = _autoUnit(num / 8)
		return '-Xms' + str(int(n)) + u + ' ' + xmx
	elif unit == 'JAVADICT' or unit == 'JDICT':
		n, u = _autoUnit(num / 8)
		return {
			'Xmx' + str(retn) + retu: True,
			'Xms' + str(n) + u: True
		}
	elif unit == '-JAVADICT' or unit == '-JDICT':
		n, u = _autoUnit(num / 8)
		return {
			'-Xmx' + str(retn) + retu: True,
			'-Xms' + str(n) + u: True
		}
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

def gztype(gzfile):
	import binascii
	with open(gzfile, 'rb') as f:
		flag = binascii.hexlify(f.read(4))
	if flag == b'1f8b0804':
		return 'bgzip'
	if flag == b'1f8b0808':
		return 'gzip'
	return 'flat'

class FileConn(object):

	def __init__(self, filename, flag = 'r'):
		self.isgz = False
		if filename.endswith('.gz'):
			self.isgz = True
		self.filename = filename
		self.flag     = flag
		self.fp       = None

	def __enter__(self):
		return self.open()

	def __exit__(self, exc_type, exc_val, exc_tb):
		self.close()

	def open(self):
		if self.isgz:
			self.fp = gzip.open(self.filename, self.flag)
		else:
			self.fp = open(self.filename, self.flag)
		return self.fp

	def close(self):
		self.fp.close()

def switch(value, cases, *args, **kwargs):
	ret = cases.get(value, 'default')
	if callable(ret) and (args or kwargs):
		return ret(*args, **kwargs)
	return ret

logger = getLogger()
