from pyppl import Box
from pyppl.utils import cmd
from collections import OrderedDict
from subprocess import list2cmdline

try:  # py3
	from shlex import quote
except ImportError:  # py2
	from pipes import quote

shquote = lambda x: x[6:] if str(x).startswith('shraw:') else quote(str(x))

TOOLS = Box()

class RuncmdException(Exception):
	pass

def cmdargs(params, dash = 'auto', equal = 'auto', duplistkey = False, ignorefalse = True):
	"""
	Convert a dict of parameters into a command line parameter
	@params:
		`params`: The dict of parameters
			- positional arguments will be with key "_"
			- positional arguments before all other arguments will be with key ""
			- argument with key "_stdout" will redirct the stdout to file:
				- {'_stdout': 'outfile'} => " > 'outfile'"
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

	ret = []
	# decide the order
	params2 = OrderedDict()
	if '' in params:
		params2['']  = params['']
		del params['']
	if isinstance(params, OrderedDict):
		for key, val in params.items():
			if key == '_':
				continue
			params2[key] = val
		if '_' in params:
			params2['_'] = params['_']
	else:
		for key in sorted(params.keys()):
			if key == '_':
				continue
			params2[key] = params[key]
		if '_' in params:
			params2['_'] = params['_']

	params = params2
	outfile = None
	if '_stdout' in params:
		outfile = params['_stdout']
		del params['_stdout']
	aoutfile = None
	if '__stdout' in params:
		aoutfile = params['__stdout']
		del params['__stdout']
	elif '_stdout_' in params:
		aoutfile = params['_stdout_']
		del params['_stdout_']
	if outfile and aoutfile:
		raise ValueError('Cannot have both out files to write and append to.')
	for key, val in params.items():
		if key in ['', '_']:
			if not isinstance(val, (list, tuple)):
				val = [val]
			ret.extend([shquote(v) for v in val])
			continue
		# allow comments for parameters for the same key
		key = key.split('#')[0].strip()
		item = dash if dash != 'auto' else '--' if len(key) > 1 else '-'
		item += key
		if isinstance(val, (tuple, list)):
			# ignore keys with no values
			if not val: continue
			item += equal if equal != 'auto' else '=' if len(key)>1 else ' '
			if duplistkey:
				ret.extend([item + shquote(v) for v in val])
			else:
				ret.append(item + shquote(str(val[0])))
				ret.extend([shquote(str(v)) for v in val[1:]])
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
			item += shquote(val)
			ret.append(item)

	if outfile:
		ret.append('>')
		ret.append(shquote(outfile))
	elif aoutfile:
		ret.append('>>')
		ret.append(shquote(aoutfile))
	return ' '.join(ret)

def runcmd(cmd2run, raiseExc = True, **kwargs):
	from . import logger
	logger.info('RUNNING: %s', list2cmdline(cmd2run) if isinstance(cmd2run, list) else cmd2run)
	# stdin to avoid blocking
	kwargs['shell'] = kwargs.get('shell', True)
	c = cmd.Cmd(cmd2run, stdin = None, **kwargs)
	logger.info('- PID: %s', c.pid)

	for line in c.readline():
		logger.info('- STDERR: %s', line.rstrip())

	logger.info('- RETURNCODE: %s', c.rc)
	logger.info('-' * 80)
	if raiseExc and c.rc != 0:
		raise RuncmdException(c.cmd)
	return c

def runcmd2(tool, params, dash = 'auto', equal = 'auto', duplistkey = False, ignorefalse = True, quit = True):
	s = Shell(tools = {'tool': tool}, dash = dash, equal = equal, duplistkey = duplistkey, ignorefalse = ignorefalse)
	r = s.tool(**params).run()
	if quit and r.rc != 0:
		raise RuncmdException(r.cmd)
	return r

# def call(command, args, quit = True):
# 	bioprocs = path.join(path.realpath(path.dirname(path.dirname(path.dirname(__file__)))), 'bin', 'bioprocs')
# 	if not isinstance(args, list):
# 		args = shlex.split(args)
# 	args = [bioprocs, command] + args
# 	return runcmd(args, quit)

class ShellResult(object):

	def __init__(self, cmdobj, tools = None):
		self.cmdobj = cmdobj
		self.tools  = tools or {}
		self.done   = False
	
	def __str__(self):
		return self.stdout

	def __repr__(self):
		return '<ShellResult {!r}>'.format(self.cmdobj)

	def run(self, raiseExc = True, save = None, logger = 'bioprocs'):
		if not self.done:
			if logger == 'bioprocs':
				from . import logger
				logit = lambda *args: logger.info(*args)
			elif logger == 'sys':
				from sys import stderr
				logit = lambda *args: stderr.write((args[0] + '\n') % args[1:])
			else:
				logit = lambda *args: None
			logit('RUNNING: %s', self.cmdobj.cmd)
			logit('- PID: %s', self.cmdobj.pid)

			for line in self.cmdobj.readline(save = save):
				logit('- STDERR: %s', line.rstrip())
			self.done = True
			logit('- RETURNCODE: %s', self.cmdobj.rc)
			logit('-' * 80)
			if raiseExc and self.cmdobj.rc != 0:
				raise RuncmdException(self.cmdobj.cmd)
		return self

	@property
	def stdout(self):
		return self.run(save = 'stdout').cmdobj.stdout
	
	@property
	def stderr(self):
		return self.run(save = 'stderr').cmdobj.stderr

	@property
	def cmd(self):
		return self.cmdobj.cmd

	@property
	def returncode(self):
		return self.run().cmdobj.rc
	
	@property
	def rc(self):
		return self.returncode

	def pipe(self, 
		tools = None,   subcmd = False,
		dash  = 'auto', equal  = 'auto', duplistkey = False, ignorefalse = True,
		base  = None,   **pargs):
		tools = tools or {}
		tools.update(self.tools)
		return Shell(
			tools       = tools,
			subcmd      = subcmd,
			dash        = dash,
			equal       = equal,
			duplistkey  = duplistkey,
			ignorefalse = ignorefalse,
			base        = None,
			stdin       = self.cmdobj.p.stdout,
			prevcmd     = self.cmd,
			**pargs
		)

class Shell(object):
	def __init__(self, 
		tools = None,   subcmd  = False,
		dash  = 'auto', equal   = 'auto', duplistkey = False, ignorefalse = True,
		base  = None,   prevcmd = None,   **pargs):
		self.dash        = dash
		self.equal       = equal
		self.duplistkey  = duplistkey
		self.ignorefalse = ignorefalse
		self.subcmd      = subcmd
		self.tools       = TOOLS
		self.base        = base
		self.prevcmd     = prevcmd

		self.tools.update(tools or {})
		# common args for tools with subcommand
		self.targs         = {}
		self._subcmdCalled = False

		self.pargs = pargs.copy()
		self.pargs['stdin'] = self.pargs.get('stdin')
		self.pargs['shell'] = True

	def __call__(self, *args, **kwargs):
		self.targs.update(kwargs)
		self.targs[''] = list(args) + self.targs.get('', [])
		return self

	def _run(self, *args, **kwargs):
		kwargs[''] = list(args) + kwargs.get('', []) 
		targs = '' if not self.targs else ' ' + cmdargs(
			self.targs, dash = self.dash, equal = self.equal, 
			duplistkey = self.duplistkey, ignorefalse = self.ignorefalse
		)
		cargs = cmdargs(
			kwargs,     dash = self.dash, equal = self.equal, 
			duplistkey = self.duplistkey, ignorefalse = self.ignorefalse
		)
		cmd2run = '{0}{1} {2}'.format(self.base, targs, cargs)
		cmdobj  = cmd.Cmd(cmd2run, **self.pargs)
		
		if self.prevcmd:
			cmdobj.cmd = '{} | {}'.format(self.prevcmd, cmdobj.cmd)
		return ShellResult(cmdobj, self.tools)

	def __getitem__(self, tool):
		return getattr(self, tool)

	def __getattr__(self, tool):
		tool = self.tools.get(tool, tool)
		# there is not subcommand or now it's already in subcommand
		if not self.subcmd:
			# in a subcommand
			if self.base:
				if not self._subcmdCalled:
					self.targs['_'] = self.targs.get('_', []) + [tool]
					self._subcmdCalled = True
				else:
					self.targs['_'].pop(-1)
					self.targs['_'].append(tool)
			# plain command
			else:
				self.base = shquote(tool)
			return self._run
		return Shell(
			subcmd = False, dash = self.dash, equal = self.equal, duplistkey = self.duplistkey, 
			ignorefalse = self.ignorefalse, base = tool, **self.pargs)

wc   = lambda *args, **kwargs: Shell().wc(*args, **kwargs).run(save = 'stdout')
wc_l = lambda filename: int(wc(l = filename).stdout.split()[0])
wcl  = wc_l

gzip      = lambda *args,  **kwargs: Shell().gzip(*args, **kwargs).run(logger = False)
bgzip     = lambda *args,  **kwargs: Shell(equal = ' ').bgzip(*args, **kwargs).run(logger = False)
gunzip    = lambda *args,  **kwargs: Shell().gunzip(*args, **kwargs).run(logger = False)
gzip_to   = lambda source, dest: gzip(_ = source, c = True, _stdout = dest)
bgzip_to  = lambda source, dest: bgzip(_ = source, c = True, _stdout = dest)
gunzip_to = lambda source, dest: gunzip(_ = source, c = True, _stdout = dest)

mv    = lambda *args,  **kwargs: Shell().mv(*args, **kwargs).run(logger = False)
cp    = lambda *args,  **kwargs: Shell().cp(*args, **kwargs).run(logger = False)
rm    = lambda *args,  **kwargs: Shell().rm(*args, **kwargs).run(logger = False)
rm_rf = lambda *args: rm(*args, r = True, f = True)
rmrf  = rm_rf
mkdir = lambda *args,  **kwargs: Shell().mkdir(*args, **kwargs).run(logger = False)
ln    = lambda *args,  **kwargs: Shell().ln(*args, **kwargs).run(logger = False)
ln_s  = lambda source, dest: ln(source, dest, s = True)

head  = lambda *args, **kwargs: Shell().head(*args, **kwargs).run(logger = False)
tail  = lambda *args, **kwargs: Shell().tail(*args, **kwargs).run(logger = False)
grep  = lambda *args, **kwargs: Shell().grep(*args, **kwargs).run(raiseExc = False, logger = False)
sort  = lambda *args, **kwargs: Shell(dash = '-', equal = ' ', duplistkey = True).sort(*args, **kwargs).run(logger = False)
uniq  = lambda *args, **kwargs: Shell().uniq(*args, **kwargs).run(logger = False)
touch = lambda *args, **kwargs: Shell().touch(*args, **kwargs).run(logger = False)

kill  = lambda *args, **kwargs: Shell().kill(*args, **kwargs).run(logger = False)
kill9 = lambda *procids: kill(**{'9': True, '': procids})

zcat  = lambda *args,  **kwargs: Shell().zcat(*args, **kwargs).run(logger = False)
cat   = lambda *args,  **kwargs: Shell().cat(*args, **kwargs).run(logger = False)
# _ means waiting for run
acat_ = lambda infile, **kwargs: Shell().zcat(infile, **kwargs) if infile.endswith('.gz') else Shell().cat(infile, **kwargs)

which = lambda *args, **kwargs: Shell().which(*args, **kwargs).run(logger = False, save = 'stdout').stdout.rstrip('\n')
