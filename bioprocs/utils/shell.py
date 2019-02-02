from os import getcwd
from pyppl import Box
from pyppl.utils import cmd
from collections import OrderedDict

try:  # py3
	from shlex import quote as shquote
except ImportError:  # py2
	from pipes import quote as shquote

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
	if outfile and aoutfile:
		raise ValueError('Cannot have both out files to write and append to.')
	for key, val in params.items():
		if key in ['', '_']:
			if not isinstance(val, (list, tuple)):
				val = [val]
			ret.extend([shquote(str(v)) for v in val])
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
				ret.extend([item + shquote(str(v)) for v in val])
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
			item += shquote(str(val))
			ret.append(item)

	if outfile:
		ret.append('>')
		ret.append(shquote(outfile))
	elif aoutfile:
		ret.append('>>')
		ret.append(shquote(aoutfile))
	return ' '.join(ret)

def runcmd(cmd2run, shell = True, quit = True):
	from . import logger
	# stdin to avoid blocking
	c = cmd.Cmd(cmd2run, shell = shell, stdin = None)
	cmdstr = ' '.join(c.cmd) if isinstance(c.cmd, list) else c.cmd
	logger.info('Running command at PID: %s' % c.pid)
	logger.info(cmdstr)

	while c.p.poll() is None:
		line = c.p.stderr.readline()
		if not line:
			continue
		logger.info('STDERR: {}'.format(line.rstrip()))
	for line in c.p.stderr:
		logger.info('STDERR: {}'.format(line.rstrip()))

	c.rc = c.p.wait()
	logger.info('-'*80)
	logger.info('Return code: %s' % c.rc)
	if quit and c.rc != 0:
		raise RuncmdException(cmdstr)
	return c.rc == 0

def runcmd2(tool, params, dash = 'auto', equal = 'auto', duplistkey = False, ignorefalse = True, quit = True):
	s = Shell(tools = {'tool': tool}, dash = dash, equal = equal, duplistkey = duplistkey, ignorefalse = ignorefalse)
	r = s.tool(**params).run()
	if quit and r.rc != 0:
		raise RuncmdException(r.cmd)
	return r

def call(command, args, quit = True):
	bioprocs = path.join(path.realpath(path.dirname(path.dirname(path.dirname(__file__)))), 'bin', 'bioprocs')
	if not isinstance(args, list):
		args = shlex.split(args)
	args = [bioprocs, command] + args
	return runcmd(args, quit)

class ShellResult(object):

	def __init__(self, cmdobj, tools = None):
		self.cmdobj = cmdobj
		self.tools  = tools or {}
		self.done   = False
	
	def __str__(self):
		return self.stdout

	def __repr__(self):
		return '<ShellResult {!r}>'.format(self.cmdobj)

	def run(self, quit = True):
		if not self.done:
			from . import logger
			logger.info('Running command at PID: %s' % self.cmdobj.pid)
			logger.info(self.cmdobj.cmd)
			while self.cmdobj.p.poll() is None:
				line = self.cmdobj.p.stderr.readline()
				if not line:
					break
				logger.info('STDERR: {}'.format(line.rstrip()))
			# if there are still something
			for line in self.cmdobj.p.stderr:
				logger.info('STDERR: {}'.format(line.rstrip()))
			self.cmdobj.rc     = self.cmdobj.p.wait()
			self.cmdobj.stdout = self.cmdobj.p.stdout.read()
			self.done          = True
			logger.info('Return code: %s' % self.cmdobj.rc)
			logger.info('-' * 80)
		if quit and self.cmdobj.rc != 0:
			raise RuncmdException(self.cmdobj.cmd)
		return self

	@property
	def stdout(self):
		return self.run().cmdobj.stdout
	
	@property
	def stderr(self):
		return self.run().cmdobj.stderr

	@property
	def cmd(self):
		return self.cmdobj.cmd

	@property
	def returncode(self):
		return self.run().cmdobj.rc
	
	@property
	def rc(self):
		return self.returncode

	def pipe(self, tools = None, subcmd = False, dash = 'auto', equal = 'auto', duplistkey  = False, ignorefalse = True, base = None):
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
			prevcmd     = self.cmd
		)

class Shell(object):
	def __init__(self, 
		tools       = None,
		subcmd      = False,
		dash        = 'auto',
		equal       = 'auto',
		duplistkey  = False,
		ignorefalse = True,
		base        = None,
		stdin       = None,
		envs        = None,
		env         = None,
		cwd         = None,
		prevcmd     = None):
		self.dash          = dash
		self.equal         = equal
		self.duplistkey    = duplistkey
		self.ignorefalse   = ignorefalse
		self.subcmd        = subcmd
		self.tools         = TOOLS
		self.base          = base
		self.stdin         = stdin
		self.cwd           = cwd or getcwd()
		self.envs          = envs or env or {}
		self.prevcmd       = prevcmd
		self.targs         = {'': []}
		self._subcmdCalled = False
		self.tools.update(tools or {})

	def __call__(self, *args, **kwargs):
		self.targs.update(kwargs)
		self.targs[''] = list(args) + self.targs['']
		return self

	def _run(self, *args, **kwargs):
		kwargs[''] = list(args) + kwargs.get('', []) 
		cmd2run = '{}{} {} '.format(
			self.base, 
			'' if not self.targs else ' ' + cmdargs(
				self.targs, dash = self.dash, equal = self.equal, duplistkey = self.duplistkey, ignorefalse = self.ignorefalse),
			cmdargs(
				kwargs, dash = self.dash, equal = self.equal, duplistkey = self.duplistkey, ignorefalse = self.ignorefalse)
		)
		
		cmdobj = cmd.Cmd(cmd2run, stdin = self.stdin, shell = True, env = self.envs, cwd = self.cwd)
		
		if self.prevcmd:
			cmdobj.cmd = '{} | {}'.format(self.prevcmd, cmdobj.cmd)
		return ShellResult(cmdobj, self.tools)

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
		return Shell({}, subcmd = False, dash = self.dash, equal = self.equal, duplistkey = self.duplistkey, ignorefalse = self.ignorefalse, base = tool)

wc   = lambda *args, **kwargs: Shell().wc(*args, **kwargs).run()
wc_l = lambda filename: int(wc(l = filename).stdout.split()[0])
wcl  = wc_l

gzip      = lambda *args,  **kwargs: Shell().gzip(*args, **kwargs).run()
gunzip    = lambda *args,  **kwargs: Shell().gunzip(*args, **kwargs).run()
gzip_to   = lambda source, dest: gzip(_ = source, c = True, _stdout = dest)
gunzip_to = lambda source, dest: gunzip(_ = source, c = True, _stdout = dest)

mv    = lambda *args,  **kwargs: Shell().mv(*args, **kwargs).run()
cp    = lambda *args,  **kwargs: Shell().cp(*args, **kwargs).run()
rm    = lambda *args,  **kwargs: Shell().rm(*args, **kwargs).run()
rm_rf = lambda *args: rm(*args, r = True, f = True)
mkdir = lambda *args,  **kwargs: Shell().mkdir(*args, **kwargs).run()
ln    = lambda *args,  **kwargs: Shell().ln(*args, **kwargs).run()
ln_s  = lambda source, dest: ln(source, dest, s = True)

head  = lambda *args, **kwargs: Shell().head(*args, **kwargs).run()
tail  = lambda *args, **kwargs: Shell().tail(*args, **kwargs).run()
grep  = lambda *args, **kwargs: Shell().grep(*args, **kwargs).run()
sort  = lambda *args, **kwargs: Shell(dash = '-', equal = ' ').grep(*args, **kwargs).run()
uniq  = lambda *args, **kwargs: Shell().uniq(*args, **kwargs).run()
touch = lambda *args, **kwargs: Shell().touch(*args, **kwargs).run()

kill  = lambda *args, **kwargs: Shell().kill(*args, **kwargs).run()
kill9 = lambda *procids: kill(**{'9': True, '': procids})

zcat = lambda *args, **kwargs: Shell().zcat(*args, **kwargs).run()
cat  = lambda *args, **kwargs: Shell().cat(*args, **kwargs).run()
