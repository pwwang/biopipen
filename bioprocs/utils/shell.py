from . import cmdargs, runcmd
from functools import partial

class Shell(object):
	def __init__(self, tools = None, subcommand = False, dash = 'auto', equal = 'auto', base = None):
		self.dash       = dash
		self.equal      = equal
		self.subcommand = subcommand
		self.tools      = tools or {}
		self.base       = base

	def _run(self, *args, **kwargs):
		if '_' in kwargs and not '' in kwargs:
			kwargs[''] = kwargs['_']
			del kwargs['_']
		if args:
			kwargs[''] = kwargs.get('', []) + list(args)
		if '_stdout' in kwargs:
			outfile = kwargs['_stdout']
			del kwargs['_stdout']
		else:
			outfile = None

		cmd = '{} {} '.format(self.base, cmdargs(kwargs, dash = self.dash, equal = self.equal))
		if outfile:
			cmd = '{} > {!r}'.format(cmd, outfile)
		return runcmd(cmd, ret = 'stdout')

	def __getattr__(self, tool):
		if not self.subcommand:
			self.base = '{!r} {!r}'.format(self.base, tool) if self.base else '{!r}'.format(tool)
			return self._run
		return Shell({}, subcommand = False, dash = self.dash, equal = self.equal, base = tool)

wc        = lambda *args, **kwargs: Shell().wc(*args, **kwargs)
wc_l      = lambda filename: int(wc(l = filename).split()[0])
gzip      = lambda *args, **kwargs: Shell().gzip(*args, **kwargs)
gunzip    = lambda *args, **kwargs: Shell().gunzip(*args, **kwargs)
gzip_to   = lambda source, dest: gzip(_ = source, c = True, _stdout = dest)
gunzip_to = lambda source, dest: gunzip(_ = source, c = True, _stdout = dest)
mv        = lambda *args, **kwargs: Shell().mv(*args, **kwargs)
cp        = lambda *args, **kwargs: Shell().cp(*args, **kwargs)
head      = lambda *args, **kwargs: Shell().head(*args, **kwargs)
tail      = lambda *args, **kwargs: Shell().tail(*args, **kwargs)
grep      = lambda *args, **kwargs: Shell().grep(*args, **kwargs)
kill      = lambda *args, **kwargs: Shell().kill(*args, **kwargs)
kill9     = lambda *procids: kill(**{'9': True, '': procids})
