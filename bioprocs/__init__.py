__version__ = '0.1.2'

import inspect
from pathlib import Path
from tempfile import gettempdir
from sys import executable, modules
from pyppl import Proc
from pyparam import Params
from . import _envs

# open to R (reticulate) to get the path of r util scripts
HERE    = Path(__file__).resolve().parent

params = Params()
params._loadFile(HERE / 'params.toml')

params.cachedir.value = Path(params.cachedir.value).expanduser().as_posix()
params.tmpdir.value = gettempdir()	if Path(params.tmpdir.value) == Path('/tmp') \
									else params.tmpdir.value

cfgfiles = [
	Path.home() / '.bioprocs.config', # values overwritten
	Path.home() / '.bioprocs.json',
	Path.home() / '.bioprocs.toml',
	Path('.') / '.bioprocs.config',
	Path('.') / '.bioprocs.json',
	Path('.') / '.bioprocs.toml',
	'bioprocs.osenv'
]
for cfgfile in cfgfiles:
	if isinstance(cfgfile, Path) and not cfgfile.exists():
		continue
	params._loadFile (cfgfile)

# lock the params in case the options are overwritten unintentionally.
params._locked = True

cachedir = Path(params.cachedir.value)
if not cachedir.exists():
	cachedir.mkdir()

EXT_MAP = {
	'Rscript': 'R',
	'python' : 'py',
	'python2': 'py',
	'python3': 'py',
}

def delefactory():
	"""The factory to give the delegator for modkit"""
	frame  = inspect.currentframe().f_back
	module = inspect.getmodule(frame)
	def delegator(proc):
		try:
			procfac =  module._mkenvs['_' + proc]
			if not procfac:
				raise KeyError
		except KeyError as exc:
			raise ImportError('No such process: {!r}'.format(proc)) from exc
		if not callable(procfac):
			raise ImportError('Wrong type of process factory: {!r} in module {!r}'.format(
				'_' + proc, module.__name__))
		return procfac()
	return delegator

def _procfactory(procfunc, pid, alias, mdname, doc):
	caller = inspect.getframeinfo(inspect.stack()[2][0])
	callerdir = Path(caller.filename).parent.resolve()
	def findscript(script):
		if not script or not script.startswith ('file:'):
			return script
		scriptfile = Path(script[5:])
		if scriptfile.is_absolute():
			return script
		scriptfile = callerdir.joinpath(scriptfile)
		return 'file:{}'.format(scriptfile)

	def factory():
		proc = procfunc()
		if isinstance(proc, dict):
			script = report = None
			if 'script' in proc:
				script = proc.pop('script')
			if 'report' in proc:
				report = proc.pop('report')
			proc = Proc(**proc)
			if script:
				proc.config.script = script
			if report:
				proc.config.report = report
		proc.id = alias
		proc.props.origin = pid
		lang = Path(proc.lang).name
		ext  = '.' + EXT_MAP.get(lang, lang)
		proc.config.script = proc.config.script or 'file:scripts/{}/{}{}'.format(mdname, pid, ext)
		proc.script = findscript(proc.config.script)
		report = proc.config.report or 'file:reports/{}/{}.md'.format(mdname, pid)
		report = findscript(report)
		if Path(report[5:]).is_file():
			proc.report = report
		return proc
	factory.__doc__ = doc
	return factory

def procfactory(procfunc):
	mdname = procfunc.__module__.split('.')[-1]
	pid    = procfunc.__name__.lstrip('_')
	module = modules[procfunc.__module__]
	args   = inspect.signature(procfunc).parameters
	alias  = args.get('alias')
	if alias:
		alias = alias.default
		if alias[0] == '_':
			alias = alias[1:]
		module._mkenvs['_' + alias] = _procfactory(procfunc, pid, alias, mdname, procfunc.__doc__)
	return _procfactory(procfunc, pid, pid, mdname, procfunc.__doc__)
