"""Utilities for bioprocs script"""
import re
import sys
import yaml
import cmdy
from pathlib import Path
from tempfile import gettempdir
from contextlib import contextmanager
from colorama import Back
from pyparam import Helps, HelpAssembler
from diot import Diot, OrderedDiot
from pyppl import utils, Channel, PyPPL, Proc
from pyppl.template import TemplateLiquid
from pyppl._proc import IN_FILESTYPE
from liquid.stream import LiquidStream
import bioprocs

def _split(string, delimit, trim = True):
	return LiquidStream.from_string(string).split(delimit, trim = trim)

def substr_replace(string, starts, lengths, replace):
	"""Replace substrings"""
	if not isinstance(starts, (tuple, list)):
		starts = [starts]
	if not isinstance(lengths, (tuple, list)):
		lengths = [lengths]
	assert len(starts) == len(lengths)
	if not isinstance(replace, (tuple, list)):
		replace = [replace] * len(starts)

	delta = 0
	for i, start in enumerate(starts):
		# adjust starts
		string = string[:start + delta] + replace[i] + string[start + lengths[i] + delta:]
		delta += len(replace[i]) - lengths[i]
	return string

def highlight(origin, query, incase = True, hicolor = Back.LIGHTRED_EX):
	"""Highlight string with query string"""
	# get all occurrences of q
	if incase:
		occurs = [m.start() for m in re.finditer(query.lower(), origin.lower())]
	else:
		occurs = [m.start() for m in re.finditer(query, origin)]
	lengths = [len(query)] * len(occurs)
	return substr_replace(origin, occurs, lengths, [
		'{}{}{}'.format(hicolor, origin[occur:occur+length], Back.RESET)
		for occur, length in zip(occurs, lengths)])

def highlight_multi(line, queries):
	"""Highlight a string with multiple queries"""
	for query in queries:
		line = highlight(line, query)
	return line

def subtract_dict(bigger, smaller, prefix = ''):
	"""Subtract a dict from another"""
	ret = bigger.copy()
	for key, val in smaller.items():
		if key not in ret:
			continue
		if isinstance(ret[key], dict) and isinstance(val, dict) and ret[key] != val:
			ret[key] = subtract_dict(ret[key], val, prefix + '  ')
		elif ret[key] == val:
			del ret[key]
	return ret

class Module:
	"""A module of bioprocs"""
	@staticmethod
	def modules():
		"""Get all modules"""
		return [module.stem
			for module in Path(bioprocs.__file__).parent.glob('*.py')
			if not module.stem.startswith('_')]

	@staticmethod
	def deindent(lines):
		"""Remove indent based on the first line"""
		indention = ''
		for line in lines:
			if not line:
				continue
			indention = line[:-len(line.lstrip())]
			break
		ret = []
		for line in lines:
			if not line:
				continue
			if not line.startswith(indention):
				raise ValueError('Unexpected indention at doc line:\n' + repr(line))
			ret.append(line[len(indention):].replace('\t', '  '))
		return ret

	def __init__(self, name):
		self.name   = name
		try:
			self.module = getattr(__import__('bioprocs', fromlist = [name]), name)
		except AttributeError as ex:
			if "has no attribute '%s'" % name in str(ex):
				raise AttributeError('No such module: %s' % name) from None
			else:
				raise
		self.desc   = self.module.__doc__ and self.module.__doc__.strip() or '[ Not documented ]'
		self._procs = {}

	def procs(self):
		"""Get the processes of the module"""
		if self._procs:
			return self._procs
		for attr in dir(self.module):
			proc = getattr(self.module, attr)
			if not isinstance(proc, Proc):
				continue
			doc = proc.config.long.splitlines()
			self._procs[attr] = Process(proc, self.module, doc)

		return self._procs

	def to_helps(self, helpsec):
		"""Send me to a Helps section"""
		helpsec.prefix = ''
		helpsec.add((self.name, '', self.desc))

	@contextmanager
	def load_procs(self, index = None):
		"""Loading procs with indicators"""
		info = 'Collecting module%s: %s ...' % (index and (' ' + index) or '', self.name)
		print('\r' + info, end = '')
		yield self.procs()
		print(end = '\r')
		print(' ' * len(info), end = '')
		if index:
			idx, total = index.split('/')
			if idx == total:
				print(end = '\r')

	def to_helps_as_sec(self, helps):
		"""Add me to helps as a section"""
		helps.add(self.name + ': ' + self.desc, sectype = 'option', prefix = '')

	def to_helps_with_procs(self, helps, nmods):
		"""Add me to helps with procs information"""
		self.to_helps_as_sec(helps)
		with self.load_procs("%s/%s" % (len(helps), nmods)) as procs:
			for proc in procs.values():
				proc.to_helps(helps.select(self.name + ': '))

	def add_to_completions(self, comp):
		"""Add me to completions"""
		comp.addCommand(self.name, self.desc)
		with self.load_procs() as procs:
			for pname, proc in procs.items():
				comp.addCommand(self.name + '.' + pname, proc.desc[0])
				proc.add_to_completions(comp.command(self.name + '.' + pname))

class Process:
	"""A bioprocs process"""
	def __init__(self, proc, module, doc):
		self.name    = proc.id
		self.module  = module
		self.proc    = proc
		self.desc    = self.proc.desc
		self.doc     = Module.deindent(doc) # list
		self._helps  = Helps()
		self._parsed = {}

	def to_helps(self, helpsec):
		"""Send me to a Helps section"""
		helpsec.prefix = ''
		if self.proc.origin == self.name:
			helpsec.add((self.name, '', self.desc))
		else:
			helpsec.add((self.name, '', 'Alias of: %s' % self.proc.origin))

	@staticmethod
	def _parse_option_sec(sec):
		ret      = {}
		sec      = Module.deindent(sec)
		lastdesc = []
		for line in sec:
			if not line:
				continue
			if line[0] in ('\t', ' '):
				lastdesc.append(line)
				continue
			parts = _split(line, ':', trim = False)
			nametype = parts.pop(0)
			lastdesc = [':'.join(parts).strip()]
			nametype = nametype.strip('`) ')
			if ' (' in nametype:
				name, typ = nametype.split(' (', 1)
			elif ':' in nametype:
				name, typ = nametype.split(':', 1)
			else:
				name, typ = nametype.strip(), ''
			ret[name.strip()] = (typ.strip(), lastdesc)
		return ret

	@staticmethod
	def _parse_plain_sec(sec):
		return Module.deindent(sec)

	def parsed(self):
		"""Get parsed doc"""
		if self._parsed:
			return self._parsed
		if not self.doc:
			self._parsed = {}
			return self._parsed
		self._parsed = {}
		lastsec = []
		for line in self.doc:
			if line.startswith('@'):
				secname = line.strip('@: ')
				lastsec = []
				self._parsed[secname] = lastsec
			else:
				lastsec.append(line)
		for key, sec in self._parsed.items():
			if key in ('input', 'output', 'args'):
				self._parsed[key] = Process._parse_option_sec(sec)
			else:
				self._parsed[key] = Process._parse_plain_sec(sec)
		return self._parsed

	@staticmethod
	def default_val(val, prefix = 'Default: ', indent = 0):
		"""Formatted value in help"""
		ret = []
		lines = utils.formatDict(val if val != '' else "''", 0).splitlines()
		for i, line in enumerate(lines):
			if i == 0:
				line = prefix + line
				if len(lines) > 1:
					line += ' \\'
				ret.append(line)
			else:
				ret.append(' ' * indent + line + ' \\')
		return ret

	def inputs(self):
		"""Get the input keys and types by definitions"""
		if isinstance(self.proc._input, dict):
			inkeys = self.proc._input.keys()
		elif isinstance(self.proc._input, str):
			inkeys = _split(self.proc._input, ',')
		else:
			inkeys = self.proc._input
		ret = OrderedDiot()
		for inkey in inkeys:
			if ':' not in inkey:
				inkey += ':var'
			inname, intype = inkey.split(':', 1)
			ret[inname] = intype
		return ret

	def outputs(self):
		"""Get the input keys, types and default values by definitions"""
		if isinstance(self.proc.config.output, dict):
			outs = list(self.proc.config.output.keys())
		elif isinstance(self.proc.config.output, str):
			outs = _split(self.proc.config.output, ',')
		else:
			outs = self.proc.config.output
		ret = OrderedDiot()
		for out in outs:
			parts = _split(out, ':')
			if len(parts) == 2:
				outname, outype, default = parts[0], 'var', parts[1]
			else:
				outname, outype, default = parts
			ret[outname] = (outype, default)
		return ret

	def requirements(self):
		req = self.parsed().get('requires')
		if not req:
			return '[ Not documented ]'
		data4render = { 'proc': self.proc,
						'args': self.proc.args}
		requires = yaml.safe_load(TemplateLiquid('\n'.join(req)).render(data4render))
		return yaml.dump(requires, default_flow_style = False)

	def add_to_completions(self, comp):
		"""Add me to completions"""
		for inname in self.inputs():
			docdesc = self.parsed().get('input', {}).get(inname, (None, ['Input %s' % inname]))[1]
			comp.addOption('-i.' + inname, docdesc[0])
		for outname, outypedeft in self.outputs().items():
			docdesc = self.parsed().get('output', {}).get(
				outname, ('var', Process.default_val(outypedeft[1])[0]))[1]
			comp.addOption('-o.' + outname, docdesc[0])
		for key in self.proc.args:
			docdesc = self.parsed().get('args', {}).get(key, ('auto', ['[ Not documented. ]']))[1]
			comp.addOption('-args.' + key, docdesc[0])
		# add -config.
		comp.addOption('-config.', 'Pipeline configrations, such as -config._log.file')

	def helps(self):
		"""Construct help page using doc"""
		if self._helps:
			return self._helps

		self._helps.add('Name',
			Path(self.module.__file__).stem + '.' + self.name + '%s [lang = %s]' % (
				'(%s)' % self.proc.origin if self.name != self.proc.origin else '',
				self.proc.lang))
		self._helps.add('Description', self.parsed().get('description') or self.desc)

		# input
		self._helps.add('Input options (Use \':list\' for multi-jobs)',
			sectype = 'option', prefix = '-')
		for inname, intype in self.inputs().items():
			doctype, docdesc = self.parsed().get('input', {}).get(inname, ('var', 'Input %s' % inname))
			intype = intype or doctype or 'var'
			self._helps.select('Input options').add(('-i.' + inname, '<%s>' % intype, docdesc))

		# output
		self._helps.add('Output options (\'exdir\' implied if path specified)',
			sectype = 'option', prefix = '-')

		for outname, outypedeft in self.outputs().items():
			outype, outdeft = outypedeft
			doctype, docdesc = self.parsed().get('output', {}).get(
				outname, ('var', Process.default_val(outdeft)))
			outype = outype or doctype or 'var'

			if not docdesc or ('default: ' not in docdesc[-1].lower() and len(docdesc[-1]) > 20):
				docdesc.extend(Process.default_val(outdeft, indent = 5))
			elif 'default: ' not in docdesc[-1].lower():
				defaults = Process.default_val(outdeft, indent = len(docdesc[-1]) + 6)
				docdesc[-1] += ' ' + defaults.pop(0)
				docdesc.extend(defaults)
			self._helps.select('Output options').add(('-o.' + outname, '<%s>' % outype, docdesc))

		# args
		self._helps.add('Process arguments', sectype = 'option', prefix = '-')
		for key, val in self.proc.args.items():
			doctype, docdesc = self.parsed().get('args', {}).get(key, (None, ['[ Not documented ]']))
			doctype = doctype or type(val).__name__
			if doctype == 'NoneType':
				doctype = 'auto'
			if not docdesc or ('default: ' not in docdesc[-1].lower() and len(docdesc[-1]) > 20):
				docdesc.extend(Process.default_val(val, indent = 5))
			elif 'default: ' not in docdesc[-1].lower():
				defaults = Process.default_val(val, indent = len(docdesc[-1]) + 6)
				docdesc[-1] += ' ' + defaults.pop(0)
				docdesc.extend(defaults)

			self._helps.select('Process arguments').add((
				'-args.' + key,
				'<%s>' % doctype if str(doctype).lower() != 'bool' else '[bool]',
				docdesc))

		self._helps.add('Other options', sectype = 'option', prefix = '-')
		# install requirements
		self._helps.select('Other options').add(
			('-install', '[BINDIR]', 'Install process requirements to BINDIR. Default: /usr/bin'))
		# validate requirements
		self._helps.select('Other options').add(
			('-validate', '', 'Validate process requirements.'))

		# process properties
		self._helps.select('Other options').add(
			('-<prop>', '', 'Process properties, such as -forks, -exdir, -cache ...'))

		# pipeline configurations
		self._helps.select('Other options').add(
			('-config.<subconf>[.subconf]', '',
			 'Pipeline configrations, such as -config._log.file'))

		# help
		self._helps.select('Other options').addParam(
			bioprocs.params[bioprocs.params._hopts[0]], bioprocs.params._hopts, ishelp = True)

		self._helps.add('Requirements', self.requirements())
		return self._helps

	@staticmethod
	def _update_args(args):
		# replace ',' with '.' in key
		# a bug of python-box, copy lost metadata
		#ret = args.copy() # try to keep the type
		ret = Diot()
		for key, val in args.items():
			ret[key.replace(',', '.')] = Process._update_args(val) \
				if isinstance(val, dict) else val
		return ret

	def print_helps(self, error = None, halt = True):
		"""Print helps"""
		assembler = HelpAssembler()
		error = error or []
		if isinstance(error, str):
			error = [error]
		for err in error:
			print(assembler.error(err))
		print('\n'.join(assembler.assemble(self.helps())), end = '')
		if halt:
			sys.exit(1)

	def _logger(self, msg, hiword = '', hicolor = Back.GREEN, end = '\n', prefix = True):
		if prefix:
			modname = Path(self.module.__file__).stem
			prefix_str = f'[{modname}.{self.name}] '
			prefix_str = highlight(prefix_str, prefix_str[1:-2],
				incase = False, hicolor = Back.MAGENTA)
		else:
			prefix_str = ''
		if not hiword:
			print(f'{prefix_str}{msg}', end = end)
		else:
			himsg = highlight(msg, hiword, hicolor = hicolor)
			print(f'{prefix_str}{himsg}', end = end)

	def _run_validate(self):
		req = self.requirements()
		if req == '[ Not documented ]':
			raise ValueError("Requirements for process %r is not documented, "
				"I have no idea about the requirements." % self.proc.id)
		requires = yaml.safe_load(req)
		failed = []
		for tool, info in requires.items():
			self._logger(f"Validating {tool} ... ", tool, end = '', hicolor = Back.YELLOW)
			cmd = cmdy.bash(c = info['validate'], _raise = False)
			if cmd.rc != 0:
				failed.append((tool, info))
				self._logger("Failed.", "Failed", hicolor = Back.RED, prefix = False)
				self._logger("Validation command: " + cmd.cmd)
				for err in cmd.stderr.splitlines():
					self._logger(f"  {err}")
			else:
				self._logger("Installed.", "Installed", hicolor = Back.GREEN, prefix = False)
		return failed

	def _run_install(self, failed, bindir):
		if not failed:
			self._logger('All requirements met, nothing to install.')
		else:
			for tool, info in failed:
				self._logger('Installing %s ...' % tool, tool, hicolor = Back.YELLOW)
				cmd = cmdy.bash(c = info['install'].replace('$bindir$', bindir),
					_iter = 'err', _raise = True)
				for line in cmd.__iter__():
					self._logger(f'  {line}'.rstrip())
				if cmd.rc != 0:
					self._logger("Failed to install, please intall it manually.", "Failed",
						hicolor = Back.RED)
					self._logger("  " + cmd.cmd)
				else:
					self._logger("Succeeded!", "succeeded")
			self._run_validate()

	def run(self, opts):
		"""Construct a pipeline with the process and run it"""
		if any(opts.get(h) for h in bioprocs.params._hopts) or \
			all(key in bioprocs.params._hopts for key in opts):
			self.print_helps()

		# install requirements
		if 'install' in opts or 'validate' in opts:
			failed = self._run_validate()
			if 'install' in opts:
				self._run_install(failed, '/usr/bin' if opts.install is True else opts.install)
			sys.exit(0)

		if not opts.get('i') or not isinstance(opts.i, dict):
			self.print_helps(error = 'No input specified.')

		indata = OrderedDiot()
		for inkey, intype in self.inputs().items():
			# We should allow some input options to be empty
			#if opts.i.get(inkey) is None:
			#	self.print_helps(error = 'Input "[-i.]%s" is not specified.' % inkey)

			# to be smart on "files"
			if intype in IN_FILESTYPE and inkey in opts.i \
				and not isinstance(opts.i[inkey], list):
				indata[inkey + ':' + intype] = Channel.create([[opts.i[inkey]]])
			else:
				indata[inkey + ':' + intype] = Channel.create(opts.i.get(inkey))
		self.proc.input = indata

		if opts.get('o') is not None:
			if not isinstance(opts.o, dict):
				self.print_helps(error = 'Malformat output specification.')

			outdata = OrderedDiot()
			for outkey, outypedeft in self.outputs().items():
				outype, outdeft = outypedeft
				if not opts.o.get(outkey):
					outdata[outkey + ':' + outype] = outdeft
					continue
				if outype not in ('file', 'dir'):
					outdata[outkey + ':' + outype] = opts.o[outkey]
					continue
				# try to extract exdir from output
				if '/' in opts.o[outkey]:
					out = Path(opts.o[outkey])
					if self.proc.exdir and out.parent != self.proc.exdir:
						raise ValueError('Cannot have output files/dirs with different parents as exdir.')
					self.proc.exdir = str(out.parent)
					outdata[outkey + ':' + outype] = out.name
				else:
					outdata[outkey + ':' + outype] = opts.o[outkey]
			self.proc.output = outdata

		args = opts.pop('args', {})
		if not isinstance(args, dict):
			self.print_helps(error = 'Malformat args specification.')
		self.proc.args.update(Process._update_args(args))

		self.proc.config.update(opts.pop('config', {}))
		# config = {
		# 	'logger'  : {'file': None },
		# 	'ppldir': Path(gettempdir()) / 'bioprocs.workdir'
		# }
		# config.update(Process._update_args(opts.get('config', {})))

		for key, val in opts.items():
			if key in bioprocs.params._hopts + ['i', 'o']:
				continue

			setattr(self.proc, key, val)

		PyPPL(logger_file = False,
			ppldir = Path(gettempdir()) / 'bioprocs.workdir').start(self.proc).run()


class Pipeline:
	"""Assembled pipeline"""

	@staticmethod
	def pipelines():
		"""Get all available pipelines"""
		return [pplfile.stem[9:]
			for pplfile in (
				Path(bioprocs.__file__).parent / 'console' / 'pipeline').glob('bioprocs_*.py')
			if not pplfile.stem.startswith('_')]

	def __init__(self, name):
		self.name = name
		self.module = __import__('bioprocs.console.pipeline', fromlist = ['bioprocs_' + name])
		self.module = getattr(self.module, 'bioprocs_' + name)
		self.desc = self.module.__doc__ and self.module.__doc__.strip() or '[ Not documented ]'

	def run(self):
		"""Run the pipeline"""
		prog = 'bioprocs ' + self.name
		sys.argv = [prog] + sys.argv[2:]
		bioprocs.params._prog = prog
		bioprocs.params._assembler.progname = prog

		self.module.main()

	def add_to_completions(self, comp):
		"""Add me to completions"""
		comp.addCommand(self.name, self.desc)
		self.module.params._addToCompletions(
			comp.command(self.name), withtype = False, alias = True)
