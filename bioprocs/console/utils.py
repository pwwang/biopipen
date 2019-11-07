"""Utilities for bioprocs script"""
import re
import sys
from pathlib import Path
from tempfile import gettempdir
from contextlib import contextmanager
from colorama import Back
from pyparam import Helps, HelpAssembler
from pyppl import utils, Channel, PyPPL, Box
import bioprocs

def substrReplace(string, starts, lengths, replace):
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

def highlight(origin, query, incase = True):
	"""Highlight string with query string"""
	# get all occurrences of q
	if incase:
		occurs = [m.start() for m in re.finditer(query.lower(), origin.lower())]
	else:
		occurs = [m.start() for m in re.finditer(query, origin)]
	lengths = [len(query)] * len(occurs)
	return substrReplace(origin, occurs, lengths, [
		'{}{}{}'.format(Back.LIGHTRED_EX, origin[occur:occur+length], Back.RESET)
		for occur, length in zip(occurs, lengths)])

def highlightMulti(line, queries):
	"""Highlight a string with multiple queries"""
	for query in queries:
		line = highlight(line, query)
	return line

def subtractDict(bigger, smaller, prefix = ''):
	"""Subtract a dict from another"""
	ret = bigger.copy()
	for key, val in smaller.items():
		if key not in ret:
			continue
		if isinstance(ret[key], dict) and isinstance(val, dict) and ret[key] != val:
			ret[key] = subtractDict(ret[key], val, prefix + '  ')
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
		self.desc   = self.module.__doc__ and self.module.__doc__.strip() or '[ Not documented. ]'
		self._procs = {}

	def procs(self):
		"""Get the processes of the module"""
		if self._procs:
			return self._procs
		for proc, factory in self.module._mkenvs.items():
			if len(proc) < 3 or proc[0] != '_' or proc[1] != 'p' \
				or not (proc[2].isdigit() or proc[2].isupper()):
				continue
			procobj = factory()
			self._procs[proc[1:]] = Process(procobj, self.module, factory.__doc__ and factory.__doc__.splitlines() or [])

		return self._procs

	def toHelps(self, helpsec):
		"""Send me to a Helps section"""
		helpsec.prefix = ''
		helpsec.add((self.name, '', self.desc))

	@contextmanager
	def loadProcs(self, index = None):
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

	def toHelpsAsSec(self, helps):
		"""Add me to helps as a section"""
		helps.add(self.name + ': ' + self.desc, sectype = 'option', prefix = '')

	def toHelpsWithProcs(self, helps, nmods):
		"""Add me to helps with procs information"""
		self.toHelpsAsSec(helps)
		with self.loadProcs("%s/%s" % (len(helps), nmods)) as procs:
			for proc in procs.values():
				proc.toHelps(helps.select(self.name + ': '))

	def addToCompletions(self, comp):
		"""Add me to completions"""
		comp.addCommand(self.name, self.desc)
		with self.loadProcs() as procs:
			for pname, proc in procs.items():
				comp.addCommand(self.name + '.' + pname, proc.desc[0])
				proc.addToCompletions(comp.command(self.name + '.' + pname))

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

	def toHelps(self, helpsec):
		"""Send me to a Helps section"""
		helpsec.prefix = ''
		if self.proc.origin == self.name:
			helpsec.add((self.name, '', self.desc))
		else:
			helpsec.add((self.name, '', 'Alias of: %s' % self.proc.origin))

	@staticmethod
	def _parseOptionSec(sec):
		ret      = {}
		sec      = Module.deindent(sec)
		lastdesc = []
		for line in sec:
			if not line:
				continue
			if line[0] in ('\t', ' '):
				lastdesc.append(line)
				continue
			parts = utils.split(line, ':', trim = False)
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
	def _parsePlainSec(sec):
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
				self._parsed[key] = Process._parseOptionSec(sec)
			else:
				self._parsed[key] = Process._parsePlainSec(sec)
		return self._parsed

	@staticmethod
	def defaultVal(val, prefix = 'Default: ', indent = 0):
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
		if isinstance(self.proc.config.input, dict):
			inkeys = self.proc.config.input.keys()
		elif isinstance(self.proc.config.input, str):
			inkeys = utils.split(self.proc.config.input, ',')
		else:
			inkeys = self.proc.config.input
		ret = utils.OBox()
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
			outs = utils.split(self.proc.config.output, ',')
		else:
			outs = self.proc.config.output
		ret = utils.OBox()
		for out in outs:
			parts = utils.split(out, ':')
			if len(parts) == 2:
				outname, outype, default = parts[0], 'var', parts[1]
			else:
				outname, outype, default = parts
			ret[outname] = (outype, default)
		return ret

	def addToCompletions(self, comp):
		"""Add me to completions"""
		for inname in self.inputs():
			docdesc = self.parsed().get('input', {}).get(inname, (None, ['Input %s' % inname]))[1]
			comp.addOption('-i.' + inname, docdesc[0])
		for outname, outypedeft in self.outputs().items():
			docdesc = self.parsed().get('output', {}).get(
				outname, ('var', Process.defaultVal(outypedeft[1])[0]))[1]
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
				outname, ('var', Process.defaultVal(outdeft)))
			outype = outype or doctype or 'var'

			if not docdesc or ('default: ' not in docdesc[-1].lower() and len(docdesc[-1]) > 20):
				docdesc.extend(Process.defaultVal(outdeft, indent = 5))
			elif 'default: ' not in docdesc[-1].lower():
				defaults = Process.defaultVal(outdeft, indent = len(docdesc[-1]) + 6)
				docdesc[-1] += ' ' + defaults.pop(0)
				docdesc.extend(defaults)
			self._helps.select('Output options').add(('-o.' + outname, '<%s>' % outype, docdesc))

		# args
		self._helps.add('Process arguments', sectype = 'option', prefix = '-')
		for key, val in self.proc.args.items():
			doctype, docdesc = self.parsed().get('args', {}).get(key, (None, ['[ Not documented. ]']))
			doctype = doctype or type(val).__name__
			if doctype == 'NoneType':
				doctype = 'auto'
			if not docdesc or ('default: ' not in docdesc[-1].lower() and len(docdesc[-1]) > 20):
				docdesc.extend(Process.defaultVal(val, indent = 5))
			elif 'default: ' not in docdesc[-1].lower():
				defaults = Process.defaultVal(val, indent = len(docdesc[-1]) + 6)
				docdesc[-1] += ' ' + defaults.pop(0)
				docdesc.extend(defaults)

			self._helps.select('Process arguments').add((
				'-args.' + key,
				'<%s>' % doctype if str(doctype).lower() != 'bool' else '[bool]',
				docdesc))

		self._helps.add('Other options', sectype = 'option', prefix = '-')
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

		return self._helps

	@staticmethod
	def _updateArgs(args):
		# replace ',' with '.' in key
		# a bug of python-box, copy lost metadata
		#ret = args.copy() # try to keep the type
		ret = Box()
		for key, val in args.items():
			ret[key.replace(',', '.')] = Process._updateArgs(val) \
				if isinstance(val, dict) else val
		return ret

	def printHelps(self, error = None, halt = True):
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

	def run(self, opts):
		"""Construct a pipeline with the process and run it"""
		if any(opts.get(h) for h in bioprocs.params._hopts) or \
			all(key in bioprocs.params._hopts for key in opts):
			self.printHelps()
		if not opts.get('i') or not isinstance(opts.i, dict):
			self.printHelps(error = 'No input specified.')

		indata = {}
		for inkey, intype in self.inputs().items():
			# We should allow some input options to be empty
			#if opts.i.get(inkey) is None:
			#	self.printHelps(error = 'Input "[-i.]%s" is not specified.' % inkey)

			# to be smart on "files"
			if intype in self.proc.IN_FILESTYPE and inkey in opts.i \
				and not isinstance(opts.i[inkey], list):
				indata[inkey + ':' + intype] = Channel.create([[opts.i[inkey]]])
			else:
				indata[inkey + ':' + intype] = Channel.create(opts.i.get(inkey))
		self.proc.input = indata

		if opts.get('o') is not None:
			if not isinstance(opts.o, dict):
				self.printHelps(error = 'Malformat output specification.')

			outdata = utils.OBox()
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
					outdata[outkey] = out.name
				else:
					outdata[outkey] = opts.o[outkey]
			self.proc.output = outdata

		if opts.get('args'):
			if not isinstance(opts.args, dict):
				self.printHelps(error = 'Malformat args specification.')
			self.proc.args.update(Process._updateArgs(opts.args))

		config = {
			'_log'  : {'file': None },
			'ppldir': Path(gettempdir()) / 'bioprocs.workdir'
		}
		config.update(Process._updateArgs(opts.get('config', {})))

		for key, val in opts.items():
			if key in bioprocs.params._hopts + ['i', 'o', 'args', 'config']:
				continue
			setattr(self.proc, key, val)

		PyPPL(config).start(self.proc).run()


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
		self.desc = self.module.__doc__ and self.module.__doc__.strip() or '[ Not documented. ]'

	def run(self):
		"""Run the pipeline"""
		prog = 'bioprocs ' + self.name
		sys.argv = [prog] + sys.argv[2:]
		bioprocs.params._prog = prog
		bioprocs.params._assembler.progname = prog

		self.module.main()

	def addToCompletions(self, comp):
		"""Add me to completions"""
		comp.addCommand(self.name, self.desc)
		self.module.params._addToCompletions(
			comp.command(self.name), withtype = False, alias = True)
