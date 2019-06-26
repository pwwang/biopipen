
import re
import sys
from pathlib import Path
from tempfile import gettempdir
from contextlib import contextmanager
from colorama import Back
from pyparam import Helps, HelpAssembler
from pyppl import utils, Channel, PyPPL
import bioprocs

def substrReplace(string, starts, lengths, replace):
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
	for query in queries:
		line = highlight(line, query)
	return line

def subtractDict(bigger, smaller, prefix = ''):
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

	@staticmethod
	def modules():
		return [module.stem
			for module in Path(bioprocs.__file__).parent.glob('*.py')
			if not module.stem.startswith('_')]

	@staticmethod
	def deindent(lines):
		indention = lines[0][:-len(lines[0].lstrip())]
		ret = []
		for line in lines:
			if not line:
				continue
			if not line.startswith(indention):
				raise ValueError('Unexpected indention at doc line:\n' + repr(line))
			ret.append(line[len(indention):].replace('\t', '  '))
		return ret

	@staticmethod
	def procnameFromDoc(docstr):
		doclines       = docstr.splitlines()
		doclines       = Module.deindent(doclines)
		nameindex      = []
		descindex      = []
		nameindex_flag = descindex_flag = False
		for i, line in enumerate(doclines):
			if line.startswith('@'):
				nameindex_flag = descindex_flag = False
			if line.startswith('@name'):
				nameindex_flag = True
			elif line.startswith('@description'):
				descindex_flag = True
			elif nameindex_flag:
				nameindex.append(i)
			elif descindex_flag:
				descindex.append(i)
		if len(nameindex) != 1:
			raise ValueError('Proc name should be on a separated line.')
		name = doclines[nameindex[0]].strip()
		desc = Module.deindent([doclines[i] for i in descindex])
		doc  = doclines[(max(nameindex + descindex)+1):]
		return name, desc, doc

	@staticmethod
	def scanAlias(proc, modsrc):
		"""Document of proc not found in this module
		So try to find it in the contexts of
		'from bedtools import pBedIntersect' or
		'from bedtools import *'
		"""
		# scan 'from bedtools import pBedIntersect' first
		match = re.search(r'^from\s+(?:bioprocs)?\.(\w+?)\s+import\s+%s$' % proc, modsrc, re.M)
		if not match:
			match = re.search(r'^from\s+(?:bioprocs)?\.(\w+?)\s+import\s+\*$', modsrc, re.M)
		if not match:
			return None
		module = Module(match.group(1))
		aliasprocs = module.procs()
		if proc in aliasprocs:
			return aliasprocs[proc].desc, aliasprocs[proc].doc
		return None

	def __init__(self, name):
		self.name   = name
		try:
			self.module = getattr(__import__('bioprocs', fromlist = [name]), name)
		except AttributeError:
			raise AttributeError('No such module: %s' % name) from None
		self.desc   = self.module.__doc__ and self.module.__doc__.strip() or '[ Not documented. ]'
		self._procs = {}

	def procs(self):
		if self._procs:
			return self._procs

		# load docs
		modulefile = Path(self.module.__file__)
		modulefile = modulefile.with_suffix('.py')
		modulesrc  = modulefile.read_text()
		regex      = re.compile(
			r'(\"\"\"|\'\'\')\s*\n(\s*@name:\s+p[_A-Z0-9].+\s+@[\s\S]+?)\1', re.M)
		docstrs    = regex.finditer(modulesrc)
		docs       = {}
		for docstr in docstrs:
			name, desc, doc = Module.procnameFromDoc(docstr.group(2))
			docs[name] = (desc, doc)

		self._procs = {}
		for pname in dir(self.module):
			if len(pname) < 2 or pname[0] != 'p' or (
				pname[1] != '_' and not pname[1].isdigit() and not pname[1].isupper()):
				continue
			if pname in docs:
				desc, doc = docs[pname]
			else:
				descdoc = Module.scanAlias(pname, modulesrc)
				desc, doc = descdoc if descdoc else (['[ Not documented. ]'], [])
			self._procs[pname] = Process(pname, self.module, desc, doc)
		return self._procs

	def toHelps(self, helpsec):
		helpsec.prefix = ''
		helpsec.add((self.name, '', self.desc))

	@contextmanager
	def loadProcs(self, index = None):
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
		helps.add(self.name + ': ' + self.desc, sectype = 'option', prefix = '')

	def toHelpsWithProcs(self, helps, nmods):
		self.toHelpsAsSec(helps)
		with self.loadProcs("%s/%s" % (len(helps), nmods)) as procs:
			for name, proc in procs.items():
				proc.toHelps(helps.select(self.name + ': '))

	def addToCompletions(self, comp):
		comp.addCommand(self.name, self.desc)
		with self.loadProcs() as procs:
			for pname, proc in procs.items():
				comp.addCommand(self.name + '.' + pname, proc.desc[0])
				proc.addToCompletions(comp.command(self.name + '.' + pname))

class Process:

	def __init__(self, proc, module, desc, doc):
		self.name    = proc
		self.module  = module
		self.proc    = getattr(module, proc)
		self.desc    = desc
		self.doc     = doc
		self._helps  = Helps()
		self._parsed = {}

	def toHelps(self, helpsec):
		helpsec.prefix = ''
		helpsec.add((self.name, '', self.desc))

	@staticmethod
	def _parseSec(sec):
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
				name, typ = nametype.split(' (')
			elif ':' in nametype:
				name, typ = nametype.split(':', 1)
			else:
				name, typ = nametype.strip(), ''
			ret[name.strip()] = (typ.strip(), lastdesc)
		return ret

	def parsed(self):
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
			self._parsed[key] = Process._parseSec(sec)
		return self._parsed

	@staticmethod
	def defaultVal(val):
		return str(val) if val != '' else "''"

	def inputs(self):
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
		for inname in self.inputs():
			docdesc = self.parsed().get('input', {}).get(inname, (None, ['Input %s' % inname]))[1]
			comp.addOption('-i.' + inname, docdesc[0])
		for outname, outypedeft in self.outputs().items():
			docdesc = self.parsed().get('output', {}).get(
				outname, ('var', ['Default: ' + Process.defaultVal(outypedeft[1])]))[1]
			comp.addOption('-o.' + outname, docdesc[0])
		for key, val in self.proc.args.items():
			docdesc = self.parsed().get('args', {}).get(key, ('auto', ['[ Not documented. ]']))[1]
			comp.addOption('-args.' + key, docdesc[0])
		# add -config.
		comp.addOption('-config.', 'Pipeline configrations, such as -config._log.file')

	def helps(self):
		"""Construct help page using doc"""
		if self._helps:
			return self._helps

		self._helps.add('Name', Path(self.module.__file__).stem + '.' + self.name)
		self._helps.add('Description', self.desc)

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
				outname, ('var', ['Default: ' + Process.defaultVal(outdeft)]))
			outype = outype or doctype or 'var'

			if not docdesc or ('default: ' not in docdesc[-1].lower() and len(docdesc[-1]) > 20):
				docdesc.append('Default: ' + Process.defaultVal(outdeft))
			elif 'default: ' not in docdesc[-1].lower():
				docdesc[-1] += ' Default: ' + Process.defaultVal(outdeft)
			self._helps.select('Output options').add(('-o.' + outname, '<%s>' % outype, docdesc))

		# args
		self._helps.add('Process arguments', sectype = 'option', prefix = '-')
		for key, val in self.proc.args.items():
			doctype, docdesc = self.parsed().get('args', {}).get(key, ('auto', ['[ Not documented. ]']))
			doctype = doctype or 'auto'
			if not docdesc or ('default: ' not in docdesc[-1].lower() and len(docdesc[-1]) > 20):
				docdesc.append('Default: ' + Process.defaultVal(val))
			elif 'default: ' not in docdesc[-1].lower():
				docdesc[-1] += ' Default: ' + Process.defaultVal(val)

			self._helps.select('Process arguments').add(('-args.' + key, '<%s>' % doctype, docdesc))

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
		ret = args.copy() # try to keep the type
		ret.clear()
		for key, val in args.items():
			ret[key.replace(',', '.')] = Process._updateArgs(val) \
				if isinstance(val, dict) else val
		return ret

	def printHelps(self, error = None, halt = True):
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
		if any(opts.get(h) for h in bioprocs.params._hopts) or \
			all(key in bioprocs.params._hopts for key in opts):
			self.printHelps()
		if not opts.get('i') or not isinstance(opts.i, dict):
			self.printHelps(error = 'No input specified.')

		indata = {}
		for inkey, intype in self.inputs().items():
			if opts.i.get(inkey) is None:
				self.printHelps(error = 'Input "[-i.]%s" is not specified.' % inkey)
			indata[inkey + ':' + intype] = Channel.create(opts.i[inkey])
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
				out = Path(opts.o[outkey])
				if out.name != opts.o[outkey]: # we have path
					if self.proc.exdir and out.parent != self.proc.exdir:
						raise ValueError('Cannot have output files/dirs with different parents.')
					self.proc.exdir = str(out.parent)
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
		return [pplfile.stem[9:]
			for pplfile in (
				Path(bioprocs.__file__).parent / 'bin' / 'pipeline').glob('bioprocs_*.py')
			if not pplfile.stem.startswith('_')]

	def __init__(self, name):
		self.name = name
		self.module = __import__('bioprocs.bin.pipeline', fromlist = ['bioprocs_' + name])
		self.module = getattr(self.module, 'bioprocs_' + name)
		self.desc = self.module.__doc__ and self.module.__doc__.strip() or '[ Not documented. ]'

	def run(self):
		prog = 'bioprocs ' + self.name
		sys.argv = [prog] + sys.argv[2:]
		bioprocs.params._prog = prog
		bioprocs.params._assembler.progname = prog

		self.module.main()

	def addToCompletions(self, comp):
		comp.addCommand(self.name, self.desc)
		self.module.params._addToCompletions(
			comp.command(self.name), withtype = False, alias = True)
