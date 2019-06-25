
import re
import sys
import math
from pathlib import Path
from contextlib import contextmanager
from colorama import Back
from pyparam import Helps, HelpAssembler, params
from pyppl import utils
import bioprocs

def substrReplace(s, starts, lengths, replace):
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
		s = s[:start + delta] + replace[i] + s[start + lengths[i] + delta:]
		delta += len(replace[i]) - lengths[i]
	return s

def highlight(origin, q, incase = True):
	# get all occurrences of q
	if incase:
		occurs = [m.start() for m in re.finditer(q.lower(), origin.lower())]
	else:
		occurs = [m.start() for m in re.finditer(q, origin)]
	lengths = [len(q)] * len(occurs)
	return substrReplace(origin, occurs, lengths, [
		'{}{}{}'.format(Back.LIGHTRED_EX, origin[occur:occur+length], Back.RESET)
		for occur, length in zip(occurs, lengths)])

def highlightMulti(line, queries):
	for query in queries:
		line = highlight(line, query)
	return line

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
			if doclines[i].startswith('@'):
				nameindex_flag = descindex_flag = False
			if doclines[i].startswith('@name'):
				nameindex_flag = True
			elif doclines[i].startswith('@description'):
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
		m = re.search(r'^from\s+(?:bioprocs)?\.(\w+?)\s+import\s+%s$' % proc, modsrc, re.M)
		if not m:
			m = re.search(r'^from\s+(?:bioprocs)?\.(\w+?)\s+import\s+\*$', modsrc, re.M)
		if not m:
			return None
		module = Module(m.group(1))
		aliasprocs = module.procs()
		if proc in aliasprocs:
			return aliasprocs[proc].desc, aliasprocs[proc].doc
		return None

	def __init__(self, name):
		self.name   = name
		self.module = getattr(__import__('bioprocs', fromlist = [name]), name)
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
				desc, doc = descdoc if descdoc else ('[ Not documented. ]', [])
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

class Process:

	def __init__(self, proc, module, desc, doc):
		self.name   = proc
		self.module = module
		self.proc   = getattr(module, proc)
		self.desc   = desc
		self.doc    = doc
		self._helps = Helps()

	def toHelps(self, helpsec):
		helpsec.prefix = ''
		helpsec.add((self.name, '', self.desc))

	@staticmethod
	def _parseSec(sec):
		ret      = {}
		sec      = Module.deindent(sec)
		lastdesc = []
		for s in sec:
			if not s:
				continue
			if s[0] in ('\t', ' '):
				lastdesc.append(s)
				continue
			parts = utils.split(s, ':', trim = False)
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

	@staticmethod
	def _parseDoc(doc):
		if not doc:
			return {}
		ret = {}
		lastsec = []
		for line in doc:
			if line.startswith('@'):
				secname = line.strip('@: ')
				lastsec = []
				ret[secname] = lastsec
			else:
				lastsec.append(line)
		for key, sec in ret.items():
			ret[key] = Process._parseSec(sec)
		return ret

	@staticmethod
	def defaultVal(val):
		return str(val) if val != '' else "''"

	def helps(self):
		"""Construct help page using doc"""
		if self._helps:
			return self._helps

		docs = Process._parseDoc(self.doc)
		self._helps.add('Name', Path(self.module.__file__).stem + '.' + self.name)
		self._helps.add('Description', self.desc)

		# input
		self._helps.add('Input options (Use \':list\' for multi-jobs)',
			sectype = 'option', prefix = '-')
		if isinstance(self.proc.config.input, dict):
			inputs = list(self.proc.config.input.keys())
		elif isinstance(self.proc.config.input, str):
			inputs = utils.split(self.proc.config.input, ',')
		else:
			inputs = self.proc.config.input
		for inopt in inputs:
			if ':' not in inopt:
				inopt += ':'
			inname, intype = inopt.split(':', 1)
			docinfo = docs.get('input', {}).get(inname, ('', '[ Not documented. ]'))
			self._helps.select('Input options').add((
				'-i.' + inname, '<%s>' % (intype or docinfo[0]) if intype or docinfo[0] \
					else '<var>',
				docinfo[1]))

		# output
		self._helps.add('Output options (\'exdir\' implied if path specified)',
			sectype = 'option', prefix = '-')
		if isinstance(self.proc.config.output, dict):
			outputs = list(self.proc.config.output.keys())
		elif isinstance(self.proc.config.input, str):
			outputs = utils.split(self.proc.config.output, ',')
		else:
			outputs = self.proc.config.output
		for outopt in outputs:
			parts = utils.split(outopt, ':')
			if len(parts) == 2:
				outname, outype, default = parts[0], 'var', parts[1]
			else:
				outname, outype, default = parts

			docinfo = docs.get('output', {}).get(
				outname, ('', ['Default: ' + Process.defaultVal(default)]))
			desc = docinfo[1]

			if not desc or ('default: ' not in desc[-1].lower() and len(desc[-1]) > 20):
				desc.append('Default: ' + Process.defaultVal(default))
			elif 'default: ' not in desc[-1].lower():
				desc[-1] += ' Default: ' + Process.defaultVal(default)
			self._helps.select('Output options').add((
				'-o.' + outname, '<%s>' % (outype or docinfo[0]) if outype or docinfo[0] \
					else '<var>',
				desc))

		# args
		self._helps.add('Process arguments', sectype = 'option', prefix = '-')
		for key in self.proc.args:
			docinfo = docs.get('args', {}).get(key, ('', '[ Not documented. ]'))
			self._helps.select('Process arguments').add((
				'-args.' + key, '<%s>' % docinfo[0] if docinfo[0] else '<auto>',
				docinfo[1]))

		self._helps.add('Other options', sectype = 'option', prefix = '-')
		# process properties
		self._helps.select('Other options').add(
			('-<prop>', '', 'Process properties, such as -forks, -exdir, -cache ...'))

		# pipeline configurations
		self._helps.select('Other options').add(('-config.<subconf>[.subconf]', '', 'Pipeline configrations, such as -config._log.file'))

		# help
		self._helps.select('Other options').addParam(
			params[params._hopts[0]], params._hopts, ishelp = True)

		return self._helps

	def printHelps(self, quit = True):
		print('\n'.join(HelpAssembler().assemble(self.helps())))
		if quit:
			sys.exit(1)

	def run(self, opts):
		if any(opts.get(h) for h in params._hopts) or \
			all(key in params._hopts for key in opts):
			self.printHelps()



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

	def run(self, args, prog = sys.argv[2:]):
		self.module.main(args, prog)
