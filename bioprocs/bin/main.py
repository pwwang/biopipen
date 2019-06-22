from bioprocs.bin.utils import highlight, Module
from bioprocs.bin.arguments import commands
from bioprocs import params
from pyparam import HelpAssembler, Helps
from pyppl import Box

def showParams(opts):
	"""Show and query a parameter"""
	# get is preserved for Box
	if opts['get']:
		for pname in opts._:
			if pname not in params:
				raise ValueError('Parameter %r is not defined.' % pname)
			print(params[pname].value)
	else:
		pminfo  = Helps()
		pminfo.add('params', sectype = 'option', prefix = '')
		queries = [query.lower() for query in opts._]
		for key in sorted(params._dict()):
			if key in params._hopts:
				continue
			for query in queries:
				param = params[key]
				if query not in param.str().lower() and query not in key.lower() \
					and not any(query in desc for desc in param.desc):
					continue
				pminfo['params'].addParam(param)

		def hiq(line):
			for query in queries:
				line = highlight(line, query)
			return line
		print('\n'.join(hiq(line) for line in  HelpAssembler().assemble(pminfo)))

def listProcs(opts):
	pminfo  = Helps()
	pminfo.add('modules', sectype = 'option', prefix = '')
	if not opts._ and not opts.all:
		for module in Module.modules():
			module = Module(module)
			pminfo['modules'].add((module.name, '', module.desc))
		print('\n'.join(HelpAssembler().assemble(pminfo)))

def main():
	command, opts, _ = commands._parse(arbi = True, dict_wrapper = Box)
	if command == 'params':
		showParams(opts)
	elif command == 'list':
		listProcs(opts)

if __name__ == "__main__":
	main()
