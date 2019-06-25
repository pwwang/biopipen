import sys
from bioprocs.bin.utils import highlightMulti, Module, Pipeline
from bioprocs.bin.arguments import commands
from bioprocs import params
from pyparam import HelpAssembler, Helps
from pyppl import Box

HELP_ASSEMBLER = HelpAssembler()

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

		print('\n'.join(highlightMulti(line, queries) for line in  HELP_ASSEMBLER.assemble(pminfo)))

def listProcs(opts):
	if not opts._ and not opts.all: # list modules only
		pminfo  = Helps()
		pminfo.add('modules', sectype = 'option', prefix = '')
		for module in Module.modules():
			Module(module).toHelps(pminfo['modules'])
		print('\n'.join(HELP_ASSEMBLER.assemble(pminfo)))
	elif not opts._:  # list all modules and processes
		pminfo  = Helps()
		modules = Module.modules()
		nmods   = len(modules)
		for i, module in enumerate(modules):
			Module(module).toHelpsWithProcs(pminfo, nmods)
		print('\n'.join(HELP_ASSEMBLER.assemble(pminfo)))
	else: # list specific modules and their processes
		pminfo  = Helps()
		modules = Module.modules()
		nmods   = len(opts._)
		for module in opts._:
			if module not in modules:
				raise ValueError('No such module: %s' % module)
			Module(module).toHelpsWithProcs(pminfo, nmods)
		print('\n'.join(HELP_ASSEMBLER.assemble(pminfo)))

def proc(opts):
	modules = Module.modules()
	nmods   = len(modules)
	helps   = Helps()
	for i, module in enumerate(modules):
		module = Module(module)
		if any(q.lower() in module.name.lower() for q in opts._):
			module.toHelpsWithProcs(helps, nmods)
		else:
			module.toHelpsAsSec(helps)
			with module.loadProcs("%s/%s" % (i+1, nmods)) as procs:
				for name, proc in procs.items():
					if any(q.lower() in name.lower() for q in opts._):
						proc.toHelps(helps.select(module.name + ':'))
	print('\n'.join(highlightMulti(line, opts._) for line in  HELP_ASSEMBLER.assemble(helps)))


def main():
	pipelines = Pipeline.pipelines()
	command = sys.argv[1] if len(sys.argv) > 1 else None
	if command == 'params':
		command, opts, _ = commands._parse(dict_wrapper = Box)
		showParams(opts)
	elif command == 'list':
		command, opts, _ = commands._parse(dict_wrapper = Box)
		listProcs(opts)
	elif command == 'proc':
		command, opts, _ = commands._parse(dict_wrapper = Box)
		proc(opts)
	elif command in pipelines:
		Pipeline(command).run(sys.argv[2:])
	else:
		command, opts, _ = commands._parse(arbi = True, dict_wrapper = Box)
		if '.' in command:
			module, proc = command.split('.')
			module = Module(module)
			if proc not in module.procs():
				raise ValueError('No such process: %s' % command)
			module.procs()[proc].run(opts)
		else: # listing module processes
			listProcs(Box(_ = [command]))

if __name__ == "__main__":
	main()
