from pyparam import commands

commands._prefix       = '-'

# command: params
commands.params            = 'Show or query a parameter'
commands.params._.type     = list
commands.params._.required = True
commands.params._.desc     = 'The name of the parameter.'
commands.params.get        = False
commands.params.get.desc   = 'Get the value.'

# command: list
commands.list          = 'List modules or process'
commands.list._hbald   = False
commands.list.all      = False
commands.list.a        = commands.list.all
commands.list.all.desc = 'List all processes if no module specified.'
commands.list._        = []
commands.list._.desc   = 'Module name.'

# command: proc
commands.proc = 'Search modules or processes.'
commands.proc._        = []
commands.proc._.desc   = 'Part of the module or process name.'

# command: completion
commands.completion = 'Generate completions for bioprocs.'
commands.completion.shell = 'auto'
commands.completion.shell.desc = 'The shell'
commands.completion.auto = False
commands.completion.auto.desc = 'Whether automatically deploy the completion script.'
commands.completion.s = commands.completion.shell
commands.completion.a = commands.completion.auto
commands.completions  = commands.completion

# command: profile
commands.profile        = 'List available running profiles.'
commands.profile._hbald = False
commands.profile._.desc = 'Additional configurations.'

def helpx(helps):
	from bioprocs.console.utils import Pipeline
	#helps.add('Parameter/Process query commands', helps.select('Available commands'))
	helps.add('Assembled pipelines', sectype = 'option', prefix = '')
	helps.add('Processes', sectype = 'option', prefix = '')
	helps.add('Other commands', sectype = 'option', prefix = '')
	help_avail    = helps.select('Available')
	help_pipeline = helps.select('Assembled')
	help_process  = helps.select('Processes')
	help_other    = helps.select('Other')

	for ppl in Pipeline.pipelines():
		help_pipeline.add((ppl, '', Pipeline(ppl).desc))

	help_process.add(help_avail.select('list'))
	help_process.add(help_avail.select('proc'))
	help_process.add(('<module>', '', [
		'List the process of module.', 'Same as "bioprocs list <module>"']))
	help_process.add(('<module.proc>', '', 'Run a process in command line.'))
	help_other.add(help_avail.select('completion'))
	help_other.add(help_avail.select('params'))
	help_other.add(help_avail.select('profile'))
	help_other.add(help_avail.select('help'))

	helps.remove('Available')


# modify help page
commands._helpx = helpx
