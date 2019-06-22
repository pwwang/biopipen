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
commands.list.all.desc = 'List all processes if no module specified.'
commands.list._        = []
commands.list._.desc   = 'Module name'
