from pyparam import commands

commands._prefix       = '-'

# command: params
commands.params            = 'Show or query a parameter'
print(commands.params._prefix)
commands.params._hbald     = False
commands.params.get        = False
commands.params.get.desc   = 'Get the value.'
commands.params.query.desc = 'The keyword used to query parameters.'


