
import re
import bioprocs
from pathlib import Path
from colorama import Back

def substr_replace(s, starts, lengths, replace):
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
	return substr_replace(origin, occurs, lengths, [
		'{}{}{}'.format(Back.LIGHTRED_EX, origin[occur:occur+length], Back.RESET)
		for occur, length in zip(occurs, lengths)])

class Module:

	@staticmethod
	def modules():
		return [module.stem
			for module in Path(bioprocs.__file__).parent.glob('*.py')
			if not module.stem.startswith('_')]

	def __init__(self, name):
		self.name   = name
		self.module = __import__('bioprocs.%s' % name, fromlist = ['bioprocs'])
		self.desc   = self.module.__doc__ and self.module.__doc__.strip() or '[ Not documented. ]'

	def procs(self):
		return {proc: Proc(self.name + '.' + proc)
				for proc in dir(self.module)
				if len(proc) > 1 and proc[0] == 'p' and (
					proc[1] == '_' or proc[1].isdigit() or proc[1].isupper()
				)}

class Process:

	def __init__(self, modproc):
		pass
