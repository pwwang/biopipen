from concurrent.futures import ThreadPoolExecutor
from loky import ProcessPoolExecutor
from bioprocs.utils import runcmd
from traceback import format_exc

class Parallel(object):

	def __init__(self, nthread = 1, backend = 'process', raiseExc = False):
		PoolExecutor   = ProcessPoolExecutor if backend.lower() in 'multiprocessing' else ThreadPoolExecutor
		self.executor  = PoolExecutor(max_workers = nthread)
		self.raiseExc  = raiseExc

	def run(self, func, args):
		def _func(arg):
			if callable(func):
				return func(*arg)
			else:
				runcmd(func.format(*arg))

		submits   = []
		results   = []
		exception = None
		excno     = 0
		for arg in args:
			submits.append(self.executor.submit(_func, arg))

		for submit in submits:
			try:
				results.append(submit.result())
			except Exception as ex:
				#results.append(None)
				exception = type(ex), format_exc()
				excno    += 1

		if excno > 0 and self.raiseExc:
			raise exception[0](exception[1])

		return results
