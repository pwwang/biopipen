from concurrent.futures import ThreadPoolExecutor
from loky import ProcessPoolExecutor
from bioprocs.utils.shell import runcmd, ShellResult
from traceback import format_exc

def distribute(total, nthread):
	"""
	Try to distribute jobs into N threads as equal as possible.
	For example: distributing 10 jobs on 3 threads, we prefer (4, 3, 3) than (4, 4, 2)
	How to do it? 
	1. get the ceiling size for each thread, that should be (3, 3, 3), from `divmod(10, 3)[0]`
	2. get the modulo by `divmod(10, 3)[1]`, which is the first # nthreads to add one more job to
	@params:
		`total`: The total # jobs
		`nthread`: The # threads
	@returns:
		A list of # jobs distribute on each thread.
	"""
	(m, d) = divmod(total, nthread)
	ret = [m] * nthread
	for i in range(d):
		ret[i] += 1
	return ret

class Parallel(object):

	def __init__(self, nthread = 1, backend = 'thread', raiseExc = False):
		PoolExecutor   = ProcessPoolExecutor if backend.lower() in 'multiprocessing' else ThreadPoolExecutor
		self.executor  = PoolExecutor(max_workers = nthread)
		self.raiseExc  = raiseExc

	@staticmethod
	def _run(func, arg):
		print 3234, repr(func)
		if isinstance(func, ShellResult):
			return func(*arg).run()
		elif callable(func):
			print 999
			ret = func(*arg)
			print ret
			return ret
		else:
			return runcmd(func.format(*arg))

	def run(self, func, args):
		submits   = []
		results   = []
		exception = None
		excno     = 0
		for arg in args:
			submits.append(self.executor.submit(Parallel._run, func, arg))

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
