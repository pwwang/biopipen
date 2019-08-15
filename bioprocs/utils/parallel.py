import threading
from time import sleep
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from bioprocs.utils.shell2 import runcmd
from traceback import format_exc
from queue import Queue

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
	nthread = min(total, nthread)
	(m, d) = divmod(total, nthread)
	ret = [m] * nthread
	for i in range(d):
		ret[i] += 1
	return ret

def distributeList(joblist, nthread):
	lists = distribute(len(joblist), nthread)
	start = 0
	for l in lists:
		yield joblist[start:l+start]
		start += l

class Parallel(object):

	def __init__(self, nthread = 1, backend = 'thread', raiseExc = False):
		PoolExecutor   = ProcessPoolExecutor if backend.lower() in 'multiprocessing' else ThreadPoolExecutor
		self.executor  = PoolExecutor(max_workers = nthread)
		self.raiseExc  = raiseExc

	def __del__(self):
		if self.executor:
			self.executor.shutdown()

	@staticmethod
	def _run(func, arg):
		if callable(func):
			ret = func(*arg)
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

class LargeFileHandler(object):
	HAND = -1

	def __init__(self, reader, compute = None, summarize = None, compute_multi = None, summarize_multi = None, size = 1000, readfunc = lambda r: next(r), nthread = 1, raiseExc = True):
		LargeFileHandler.HAND = -1
		self.reader           = reader
		self.readfunc         = readfunc
		self.raiseExc         = raiseExc
		self.size             = size
		self.queue            = Queue()
		self.nthread          = nthread
		self.compute          = compute
		self.summarize        = summarize
		self.compute_multi    = compute_multi
		self.summarize_multi  = summarize_multi
		self.index            = 0
		self.lock             = threading.Lock()
		self.stop             = False
		if self.compute and self.compute_multi:
			raise ValueError('Only one of compute and compute_multi is needed.')
		if self.summarize and self.summarize_multi:
			raise ValueError('Only one of summarize and summarize_multi is needed.')
		if not callable(self.compute) and not callable(self.compute_multi):
			raise ValueError('One of compute and compute_multi is needed.')
		if not callable(self.summarize) and not callable(self.summarize_multi):
			raise ValueError('One of summarize and summarize_multi is needed.')

		for _ in range(nthread):
			self.producer()

	def producer(self):
		if self.stop: return
		with self.lock:
			lines = []
			for _ in range(self.size):
				try:
					line = self.readfunc(self.reader)
					if line is False: break
					if line is None:
						self.stop = True
						break
					lines.append(line)
				except StopIteration:
					break
			if lines:
				self.queue.put((self.index, lines))
				self.index += 1

	def run(self):
		for _ in range(self.nthread):
			t = threading.Thread(target = self.worker)
			t.daemon = True
			t.start()
		self.queue.join()

	def worker(self):
		while not self.queue.empty():
			index, data = self.queue.get()
			try:
				if self.compute:
					computed = [self.compute(d) for d in data]
				else:
					computed = self.compute_multi(data)
				while LargeFileHandler.HAND + 1 != index:
					sleep(.01)
				if self.summarize:
					for c in computed:
						self.summarize(c)
				else:
					self.summarize_multi(computed)
				self.producer()
			except:
				if self.raiseExc:
					raise
			finally:
				LargeFileHandler.HAND = index
			self.queue.task_done()
