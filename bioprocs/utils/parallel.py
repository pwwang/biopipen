from bioprocs.utils.helpers import runcmd

def parallel(func, args, nthread, method = 'thread'):
	"""
	Call functions in a parallel way.
	If nthread == 1, will be running in single-threading manner.
	@params:
		`func`: The function
		`args`: The arguments, in list. Each element should be the arguments for the function in one thread.
		`nthread`: Number of threads
		`method`: use multithreading (thread) or multiprocessing (process)
	"""
	if method == 'thread':
		Queue = moves.queue.Queue
		Unit  = Thread
	else:
		Queue = JoinableQueue
		Unit  = Process

	def _parallelWorker(q):
		while True:
			if q.empty(): # pragma: no cover
				q.task_done()
				break
			arg = q.get()
			if arg is None:
				q.task_done()
				break
			
			try:
				func(*arg)
			except Exception: # pragma: no cover
				raise
			finally:
				q.task_done()

	nworkers = min(len(args), nthread)
	q = Queue()

	for arg in args: q.put(arg)
	for _ in range(nworkers): q.put(None)

	for _ in range(nworkers):
		t = Unit(target = _parallelWorker, args=(q, ))
		t.daemon = True
		t.start()
	
	q.join()
	
def parallelCmd(cmds, nthread, method = 'thread'):
	parallel(runcmd, [(cmd,) for cmd in cmds], nthread, method)