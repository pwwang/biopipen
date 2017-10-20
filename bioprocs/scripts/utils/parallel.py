
if 'parallel' not in vars() or not callable(parallel):
	from threading import Thread
	from Queue import Queue

	def _parallelWorker(sq):
		while True:
			if sq.empty(): 
				sq.task_done()
				break
			try:
				cmd2 = sq.get()
			except Exception:
				sq.task_done()
				break
			if not cmd2: 
				sq.task_done()
				break

			runcmd(cmd2)
			sq.task_done()

	def parallel(cmds, nthread):
		if nthread == 1 or len(cmds) == 1:
			for cmd in cmds: 
				runcmd(cmd)
		else:
			q = Queue()
			for cmd in cmds:
				q.put(cmd)
			for i in range(nthread):
				q.put(None)

			for i in range(nthread):
				t = Thread(target = _parallelWorker, args=(q, ))
				t.setDaemon(True)
				t.start()
			q.join()