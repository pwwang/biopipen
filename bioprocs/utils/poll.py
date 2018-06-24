import filelock
from os import path
from time import sleep
from bioprocs.utils import runcmd, log2pyppl

class HardFileLock(filelock.BaseFileLock):
	"""
	The reverse of SortFileLock.
	With HardFileLock, if file exists, wait
	Here we do it reversely, it file does not exists, wait, until it's created.
	Note that the lockfile in <outdir> will always be removed as job re-created.
	"""
	def _acquire(self):
		try:
			self._lock_file_fd = open(self._lock_file, 'r')
		except (IOError, OSError):
			# wait
			self._lock_file_fd = None
		
	def _release(self):
		self._lock_file_fd = None


class Poll(object):
	"""
	Initiate a polling isinstance
	This is specifically applied in PyPPL framework
	@params:
		workdir    : The workdir of the process
		joblen     : How many jobs are there in this process
		jobindex   : What's the index for this job
		todo       : What to do while polling
		args, kwargs: The arguments for todo
	"""
	def __init__(self, workdir, joblen, jobindex):
		self.workdir  = workdir
		self.joblen   = joblen
		self.jobindex = jobindex
		
	"""
	A helper function to poll `func`, if it is true, then wait
	@params:
		func:  The function to tell if we need to wait
		args, kwargs: The arguments for func
	"""
	@staticmethod
	def wait(func, *args, **kwargs):
		interval = .1
		if 'interval' in kwargs:
			interval = kwargs['interval']
			del kwargs['interval']
		while func(*args, **kwargs):
			sleep(interval)


	"""
	Poll and do stuff at the first job, other jobs just wait
	"""
	def first(self, todo, *args, **kwargs):
		lockfilename = 'poll.first.lock'
		if 'lockfile' in kwargs:
			lockfilename = kwargs['lockfile']
			del kwargs['lockfile']
		# make sure it's cleaned when job reset
		lockfilename = path.join(self.workdir, '1', 'output', lockfilename)
		lockfile     = HardFileLock(lockfilename)
		if self.jobindex == 0:
			#log2pyppl('JOB #0: DOING stuff ... ')
			#with lockfile:
			#lockfile.acquire()
			try:
				if callable(todo):
					todo(*args, **kwargs)
				else:
					runcmd(todo.format(*args, **kwargs))
			except:
				raise
			finally:
				open(lockfilename, 'w').close()
			#log2pyppl('JOB #0: DOING stuff ... done')
		else:
			#log2pyppl('JOB #x: waiting ... for flag file')
			#Poll.wait(lambda x: not path.exists(x), lockfilename)
			#log2pyppl('JOB #x: waiting ... for job #0')
			with lockfile: pass
			#log2pyppl('JOB #x: waiting ... done')

	"""
	Poll jobs other than the first job to do stuff, first job just wait.
	Make sure you have more than 2 jobs running (proc.forks > 1), otherwise, the first job
	will wait forever
	"""
	def non1st(self, todo, *args, **kwargs):
		lockfilename = 'poll.non1st.lock'
		if 'lockfile' in kwargs:
			lockfilename = kwargs['lockfile']
			del kwargs['lockfile']
		lockfilenames = [
			path.join(self.workdir, str(jobindex + 1), 'output', lockfilename) \
			for jobindex in range(self.joblen)
		]
		lockfiles = [HardFileLock(f) for f in lockfilenames]
		if self.jobindex == 0:
			for i, lockfilename in enumerate(lockfilenames):
				if i == 0: continue
				#Poll.wait(lambda x: not path.exists(x), lockfilename)
				with lockfiles[i]: pass
		else:
			#with lockfiles[self.jobindex]:
			#lockfiles[self.jobindex].acquire()
			try:
				if callable(todo):
					todo(*args, **kwargs)
				else:
					runcmd(todo.format(*args, **kwargs))
			except:
				raise
			finally:
				open(lockfilename[self.jobindex], 'w').close()
				#lockfiles[self.jobindex].release()

	"""
	You have to have all jobs running
	Wait for all jobs done the same thing (todo)
	"""
	def all(self, todo, *args, **kwargs):
		lockfilename = 'poll.all.lock'
		if 'lockfile' in kwargs:
			lockfilename = kwargs['lockfile']
			del kwargs['lockfile']
		lockfilenames = [
			path.join(self.workdir, str(jobindex + 1), 'output', lockfilename) \
			for jobindex in range(self.joblen)
		]
		lockfiles = [HardFileLock(f) for f in lockfilenames]
		
		#log2pyppl('Doing stuff at job #%s ... ' % self.jobindex)
		# with lockfiles[self.jobindex]:
		#lockfiles[self.jobindex].acquire()
		try:
			if callable(todo):
				todo(*args, **kwargs)
			else:
				runcmd(todo.format(*args, **kwargs))
		except:
			raise
		finally:
			open(lockfilenames[self.jobindex], 'w').close()
			#lockfiles[self.jobindex].release()
		#log2pyppl('Doing stuff at job #%s ... done ' % self.jobindex)

		for i, lockfile in enumerate(lockfiles):
			if i == self.jobindex: continue
			#log2pyppl('Waiting for lock file: %s ... ' % self.jobindex)
			#Poll.wait(lambda x: not path.exists(x), lockfilenames[i])
			#log2pyppl('Waiting for job: %s ... ' % self.jobindex)
			with lockfile: pass
			#log2pyppl('Waiting for job: %s ... done ' % self.jobindex)
