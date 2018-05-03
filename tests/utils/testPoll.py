import helpers, testly
import random
from time import time, sleep
from bioprocs.utils import runcmd, logger
from bioprocs.utils.poll import Poll
from bioprocs.utils.parallel import Parallel
from helpers import testdirs
from os import path, makedirs

class TestPoll(helpers.TestCase):
	testdir, indir, outdir = testdirs('TestPoll')

	def dataProvider_testWait(self):
		yield (lambda : sleep(.1) or False), [], {}, .1
		d = [6, 7, 3]
		yield lambda x: d.pop(0) > x, [5], dict(interval = .2), .4

	def testWait(self, func, args, kwargs, out):
		t0 = time()
		Poll.wait(func, *args, **kwargs)
		t = time() - t0
		self.assertGreater(t, out)

	def dataProvider_testFirst(self):
		workdir = path.join(self.testdir, 'testfirst')
		if not path.exists(workdir): makedirs(workdir)
		firstoutdir = path.join(workdir, '1', 'output')
		if not path.exists(firstoutdir): makedirs(firstoutdir)
		yield workdir, 5, (lambda x: sleep(x)), [.1], {}, .5

	def testFirst(self, workdir, joblen, todo, args, kwargs, out):
		p = Parallel(nthread = joblen, raiseExc = True, backend = 'thread')
		def pfunc(jobindex):
			t0 = time()
			Poll(workdir, joblen, jobindex).first(todo, *args, **kwargs)
			return time() - t0
		ret = p.run(pfunc, [(i, ) for i in range(joblen)])
		self.assertGreater(sum(ret), out) # other jobs waited
		self.assertTrue(path.exists(path.join(workdir, '1', 'output', 'poll.first.lock')))
		for i in range(joblen):
			if i == 0: continue
			self.assertFalse(path.exists(path.join(workdir, str(i+1), 'output', 'poll.first.lock')))

	def dataProvider_testNon1st(self):
		workdir = path.join(self.testdir, 'testnon1st')
		if not path.exists(workdir): makedirs(workdir)
		joblen = 5
		for i in range(joblen):
			outdir = path.join(workdir, str(i+1), 'output')
			if not path.exists(outdir): makedirs(outdir)
		yield workdir, 5, (lambda x: sleep(x)), [.1], {}, .5

	def testNon1st(self, workdir, joblen, todo, args, kwargs, out):
		p = Parallel(nthread = joblen, raiseExc = True, backend = 'thread')
		def pfunc(jobindex):
			t0 = time()
			Poll(workdir, joblen, jobindex).non1st(todo, *args, **kwargs)
			return time() - t0
		ret = p.run(pfunc, [(i, ) for i in range(joblen)])
		self.assertGreater(sum(ret), out)
		self.assertFalse(path.exists(path.join(workdir, '1', 'output', 'poll.non1st.lock')))
		for i in range(joblen):
			if i == 0: continue
			self.assertTrue(path.exists(path.join(workdir, str(i+1), 'output', 'poll.non1st.lock')))

	def dataProvider_testAll(self):
		workdir = path.join(self.testdir, 'testall')
		if not path.exists(workdir): makedirs(workdir)
		joblen = 5
		for i in range(joblen):
			outdir = path.join(workdir, str(i+1), 'output')
			if not path.exists(outdir): makedirs(outdir)
		yield workdir, 5, (lambda x: sleep(x)), [.1], {}, .5

	def testAll(self, workdir, joblen, todo, args, kwargs, out):
		p = Parallel(nthread = joblen, raiseExc = True, backend = 'thread')
		def pfunc(jobindex):
			t0 = time()
			Poll(workdir, joblen, jobindex).all(todo, *args, **kwargs)
			return time() - t0
		ret = p.run(pfunc, [(i, ) for i in range(joblen)])
		self.assertGreater(sum(ret), out)
		for i in range(joblen):
			self.assertTrue(path.exists(path.join(workdir, str(i+1), 'output', 'poll.all.lock')))

class TestPollR(helpers.TestCase):
	pollr = None
	testdir, indir, outdir = testdirs('TestPollR')

	@staticmethod
	def setUpClass():
		import inspect
		pollrfile  = path.join(path.dirname(inspect.getmodule(Poll).__file__), 'poll.r')
		with open(pollrfile, 'r') as f:
			lib = f.read()
		TestPollR.writer = lambda self, sfile, code: [f.write(lib + '\n' + code) or f.close() for f in [open(sfile, 'w')]]

	def dataProvider_testWait(self):
		yield 'function(x) {Sys.sleep(.1); F}', 'list()', .1
		yield 'function(x) {v = d[1]; d <<- d[-1]; v > x}', 'list(5, interval = .2)', .4

	def testWait(self, func, kwargs, out):
		rscript = path.join(self.testdir, self._testMethodName + '.r')
		self.writer(rscript, """
			d = c(6, 7, 3)
			do.call(wait, c(list({func}), {kwargs}))
		""".format(func = func, kwargs = kwargs))
		t0 = time()
		runcmd(['Rscript', rscript])
		t = time() - t0
		self.assertGreater(t, out)

	def dataProvider_testFirst(self):
		workdir = path.join(self.testdir, 'testfirst')
		if not path.exists(workdir): makedirs(workdir)
		firstoutdir = path.join(workdir, '1', 'output')
		if not path.exists(firstoutdir): makedirs(firstoutdir)
		yield workdir, 5, 'function(x) Sys.sleep(x)', 'list(.1)', .5

	def testFirst(self, workdir, joblen, todo, kwargs, out):
		p = Parallel(nthread = joblen, raiseExc = True)

		def pfunc(jobindex):
			rscript = path.join(self.testdir, self._testMethodName + '-{}.r'.format(jobindex))
			self.writer(rscript, """
				do.call(Poll('{workdir}', {joblen}, {jobindex})$first, c(list({todo}), {kwargs}))
			""".format(workdir = workdir, joblen = joblen, jobindex = jobindex, todo = todo, kwargs = kwargs))
			t0 = time()
			runcmd(['Rscript', rscript])
			return time() - t0
		ret = p.run(pfunc, [(i, ) for i in range(joblen)])
		self.assertGreater(sum(ret), out) # other jobs waited
		self.assertTrue(path.exists(path.join(workdir, '1', 'output', 'poll.first.lock')))
		for i in range(joblen):
			if i == 0: continue
			self.assertFalse(path.exists(path.join(workdir, str(i+1), 'output', 'poll.first.lock')))


	def dataProvider_testNon1st(self):
		workdir = path.join(self.testdir, 'testnon1st')
		if not path.exists(workdir): makedirs(workdir)
		joblen = 5
		for i in range(joblen):
			outdir = path.join(workdir, str(i+1), 'output')
			if not path.exists(outdir): makedirs(outdir)
		yield workdir, 5, 'function(x) Sys.sleep(x)', 'list(.1)', .5

	def testNon1st(self, workdir, joblen, todo, kwargs, out):
		p = Parallel(nthread = joblen, raiseExc = True)
		def pfunc(jobindex):
			rscript = path.join(self.testdir, self._testMethodName + '-{}.r'.format(jobindex))
			self.writer(rscript, """
				do.call(Poll('{workdir}', {joblen}, {jobindex})$non1st, c(list({todo}), {kwargs}))
			""".format(workdir = workdir, joblen = joblen, jobindex = jobindex, todo = todo, kwargs = kwargs))
			t0 = time()
			runcmd(['Rscript', rscript])
			return time() - t0
		ret = p.run(pfunc, [(i, ) for i in range(joblen)])
		self.assertGreater(sum(ret), out)
		self.assertFalse(path.exists(path.join(workdir, '1', 'output', 'poll.non1st.lock')))
		for i in range(joblen):
			if i == 0: continue
			self.assertTrue(path.exists(path.join(workdir, str(i+1), 'output', 'poll.non1st.lock')))

	def dataProvider_testAll(self):
		workdir = path.join(self.testdir, 'testall')
		if not path.exists(workdir): makedirs(workdir)
		joblen = 5
		for i in range(joblen):
			outdir = path.join(workdir, str(i+1), 'output')
			if not path.exists(outdir): makedirs(outdir)
		yield workdir, 5, 'function(x) Sys.sleep(x)', 'list(.1)', .5

	def testAll(self, workdir, joblen, todo, kwargs, out):
		p = Parallel(nthread = joblen, raiseExc = True)
		def pfunc(jobindex):
			rscript = path.join(self.testdir, self._testMethodName + '-{}.r'.format(jobindex))
			self.writer(rscript, """
				do.call(Poll('{workdir}', {joblen}, {jobindex})$all, c(list({todo}), {kwargs}))
			""".format(workdir = workdir, joblen = joblen, jobindex = jobindex, todo = todo, kwargs = kwargs))
			t0 = time()
			runcmd(['Rscript', rscript])
			return time() - t0
		ret = p.run(pfunc, [(i, ) for i in range(joblen)])
		self.assertGreater(sum(ret), out)
		for i in range(joblen):
			self.assertTrue(path.exists(path.join(workdir, str(i+1), 'output', 'poll.all.lock')))

if __name__ == '__main__':
	testly.main(verbosity = 2)
