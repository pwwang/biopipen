import helpers, testly
from time import time, sleep
from bioprocs.utils.parallel import Parallel

class TestParallel(helpers.TestCase):

	def dataProvider_testParallel(self):
		yield lambda x: sleep(x), [(.1)]*5, 5, 'thread', .5
		yield lambda x: sleep(x), [(.1,)]*5, 5, 'process', .5
		yield 'sleep {}', [(.1,)]*5, 5, 'thread', .5
		yield 'sleep {}', [(.1,)]*5, 5, 'process', .5

	def testParallel(self, func, args, nthread, backend, sec):
		t0 = time()
		Parallel(nthread, backend).run(func,args)
		self.assertLess(time() - t0, sec)

	def dataProvider_testParallelEx(self):
		yield lambda x: 1/0, [(1,)]*2, 1, 'process', ZeroDivisionError
		yield lambda x: 1/0, [(1,)]*2, 2, 'thread', ZeroDivisionError

	def testParallelEx(self, func, args, nthread, backend, exception):
		self.assertRaises(exception, Parallel(nthread, backend, True).run, func, args)

	def dataProvider_testParallelReturn(self):
		yield lambda : 1, [()]*5, 5, 'thread', [1]*5
		yield lambda x, y: x + y, [(1, 2)]*5, 5, 'thread', [3]*5
		yield lambda : 1, [()]*5, 5, 'process', [1]*5
		yield lambda x, y: x + y, [(1, 2)]*5, 5, 'process', [3]*5

	def testParallelReturn(self, func, args, nthread, backend, out):
		ret = Parallel(nthread, backend, False).run(func, args)
		self.assertCountEqual(ret, out)

if __name__ == '__main__':
	testly.main()
