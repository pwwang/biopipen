import addPath, unittest

from os import path
from pyppl import PyPPL, params
from bioaggrs.common import aSimRead2

class testAggrCommon (unittest.TestCase):

	def testSimRead(self):
		file1 = path.join(params.tmpdir.value, 'asimread1.txt')
		file2 = path.join(params.tmpdir.value, 'asimread2.txt')

		with open(file1, 'w') as f1, open(file2, 'w') as f2:
			f1.write('# something to be skipped\n')
			f1.write('m1\t2\n')
			f1.write('m2\t3\n')
			f1.write('m6\t4\n')
			f1.write('m3\t5\n')
			f1.write('m8\t6\n')
			f1.write('m4\t7\n')
			f1.write('m5\t8\n')
			f2.write('#\tm1\t2\n')
			f2.write('#\tm2\t8\n')
			f2.write('#\tm6\t1\n')
			f2.write('#\tm3\t3\n')
			f2.write('#\tm8\t6\n')
			f2.write('#\tm4\t9\n')
			f2.write('#\tm5\t3\n')
		aSimRead2.input                = [[file1], [file2]]
		aSimRead2.pSort1.args.params.k = '1'
		aSimRead2.pSort1.args.skip     = 1
		aSimRead2.pSort2.args.params.k = '2'
		aSimRead2.args.skip            = [1]
		aSimRead2.args.match           = 'lambda line1, line2: -1 if line1[0] == line2[1] else 0 if line1[0] < line2[1] else 1'
		aSimRead2.args.do              = 'lambda line1, line2: fout.write(line1[0] + "\\t" + str(int(line1[1]) + int(line2[2])) + "\\n")'

		PyPPL().start(aSimRead2).run()
		with open(aSimRead2.pSimRead.channel.get()) as f:
			outs = f.read().splitlines()

		self.assertIn('m1\t4', outs)
		self.assertIn('m2\t11', outs)
		self.assertIn('m3\t8', outs)
		self.assertIn('m4\t16', outs)
		self.assertIn('m5\t11', outs)
		self.assertIn('m6\t5', outs)
		self.assertIn('m8\t12', outs)


if __name__ == '__main__':
	unittest.main()