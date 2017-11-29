import unittest, json

from pyppl import PyPPL, Channel
from helpers import getfile, procOK, config
from bioaggrs.common import aSimRead2


class testCommon (unittest.TestCase):

	def testSimRead(self):
		
		aSimRead2.input                = [[getfile('asimread1.txt')], [getfile('asimread2.txt')]]
		aSimRead2.pSort1.args.params.k = '1'
		aSimRead2.pSort1.args.skip     = 1
		aSimRead2.pSort2.args.params.k = '2'
		aSimRead2.args.skip            = [1]
		aSimRead2.args.match           = 'lambda line1, line2: -1 if line1[0] == line2[1] else 0 if line1[0] < line2[1] else 1'
		aSimRead2.args.do              = 'lambda line1, line2: fout.write(line1[0] + "\\t" + str(int(line1[1]) + int(line2[2])) + "\\n")'

		PyPPL().start(aSimRead2).run()
		procOK(aSimRead2.pSimRead, 'asimread.txt', self)


if __name__ == '__main__':
	unittest.main()