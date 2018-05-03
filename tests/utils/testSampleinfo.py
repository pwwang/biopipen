import helpers, unittest
from os import path
from bioprocs.utils.sampleinfo import SampleInfo, SampleInfoException
from bioprocs.utils import runcmd

class TestSampleInfo(helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestSampleInfo')

	def dataProvider_testSampleInfo(self):
		sfile = path.join(self.indir, 'sampleinfo1.txt')
		yield sfile, 8, 4, ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8'], ['ROWNAMES', 'Patient', 'Group', 'Batch']

	def testSampleInfo(self, sfile, nrow, ncol, rownames, colnames):
		si = SampleInfo(sfile)
		self.assertEqual(si.nrow, nrow)
		self.assertEqual(si.ncol, ncol)
		self.assertEqual(si.rownames, rownames)
		self.assertEqual(si.colnames, colnames)
		self.assertEqual(si.samcol, colnames[0])

	def dataProvider_testSelect(self):
		sifile = path.join(self.indir, 'sampleinfo1.txt')
		si = SampleInfo(sifile)
		yield si, 'S2', None, None, None, None, si.data[1:2]
		yield si, None, 'P1', None, None, None, si.data[0:2]
		yield si, None, 'P1', 'G1', None, None, si.data[0:1]
		yield si, None, 'P2', None, 'B1', None, si.data[2:4]
		yield si, None, 'P2', None, 'B1', 'Sample', ['S3', 'S4']

	def testSelect(self, si, sample, patient, group, batch, get, items):
		self.maxDiff = None
		self.assertListEqual(si.select(sample, patient, group, batch, get), items)

	def dataProvider_testSampleInfoR(self):
		rfile = path.join(self.indir, 'sampleinfo.r')
		sfile = path.join(self.indir, 'sampleinfo1.txt')
		yield rfile, sfile, 8, 3

	def testSampleInfoR(self, rfile, sfile, nrow, ncol):
		with helpers.captured_output():
			r = runcmd(['Rscript', rfile, sfile, nrow, ncol])
		self.assertTrue(r)

if __name__ == '__main__':
	unittest.main(verbosity = 2)
