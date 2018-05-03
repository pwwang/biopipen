import helpers, unittest, testly
from pyppl import PyPPL
from os import path
from bioprocs.bed import pBedGetfasta, pBedFlank, pBedRandom, pBedSort, pBedClosest, pBedClosest2, pBedIntersect, pBedIntersect2, pBedMakewindows, pBedMerge, pBedMerge2
from helpers import getfile, procOK, config

class TestBed (helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestBed')

	def dataProvider_testBedGetfasta(self):
		infile = path.join(self.indir, 'test.getfasta.bed')
		outfile = path.join(self.outdir, 'test.getfasta.ret')
		args = {'ref': path.join(self.indir, 'test.getfasta.fa')}
		yield 't1', infile, args, outfile

	def testBedGetfasta(self, tag, infile, args, outfile):
		pBedGetfastaTest = pBedGetfasta.copy(tag = tag)
		pBedGetfastaTest.input = [infile]
		pBedGetfastaTest.args.update(args)
		PyPPL(config).start(pBedGetfastaTest).run()
		self.assertFileEqual(pBedGetfastaTest.channel.get(), outfile)

	def dataProvider_testBedFlank(self):
		infile = path.join(self.indir, 'bedflank.bed')
		outfile1 = path.join(self.outdir, 'bedflank1.bed')
		args1 = {'params': {'b': 5}}
		yield 't1', infile, args1, outfile1
		outfile2 = path.join(self.outdir, 'bedflank2.bed')
		args2 = {'extend': True, 'params': {'b': 5}}
		yield 't2', infile, args2, outfile2

	def testBedFlank (self, tag, infile, args, outfile):
		pBedFlankTest = pBedFlank.copy(tag = tag)
		pBedFlankTest.input = [infile]
		pBedFlankTest.args.update(args)
		PyPPL().start(pBedFlankTest).run()
		self.assertFileEqual(pBedFlankTest.channel.get(), outfile)

	def dataProvider_testBedRandom(self):
		yield 't1', 100, 10, 8525, path.join(self.outdir, 'bedrandom.bed')
		yield 't2', 1000, 10, 8525, path.join(self.outdir, 'bedrandom2.bed')


	def testBedRandom(self, tag, length, n, seed, outfile):
		pBedRandomTest = pBedRandom.copy(tag = tag)
		pBedRandomTest.input = length, n
		pBedRandomTest.args.seed = seed
		PyPPL(config).start(pBedRandomTest).run()
		self.assertFileEqual(pBedRandomTest.channel.get(), outfile)

	def dataProvider_testBedSort(self):
		infile = path.join(self.indir, 'bedsort.bed')
		yield 't1', infile, {'tool': 'sort', 'unique': False}, path.join(self.outdir, 'bedsort.bed')
		yield 't2', infile, {'tool': 'sort', 'unique': True}, path.join(self.outdir, 'bedsort-unique.bed')
		yield 't3', infile, {'tool': 'bedops', 'unique': False}, path.join(self.outdir, 'bedsort.bed')
		yield 't4', infile, {'tool': 'bedops', 'unique': True}, path.join(self.outdir, 'bedsort-unique.bed')
		yield 't5', infile, {'tool': 'bedtools', 'unique': False}, path.join(self.outdir, 'bedsort.bed')
		yield 't6', infile, {'tool': 'bedtools', 'unique': True}, path.join(self.outdir, 'bedsort-unique.bed')

	def testBedSort(self, tag, infile, args, outfile):
		pBedSortTest = pBedSort.copy(tag = tag)
		pBedSortTest.input = [infile]
		pBedSortTest.args.update(args)
		PyPPL(config).start(pBedSortTest).run()
		self.assertFileEqual(pBedSortTest.channel.get(), outfile)

	def dataProvider_testBedClosest(self):
		afile = path.join(self.indir, 'bedclosest.a.bed')
		bfile = path.join(self.indir, 'bedclosest.b.bed')
		outfile = path.join(self.outdir, 'bedclosest1.bt')
		outfile2 = path.join(self.outdir, 'bedclosest2.bt')
		yield 't1', afile, bfile, {}, outfile
		yield 't2', afile, bfile, {"io": True}, outfile2

	def testBedClosest(self, tag, afile, bfile, params, outfile):
		pBedClosestTest = pBedClosest.copy(tag = tag)
		pBedClosestTest.input = afile, bfile
		pBedClosestTest.args.params.update(params)
		PyPPL(config).start(pBedClosestTest).run()
		self.assertFileEqual(pBedClosestTest.channel.get(), outfile)

	def dataProvider_testBedClosest2(self):
		afile = path.join(self.indir, 'bedclosest.a.bed')
		bfiles = [
			path.join(self.indir, 'bedclosest.b.bed'),
			path.join(self.indir, 'bedclosest.b2.bed'),
		]
		outfile = path.join(self.outdir, 'bedclosest2_1.bt')
		outfile2 = path.join(self.outdir, 'bedclosest2_2.bt')
		yield 't1', afile, bfiles, {}, outfile
		yield 't2', afile, bfiles, {"mdb": 'all', "io": True}, outfile2

	def testBedClosest2(self, tag, afile, bfiles, params, outfile):
		pBedClosest2Test = pBedClosest2.copy(tag = tag)
		pBedClosest2Test.input = afile, bfiles
		pBedClosest2Test.args.params.update(params)
		PyPPL(config).start(pBedClosest2Test).run()
		self.assertFileEqual(pBedClosest2Test.channel.get(), outfile)

	def dataProvider_testBedIntersect(self):
		afile = path.join(self.indir, 'bedintersect.a.bed')
		bfile = path.join(self.indir, 'bedintersect.b.bed')
		outfile = path.join(self.outdir, 'bedintersect.out.bt')
		yield 't1', afile, bfile, {'wao': False}, outfile

	def testBedIntersect(self, tag, afile, bfile, params, outfile):
		pBedIntersectTest = pBedIntersect.copy(tag = tag)
		pBedIntersectTest.input = afile, bfile
		pBedIntersectTest.args.params.update(params)
		PyPPL(config).start(pBedIntersectTest).run()
		self.assertFileEqual(pBedIntersectTest.channel.get(), outfile)

	def dataProvider_testBedIntersect2(self):
		afile = path.join(self.indir, 'bedintersect2.a.bed')
		bfiles = [
			path.join(self.indir, 'bedintersect2.b1.bed'),
			path.join(self.indir, 'bedintersect2.b2.bed'),
			path.join(self.indir, 'bedintersect2.b3.bed'),
		]
		outfile = path.join(self.outdir, 'bedintersect2.out.bt')
		yield 't1', afile, bfiles, {'wao': False}, outfile

	def testBedIntersect2(self, tag, afile, bfiles, params, outfile):
		pBedIntersect2Test = pBedIntersect2.copy(tag = tag)
		pBedIntersect2Test.input = afile, bfiles
		pBedIntersect2Test.args.params.update(params)
		PyPPL(config).start(pBedIntersect2Test).run()
		self.assertFileEqual(pBedIntersect2Test.channel.get(), outfile)

	def dataProvider_testBedMakewindows(self):
		infile = path.join(self.indir, 'bedclosest.a.bed')
		args   = {'intype': 'bed', 'params': {'n': 10}}
		outfile = path.join(self.outdir, 'bedmakewindows.out.bed')
		yield 't1', infile, args, outfile

	def testBedMakewindows(self, tag, infile, args, outfile):
		pBedMakewindowsTest = pBedMakewindows.copy(tag = tag)
		pBedMakewindowsTest.input = [infile]
		pBedMakewindowsTest.args.update(args)
		PyPPL(config).start(pBedMakewindowsTest).run()
		self.assertFileEqual(pBedMakewindowsTest.channel.get(), outfile)

	def dataProvider_testBedMerge(self):
		infile = path.join(self.indir, 'bedmerge.bed')
		outfile = path.join(self.outdir, 'bedmerge.out.bed')
		yield 't1', infile, outfile

	def testBedMerge(self, tag, infile, outfile):
		pBedMergeTest = pBedMerge.copy(tag = tag)
		pBedMergeTest.input = [infile]
		PyPPL(config).start(pBedMergeTest).run()
		self.assertFileEqual(pBedMergeTest.channel.get(), outfile)

	def dataProvider_testBedMerge2(self):
		infiles = [path.join(self.indir, 'bedmerge.bed')]*3
		outfile = path.join(self.outdir, 'bedmerge.out.bed')
		yield 't1', infiles, outfile

	def testBedMerge2(self, tag, infiles, outfile):
		pBedMerge2Test = pBedMerge2.copy(tag = tag)
		pBedMerge2Test.input = [infiles]
		PyPPL(config).start(pBedMerge2Test).run()
		self.assertFileEqual(pBedMerge2Test.channel.get(), outfile)


	# TODO:
	# - pBedMultiinter
	# - pBedShift
	# - pBedShuffle
	# - pBedSubtract
	# - pBedWindow
	# - pBedGenomecov
	# - pBedCluster
if __name__ == '__main__':
	testly.main()
