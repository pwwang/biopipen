import helpers, unittest
from pyppl import Box
from os import path, makedirs
from bioprocs.utils.snp import snpinfo
from medoo import Medoo

class TestUtilsSnp(helpers.TestCase):
	
	def dataProvider_testSnpinfo(self, testdir, indir):
		infile = path.join(indir, 'snps.txt')
		yield infile, None, 'skip', 'hg19', '150', {'cnames': ['Snp', 'Index']}, {}, None, testdir, {
			'rs568942921': Box(chrom='chr1',chromStart=15073281L,chromEnd=15073282L,name='rs568942921',score=0,strand='+',refNCBI='C',refUCSC='C',avHet=0.0,avHetSE=0.0,func='|intron',alleleFreqCount=0,alleles='',alleleNs='',alleleFreqs=''),
			'rs529435559': Box(chrom='chr1',chromStart=41287734L,chromEnd=41287735L,name='rs529435559',score=0,strand=u'+',refNCBI='T',refUCSC='T',avHet=0.0,avHetSE=0.0,func='|intron',alleleFreqCount=0,alleles='',alleleNs='',alleleFreqs=''), 
			'rs868324925': Box(chrom='chr1',chromStart=28180487L,chromEnd=28180488L,name='rs868324925',score=0,strand=u'+',refNCBI='G',refUCSC='G',avHet=0.0,avHetSE=0.0,func='|unknown',alleleFreqCount=0,alleles='',alleleNs='',alleleFreqs=''), 
			'rs775809821': Box(chrom='chr1',chromStart=10019L,chromEnd=10020L,name='rs775809821',score=0,strand=u'+',refNCBI='A',refUCSC='A',avHet=0.0,avHetSE=0.0,func='|near-gene-5',alleleFreqCount=0,alleles='',alleleNs='',alleleFreqs=''), 
			'rs764544638': Box(chrom='chr1',chromStart=1966084L,chromEnd=1966085L,name='rs764544638',score=0,strand=u'+',refNCBI='A',refUCSC='A',avHet=0.0,avHetSE=0.0,func='|unknown',alleleFreqCount=0,alleles='',alleleNs='',alleleFreqs='')
		}
	
	def testSnpinfo(self, infile, outfile, notfound, genome, dbsnpver, inopts, outopts, snpcol, cachedir, output):
		self.maxDiff = None		
		ret = snpinfo(infile, outfile, notfound, genome, dbsnpver, inopts, outopts, snpcol, cachedir)
		for key, val in output.items():
			for k, v in val.items():
				self.assertEqual(getattr(ret[key], k), v)
				
	def dataProvider_testSnpinfoCached(self, testdir, indir):
		infile = path.join(indir, 'snps.txt')
		cachedir = path.join(testdir, 'cached')
		makedirs(cachedir)
		snpinfo(infile, None, 'skip', 'hg19', '150', {'cnames': ['Snp', 'Index']}, {}, None, cachedir)
		
		yield infile, None, 'skip', 'hg19', '150', {'cnames': ['Snp', 'Index']}, {}, None, cachedir, {
			'rs568942921': Box(chrom='chr1',chromStart=15073281L,chromEnd=15073282L,name='rs568942921',score=0,strand='+',refUCSC='C',alleleFreqCount=0,alleles=''),
			'rs529435559': Box(chrom='chr1',chromStart=41287734L,chromEnd=41287735L,name='rs529435559',score=0,strand=u'+',refUCSC='T',alleleFreqCount=0,alleles=''), 
			'rs868324925': Box(chrom='chr1',chromStart=28180487L,chromEnd=28180488L,name='rs868324925',score=0,strand=u'+',refUCSC='G',alleleFreqCount=0,alleles=''), 
			'rs775809821': Box(chrom='chr1',chromStart=10019L,chromEnd=10020L,name='rs775809821',score=0,strand=u'+',refUCSC='A',alleleFreqCount=0,alleles=''), 
			'rs764544638': Box(chrom='chr1',chromStart=1966084L,chromEnd=1966085L,name='rs764544638',score=0,strand=u'+',refUCSC='A',alleleFreqCount=0,alleles='')
		}
	
	def testSnpinfoCached(self, infile, outfile, notfound, genome, dbsnpver, inopts, outopts, snpcol, cachedir, output):
		self.maxDiff = None		
		ret = snpinfo(infile, outfile, notfound, genome, dbsnpver, inopts, outopts, snpcol, cachedir)
		for key, val in output.items():
			for k, v in val.items():
				self.assertEqual(getattr(ret[key], k), v)
				
	def testSnpinfo(self, infile, outfile, notfound, genome, dbsnpver, inopts, outopts, snpcol, cachedir, output):
		self.maxDiff = None		
		ret = snpinfo(infile, outfile, notfound, genome, dbsnpver, inopts, outopts, snpcol, cachedir)
		for key, val in output.items():
			for k, v in val.items():
				self.assertEqual(getattr(ret[key], k), v)
				
	def dataProvider_testSnpinfoOutfile(self, testdir, indir, outdir):
		infile = path.join(indir, 'snps.txt')
		cachedir = path.join(testdir, 'cached')
		makedirs(cachedir)
		snpinfo(infile, None, 'skip', 'hg19', '150', {'cnames': ['Snp', 'Index']}, {}, None, cachedir)
		
		outfile = path.join(testdir, 'snpinfo.txt')
		exptfile = path.join(outdir, 'snpinfo.txt')
		
		yield infile, outfile, 'skip', 'hg19', '150', {'cnames': ['Snp', 'Index']}, {}, None, cachedir, exptfile
	
	def testSnpinfoOutfile(self, infile, outfile, notfound, genome, dbsnpver, inopts, outopts, snpcol, cachedir, exptfile):
		self.maxDiff = None		
		snpinfo(infile, outfile, notfound, genome, dbsnpver, inopts, outopts, snpcol, cachedir)
		self.assertFileEqual(outfile, exptfile)
		
		
if __name__ == '__main__':
	unittest.main(verbosity = 2)