import sys, os, json, filecmp

sys.path.insert(0, "/home/m161047/tools/pyppl")
sys.path.insert(0, "/home/m161047/tools/bioprocs")

import pyppl, unittest, subprocess
from bioprocs.wxs import *
from bioprocs.gatk import *
from bioprocs.picard import *

# require samtools to compile bam files
def compfile(file1, proc, key = 'outfile', ignore = '@PG'):
	op   = __import__('gzip').open if file1.endswith (".gz") else open
	file2 = proc.output[key][0]
	ignores = [ignore] if not isinstance(ignore, list) else ignore
	def startsIgnore (line, ignores):
		for ig in ignores:
			if not line.startswith(ig):
				return False
		return True
	f1 = op(file1)
	f2 = op(file2)
	if file1.endswith('.bam'):
		p1 = subprocess.Popen(['samtools', 'view', file1], stdout=subprocess.PIPE)
		p2 = subprocess.Popen(['samtools', 'view', file2], stdout=subprocess.PIPE)
		f1 = p1.stdout
		f2 = p2.stdout
		
	while True:
		line1 = f1.readline()
		line2 = f2.readline()
		if line1 != line2 and not startsIgnore(line1, ignores) and not startsIgnore(line2, ignores):
			print "Line different:"
			print "[%s]:\n%s" % (file1, line1)
			print "[%s]:\n%s" % (file2, line2)
			return False
		if not line1 or not line2:
			break
	return True

class testWxs (unittest.TestCase):
	
	def testTrimmomaticPE (self):
		fqfile1   = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.1.fq.gz")
		fqfile2   = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.2.fq.gz")
		outfile1  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.1.clean.fq.gz")
		outfile2  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.2.clean.fq.gz")
		pTrimmomaticPE.input = {pTrimmomaticPE.input: [(fqfile1, fqfile2)]}
		pTrimmomaticPE.run()
		self.assertTrue (compfile(outfile1, pTrimmomaticPE, 'outfile1'))
		self.assertTrue (compfile(outfile2, pTrimmomaticPE, 'outfile2'))
	
	def testTrimmomaticSE (self):
		fqfile   = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.1.fq.gz")
		outfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticSE.1.clean.fq.gz")
		pTrimmomaticSE.input = {pTrimmomaticSE.input: [fqfile]}
		pTrimmomaticSE.run()
		self.assertTrue (compfile(outfile, pTrimmomaticSE))
		
	def testAlignPEByBWA (self):
		infile1  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.1.clean.fq.gz")
		infile2  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.2.clean.fq.gz")
		outfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.sam")
		pAlignPEByBWA.input = {pAlignPEByBWA.input: [(infile1, infile2)]}
		self.assertRaises(Exception, pAlignPEByBWA.run)
		# use ucsc hg19 reference genome
		pAlignPEByBWA.args['reffile'] = os.path.join ("testfiles", "wxs", "ref.fa")
		pAlignPEByBWA.run()
		self.assertTrue (compfile(outfile, pAlignPEByBWA))

		
	def testAlignSEByBWA (self):
		infile  = os.path.join ("testfiles", "wxs", "test.trimmomaticSE.1.clean.fq.gz")
		outfile = os.path.join ("testfiles", "wxs", "test.trimmomaticSE.1.sam")
		pAlignSEByBWA.input = {pAlignSEByBWA.input: [infile]}
		self.assertRaises(Exception, pAlignSEByBWA.run)
		pAlignSEByBWA.args['reffile'] = os.path.join ("testfiles", "wxs", "ref.fa")
		pAlignSEByBWA.run()
		self.assertTrue (compfile(outfile, pAlignSEByBWA))
		
	"""
	def testAlignPEByNGM (self):
		infile1  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.1.clean.fq.gz")
		infile2  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.2.clean.fq.gz")
		outfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.ngm.bam")
		pAlignPEByNGM.input = {pAlignPEByNGM.input: [(infile1, infile2)]}
		self.assertRaises(Exception, pAlignPEByNGM.run)
		# use ucsc hg19 reference genome
		pAlignPEByNGM.args['reffile'] = "/data2/junwenwang/shared/reference/hg19/ucsc_hg19.fa"
		#pAlignPEByNGM.args['nthread']  = 40
		pAlignPEByNGM.run()
		self.assertTrue (compfile(outfile, pAlignPEByNGM))

		
	def testAlignSEByNGM (self):
		infile  = os.path.join ("testfiles", "wxs", "test.trimmomaticSE.1.clean.fq.gz")
		outfile = os.path.join ("testfiles", "wxs", "test.trimmomaticSE.1.ngm.bam")
		pAlignSEByNGM.input = {pAlignSEByNGM.input: [infile]}
		self.assertRaises(Exception, pAlignSEByNGM.run)
		pAlignSEByNGM.args['reffile'] = "/data2/junwenwang/shared/reference/hg19/ucsc_hg19.fa"
		#pAlignPEByNGM.args['nthread']  = 40
		pAlignSEByNGM.run()
		self.assertTrue (compfile(outfile, pAlignSEByNGM))
	"""
	
	def testSortSam (self):
		infile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.sam")
		outfile = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.sorted.bam")
		pSortSam.input = {pSortSam.input: [infile]}
		pSortSam.run()
		self.assertTrue (compfile(outfile, pSortSam))
		
	def testMarkDup (self):
		infile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.sorted.bam")
		outfile = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.dedup.bam")
		pMarkDuplicates.input = {pMarkDuplicates.input: [infile]}
		pMarkDuplicates.run()
		self.assertTrue (compfile(outfile, pMarkDuplicates))
		
	def testIndexBam (self):
		infile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.dedup.bam")
		idxfile = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.dedup.bai")
		pIndexBam.input = {pIndexBam.input: [infile]}
		pIndexBam.run()
		self.assertTrue (filecmp.cmp(idxfile, pIndexBam.output['outfile'][0][:-3] + 'bai'))
	
	def testRealignerTargetCreator (self):
		infile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.dedup.bam")
		outfile = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.dedup.indelRealigner.intervals")
		pRealignerTargetCreator.input = {pRealignerTargetCreator.input: [infile]}
		pRealignerTargetCreator.args['reffile'] = os.path.join ("testfiles", "wxs", "ref.fa")
		try:
			pRealignerTargetCreator.run ()
			self.assertTrue (compfile(outfile, pRealignerTargetCreator))
		except Exception:
			pCreateSequenceDictionary.input = {pCreateSequenceDictionary.input: [pRealignerTargetCreator.args['reffile']]}
			pCreateSequenceDictionary.args["params"] = 'O="%s"' % os.path.join ("testfiles", "wxs", "ref.dict") 
			pCreateSequenceDictionary.run()
			pRealignerTargetCreator.run ()
			self.assertTrue (compfile(outfile, pRealignerTargetCreator))
		
	def testIndelRealigner (self):
		bamfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.dedup.bam")
		intfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.dedup.indelRealigner.intervals")
		outfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.realigned.bam")
		pIndelRealigner.input = {pIndelRealigner.input: [(bamfile, intfile)]}
		pIndelRealigner.args['reffile'] = os.path.join ("testfiles", "wxs", "ref.fa")
		pIndelRealigner.run ()
		self.assertTrue (compfile(outfile, pIndelRealigner))
		
	def testBaseRecalibrator (self):
		bamfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.realigned.bam")
		outfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.realigned.recal.table")
		pBaseRecalibrator.input = {pBaseRecalibrator.input: [bamfile]}
		pBaseRecalibrator.args['reffile']    = os.path.join ("testfiles", "wxs", "ref.fa")
		pBaseRecalibrator.args['knownSites'] = os.path.join ("testfiles", "wxs", "dbsnp.chrM.vcf")
		pBaseRecalibrator.run ()
		self.assertTrue (compfile(outfile, pBaseRecalibrator))
	
	def testPrintReads (self):
		bamfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.realigned.bam")
		rcfile   = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.realigned.recal.table")
		outfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.bam")
		pPrintReads.input = {pPrintReads.input: [(bamfile, rcfile)]}
		pPrintReads.args['reffile'] = os.path.join ("testfiles", "wxs", "ref.fa")
		pPrintReads.run ()
		self.assertTrue (compfile(outfile, pPrintReads))
	
	def testHaplotypeCaller (self):
		bamfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.bam")
		outfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.vcf")
		pHaplotypeCaller.input = {pHaplotypeCaller.input: [bamfile]}
		pHaplotypeCaller.args['reffile'] = os.path.join ("testfiles", "wxs", "ref.fa")
		pHaplotypeCaller.run ()
		self.assertTrue (compfile(outfile, pHaplotypeCaller, 'outfile', '##GATKCommandLine'))
	
	def testSelectVariants (self):
		vcffile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.vcf")
		outfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.selected.vcf")
		pSelectVariants.input = {pSelectVariants.input: [vcffile]}
		pSelectVariants.args['reffile'] = os.path.join ("testfiles", "wxs", "ref.fa")
		pSelectVariants.run ()
		self.assertTrue (compfile(outfile, pSelectVariants, 'outfile', '##GATKCommandLine'))
	
	def testVariantFiltration (self):
		vcffile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.selected.vcf")
		outfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.filtered.vcf")
		pVariantFiltration.input = {pVariantFiltration.input: [vcffile]}
		pVariantFiltration.args['reffile'] = os.path.join ("testfiles", "wxs", "ref.fa")
		pVariantFiltration.run ()
		self.assertTrue (compfile(outfile, pVariantFiltration, 'outfile', '##GATKCommandLine'))
	
	def testCNVnator (self):
		infile   = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.bam")
		outfile  = os.path.join ("testfiles", "wxs", "test.trimmomaticPE.cnv.vcf")
		pCNVnator.input = {pCNVnator.input: [infile]}
		pCNVnator.args['chrom']   = "chrM"
		pCNVnator.args['chrdir']  = os.path.join ("testfiles", "wxs")
		pCNVnator.run()
		self.assertTrue (compfile(outfile, pCNVnator))
		
		
if __name__ == '__main__':
	unittest.main()
		