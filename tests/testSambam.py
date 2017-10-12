import os, sys, unittest, shutil, addPath

from pyppl import PyPPL, Channel, Proc, Box
from bioprocs.fastx import pFastqSim, pFastq2Sam, pFastqSE2Sam
from bioprocs.sambam import pSam2Bam, pBamMarkdup, pBamRecal, pBamReadGroup, pBamReorder, pBamMerge, pBam2Gmut, pBamPair2Smut, pBam2Cnv, pBamStats, pBam2Fastq, pBam2FastqSE, pBam2Counts
from bioprocs.common import pFiles2List
from bioprocs.web import pDownloadGet
from bioprocs import params

params = params.toDict()
unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)
class TestSambam (unittest.TestCase):
	
	data = Box()
	
	def test0_getRef (self):
		pDownloadGet.tag      = 'lambdaPhage'
		pDownloadGet.input    = ['https://raw.githubusercontent.com/BenLangmead/bowtie2/master/example/reference/lambda_virus.fa']
		pChangeChrom          = Proc(desc = 'Change the chrom to regular name.')
		pChangeChrom.depends  = pDownloadGet
		pChangeChrom.input    = 'infile:file'
		pChangeChrom.output   = 'outfile:file:{{in.infile | bn}}'
		pChangeChrom.lang     = 'python'
		pChangeChrom.script   = """#
		## indent remove ##
		from random import choice
		open("{{out.outfile}}", "w").write(
			open("{{in.infile}}").read().replace(">gi|9626243|ref|NC_001416.1| Enterobacteria phage lambda, complete genome", ">chr1") + 
			">chr2\\n" +
			"\\n".join([''.join([choice(['A', 'T', 'G', 'C']) for x in range(70)]) for y in range(500)]) +
			"\\n"
		)
		"""
		pChangeChrom.callback = lambda p: setattr(params, 'ref', p.channel[0][0])
		PyPPL().start(pDownloadGet).run()
			
	def test1_pFastq2Sam (self):
		self.data.fqs = []
		
		pFastqSimWgsim  = pFastqSim.copy ()
		pFastqSimDwgsim = pFastqSim.copy ()
		
		pFastqSimWgsim.input       = [1]
		pFastqSimDwgsim.input      = [1]
		pFastqSimWgsim.args.ref    = params.ref
		pFastqSimDwgsim.args.ref   = params.ref
		pFastqSimWgsim.args.num    = 100000
		pFastqSimDwgsim.args.num   = 100000
		pFastqSimWgsim.args.tool   = 'wgsim'
		pFastqSimDwgsim.args.tool  = 'dwgsim'
		pFastqSimWgsim.callback    = lambda p: self.data.fqs.extend(p.channel)
		pFastqSimDwgsim.callback   = lambda p: self.data.fqs.extend(p.channel)
		
		PyPPL().start(pFastqSimWgsim, pFastqSimDwgsim).run()
		
	def test2_pFastq2Sam (self):
		pFastq2Sam.args.nthread     = 2
		pFastq2Sam.forks            = 4
		pFastq2Sam.input            = self.data.fqs
		pFastq2Sam.args.ref         = params.ref
		pFastq2SamBwa               = pFastq2Sam.copy('bwa')
		pFastq2SamNgm               = pFastq2Sam.copy('ngm')
		pFastq2SamBowtie2           = pFastq2Sam.copy('bowtie2')
		pFastq2SamBwa.args.tool     = 'bwa'
		pFastq2SamNgm.args.tool     = 'ngm'
		pFastq2SamBowtie2.args.tool = 'bowtie2'
		pFastq2SamBowtie2.callback  = lambda p: self.data.update({'sams': p.channel})
		
		PyPPL().start(pFastq2SamBwa, pFastq2SamNgm, pFastq2SamBowtie2).run()
	
	def test3_Sam2BamSambamba (self):
		pSam2BamSbb1  = pSam2Bam.copy('sambamba')
		pSam2BamSbb2  = pSam2Bam.copy('sambamba')
		pSam2BamSbb3  = pSam2Bam.copy('sambamba')
		pSam2BamSbb1.args.tool    = 'sambamba'
		pSam2BamSbb1.args.markdup = True
		pSam2BamSbb2.args.tool    = 'sambamba'
		pSam2BamSbb2.args.sort    = False
		pSam2BamSbb2.args.index   = False
		pSam2BamSbb2.args.markdup = False
		pSam2BamSbb2.args.rmdup   = False
		pSam2BamSbb3.args.tool    = 'sambamba'
		pSam2BamSbb3.args.sort    = True
		pSam2BamSbb3.args.index   = True
		pSam2BamSbb3.args.markdup = False
		pSam2BamSbb3.args.rmdup   = False
		
		PyPPL().start(pSam2BamSbb1, pSam2BamSbb2, pSam2BamSbb3).run()
			
	def test3_Sam2BamSamtools (self):
		pSam2BamSmt1  = pSam2Bam.copy('samtools')
		pSam2BamSmt2  = pSam2Bam.copy('samtools')
		pSam2BamSmt3  = pSam2Bam.copy('samtools')
		pSam2BamSmt1.args.tool    = 'samtools'
		pSam2BamSmt1.args.markdup = True
		pSam2BamSmt2.args.tool    = 'samtools'
		pSam2BamSmt2.args.sort    = False
		pSam2BamSmt2.args.index   = False
		pSam2BamSmt2.args.markdup = False
		pSam2BamSmt2.args.rmdup   = False
		pSam2BamSmt3.args.tool    = 'samtools'
		pSam2BamSmt3.args.sort    = True
		pSam2BamSmt3.args.index   = True
		pSam2BamSmt3.args.markdup = False
		pSam2BamSmt3.args.rmdup   = False
		
		PyPPL().start(pSam2BamSmt1, pSam2BamSmt2, pSam2BamSmt3).run()
			
	def test3_Sam2BamPicard (self):
		pSam2BamPcd1  = pSam2Bam.copy('picard')
		pSam2BamPcd2  = pSam2Bam.copy('picard')
		pSam2BamPcd3  = pSam2Bam.copy('picard')
		pSam2BamPcd1.args.tool    = 'picard'
		pSam2BamPcd1.args.markdup = True
		pSam2BamPcd2.args.tool    = 'picard'
		pSam2BamPcd2.args.sort    = False
		pSam2BamPcd2.args.index   = False
		pSam2BamPcd2.args.markdup = False
		pSam2BamPcd2.args.rmdup   = False
		pSam2BamPcd3.args.tool    = 'picard'
		pSam2BamPcd3.args.sort    = True
		pSam2BamPcd3.args.index   = True
		pSam2BamPcd3.args.markdup = False
		pSam2BamPcd3.args.rmdup   = False

		PyPPL().start(pSam2BamPcd1, pSam2BamPcd2, pSam2BamPcd3).run()

	
	def test3_1Sam2BamBiobambam (self):
		pSam2Bam.forks= 2
		pSam2Bam.input= self.data.sams
		pSam2BamBbb1  = pSam2Bam.copy('biobambam')
		pSam2BamBbb2  = pSam2Bam.copy('biobambam')
		pSam2BamBbb3  = pSam2Bam.copy('biobambam')
		pSam2BamBbb1.args.params.maxreadlen = 5000
		pSam2BamBbb2.args.params.maxreadlen = 5000
		pSam2BamBbb3.args.params.maxreadlen = 5000
		pSam2BamBbb1.args.tool    = 'biobambam'
		pSam2BamBbb1.args.markdup = True
		pSam2BamBbb2.args.tool    = 'biobambam'
		pSam2BamBbb2.args.sort    = False
		pSam2BamBbb2.args.index   = False
		pSam2BamBbb2.args.markdup = False
		pSam2BamBbb2.args.rmdup   = False
		pSam2BamBbb3.args.tool    = 'biobambam'
		pSam2BamBbb3.args.sort    = True
		pSam2BamBbb3.args.index   = True
		pSam2BamBbb3.args.markdup = False
		pSam2BamBbb3.args.rmdup   = False
		pSam2BamBbb3.callback     = lambda p: self.data.update({'bams': p.channel.colAt(0)})
		
		PyPPL().start(pSam2BamBbb1, pSam2BamBbb2, pSam2BamBbb3).run()
		
	def test4_BamMarkdup (self):
		pBamMarkdup.args.nthread = 2
		pBamMarkdup.forks = 4
		pBamMarkdup.input = self.data.bams
		pBamMarkdup1 = pBamMarkdup.copy ('sambamba')
		pBamMarkdup2 = pBamMarkdup.copy ('biobambam')
		pBamMarkdup3 = pBamMarkdup.copy ('picard')
		pBamMarkdup4 = pBamMarkdup.copy ('samtools')
		pBamMarkdup5 = pBamMarkdup.copy ('bamutil')
		pBamMarkdup1.args.tool = 'sambamba'
		pBamMarkdup2.args.tool = 'biobambam'
		pBamMarkdup3.args.tool = 'picard'
		pBamMarkdup4.args.tool = 'samtools'
		pBamMarkdup5.args.tool = 'bamutil'
		
		PyPPL().start(pBamMarkdup1,pBamMarkdup2,pBamMarkdup4,pBamMarkdup5).run()
				
	def test2_FastqSE2Sam (self):
		pFastqSE2Sam.args.nthread = 2
		pFastqSE2Sam.forks = 8
		pFastqSE2Sam.input = Channel(self.data.fqs).fold(1)
		pFastqSE2Sam.args.ref = params.ref
		pFastqSE2Sam1 = pFastqSE2Sam.copy('bwa')
		pFastqSE2Sam2 = pFastqSE2Sam.copy('ngm')
		pFastqSE2Sam3 = pFastqSE2Sam.copy('bowtie2')
		pFastqSE2Sam1.args.tool = 'bwa'
		pFastqSE2Sam2.args.tool = 'ngm'
		pFastqSE2Sam3.args.tool = 'bowtie2'
		
		PyPPL().start(pFastqSE2Sam1, pFastqSE2Sam2, pFastqSE2Sam3).run()
	
	
	def test4_BamRecal (self):
		pBamRecalTest = pBamRecal.copy()
		pBamRecalTest.input = self.data.bams
		pBamRecalTest.args.ref = params.ref
		pBamRecalTest.args.knownSites = params.dbsnp
		pBamRecalTest.args.paramsBaseRecalibrator.maximum_cycle_value = 5000
		pBamRecalTest.forks = 2
		
		PyPPL().start(pBamRecalTest).run()
		
	def test4_BamReadGroup (self):
		pBamReadGroup1 = pBamReadGroup.copy()
		pBamReadGroup2 = pBamReadGroup.copy()
		pBamReadGroup1.input = self.data.bams
		pBamReadGroup2.input = self.data.bams
		pBamReadGroup2.args.tool = 'bamutil'
		pBamReadGroup1.forks = 10
		pBamReadGroup2.forks = 10
		
		PyPPL().start(pBamReadGroup1, pBamReadGroup2).run()
	
	def test4_BamReorder (self):
		
		pBamReorder.args.ref   = params.ref
		pBamReorder.forks = 10
		pBamReorder.input      = self.data.bams
		
		PyPPL().start(pBamReorder).run()
	
	def test4_Mergebams (self):
		pFiles2List.input = [Channel.create(self.data.bams).flatten()]
		
		pBamMerge1 = pBamMerge.copy()
		pBamMerge2 = pBamMerge.copy()
		pBamMerge3 = pBamMerge.copy()
		pBamMerge4 = pBamMerge.copy()
		
		pBamMerge1.depends = pFiles2List
		pBamMerge2.depends = pFiles2List
		pBamMerge3.depends = pFiles2List
		pBamMerge4.depends = pFiles2List
		
		pBamMerge1.args.tool = 'bamutil'
		pBamMerge2.args.tool = 'samtools'
		pBamMerge3.args.tool = 'picard'
		pBamMerge4.args.tool = 'sambamba'
		
		PyPPL().start(pFiles2List).run()
		
	def test4_pBam2Gmut (self):
		pBam2Gmut.args.ref = params.ref
		pBam2Gmut.args.nthread = 16
		pBam2Gmut.forks    = 10
		pBam2Gmut.input    = self.data.bams
		
		pBam2Gmut1         = pBam2Gmut.copy()
		pBam2Gmut2         = pBam2Gmut.copy()
		pBam2Gmut3         = pBam2Gmut.copy()
		pBam2Gmut4         = pBam2Gmut.copy()
		pBam2Gmut5         = pBam2Gmut.copy()

		pBam2Gmut1.args.tool = 'gatk'
		pBam2Gmut2.args.tool = 'snvsniffer'
		pBam2Gmut3.args.tool = 'platypus'
		pBam2Gmut4.args.tool = 'vardict'
		pBam2Gmut5.args.tool = 'strelka'
		
		starts = [
			pBam2Gmut1, 
			pBam2Gmut2, 
			pBam2Gmut3, 
			#pBam2Gmut4, # too slow for test
			pBam2Gmut5
		]
		
		PyPPL().start(*starts).run()
		
	def test4_pBamPair2Smut (self):
		pBamPair2Smut.input        = Channel.create(self.data.bams).unfold(2)
		pBamPair2Smut.args.ref     = params.ref
		pBamPair2Smut.args.nthread = 16
		pBamPair2Smut.forks        = 10
		
		pBamPair2Smut1 = pBamPair2Smut.copy()
		pBamPair2Smut2 = pBamPair2Smut.copy()
		pBamPair2Smut3 = pBamPair2Smut.copy()
		pBamPair2Smut4 = pBamPair2Smut.copy()
		pBamPair2Smut5 = pBamPair2Smut.copy()
		pBamPair2Smut6 = pBamPair2Smut.copy()
		
		pBamPair2Smut1.args.tool = 'gatk'
		pBamPair2Smut2.args.tool = 'somaticsniper'
		pBamPair2Smut3.args.tool = 'snvsniffer'
		pBamPair2Smut4.args.tool = 'strelka'
		pBamPair2Smut5.args.tool = 'virmid'
		pBamPair2Smut6.args.tool = 'vardict'
		pBamPair2Smut4.args.gz   = True
		
		starts = [
			pBamPair2Smut1, 
			pBamPair2Smut2, 
			pBamPair2Smut3, 
			pBamPair2Smut4, 
			#pBamPair2Smut5, # itself has error for test
			#pBamPair2Smut6, # too slow for test
		]
		PyPPL().start(*starts).run()
		
	def test5_pBam2Cnv (self):
		pBam2Cnv.input = self.data.bams
		
		pBam2Cnv1 = pBam2Cnv.copy()
		pBam2Cnv1.forks  = 10
		pBam2Cnv1.args.tool  = 'cnvkit'
		pBam2Cnv1.args.ref   = params.ref
		
		pBam2Cnv2 = pBam2Cnv.copy()
		pBam2Cnv2.forks = 10
		pBam2Cnv2.args.tool = 'cnvnator'
		pBam2Cnv2.args.ref   = params.ref
		
		pBam2Cnv3 = pBam2Cnv.copy()
		pBam2Cnv3.forks = 10
		pBam2Cnv3.args.tool = 'wandy'
		pBam2Cnv3.args.wandy = '/data2/bsi/RandD/s115463.Aneuploidy/SinlgeSampleWandy/Wandy/Wandy'
		
		starts = [
			pBam2Cnv1,
			pBam2Cnv2,
			pBam2Cnv3,
		]
		PyPPL().start(*starts).run()
		
	def test6_pBamStats (self):
		pBamStats.input = self.data.bams
		pBamStats.forks = 10
		
		pBamStats1 = pBamStats.copy()
		pBamStats1.args.plot = True
		
		starts = [
			pBamStats1,
		]
		PyPPL().start(*starts).run()
		
	def test7_pBam2Fastq (self):
		
		pBam2FastqBmm       = pBam2Fastq.copy(newid = 'pBam2FastqTest', tag = 'biobambam')
		pBam2FastqBmm.input = self.data.bams
		pBam2FastqBmm.forks = 10
		
		pBam2FastqBtl       = pBam2FastqBmm.copy(newid = 'pBam2FastqTest', tag = 'bedtools')
		pBam2FastqBtl.args.tool = 'bedtools'
		
		pBam2FastqStl       = pBam2FastqBmm.copy(newid = 'pBam2FastqTest', tag = 'samtools')
		pBam2FastqStl.args.tool = 'samtools'
		
		pBam2FastqPcd       = pBam2FastqBmm.copy(newid = 'pBam2FastqTest', tag = 'picard')
		pBam2FastqPcd.args.tool = 'picard'
		
		PyPPL().start('pBam2FastqTest').run()
		
	def test7_pBam2FastqSE (self):
		
		pBam2FastqSEBmm       = pBam2FastqSE.copy(newid = 'pBam2FastqSETest', tag = 'biobambam')
		pBam2FastqSEBmm.input = self.data.bams
		pBam2FastqSEBmm.forks = 10
		
		pBam2FastqSEBtl       = pBam2FastqSEBmm.copy(newid = 'pBam2FastqSETest', tag = 'bedtools')
		pBam2FastqSEBtl.args.tool = 'bedtools'
		
		pBam2FastqSEStl       = pBam2FastqSEBmm.copy(newid = 'pBam2FastqSETest', tag = 'samtools')
		pBam2FastqSEStl.args.tool = 'samtools'
		
		#pBam2FastqSEPcd       = pBam2FastqSEBmm.copy(newid = 'pBam2FastqSETest', tag = 'picard')
		#pBam2FastqSEPcd.args.tool = 'picard'
		
		PyPPL().start('pBam2FastqSETest').run()

	def test8_pBam2Counts(self):
		pGenerateGTF = Proc()
		pGenerateGTF.input = {'in': 0}
		pGenerateGTF.output = "outfile:file:gene.gtf"
		pGenerateGTF.lang = "python"
		pGenerateGTF.script = """#
		## indent remove ##
		import random
		random.seed(0)
		with open("{{out.outfile}}", "w") as f:
			for i, chr in enumerate(['chr1', 'chr2']):
				gstart = 0
				for j in list(range(100)):
					gene = 'GENE' + str(i*100 + j+1)
					glen = random.randint(50, 200)
					gap  = random.randint(100, 300)
					tmp  = [chr, 'TEST', 'transcript', str(gstart + gap), str(gstart + gap + glen), '.', random.choice(['+', '-']), '.', 'gene_id "%s"; gene_name "%s"' % (gene, gene)]
					f.write("\\t".join(tmp) + "\\n")
					gstart += gap + glen
		"""

		pBam2Counts.depends = pGenerateGTF
		pBam2Counts.input = self.data.bams
		pBam2Counts.forks = 10
		pBam2Counts.callfront = lambda p: p.args.update({'refgene':pGenerateGTF.channel.get()})
		PyPPL().start(pGenerateGTF).run()
		
if __name__ == '__main__':
	unittest.main(failfast=True)