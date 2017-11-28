import unittest
from pyppl import PyPPL, Proc
from helpers import getfile, procOK, config
from bioprocs.fastx import pFastqSim
from bioprocs.vcf import pVcfFilter, pVcfAnno, pVcfSplit, pVcf2Maf, pVcfMerge
from bioprocs.vcfnext import pVcfStatsPlot, pCallRate
from bioprocs.web import pDownloadGet
from bioprocs.common import pFiles2Dir
from bioprocs.tabix import pTabix, pTabixIndex


class TestVcf (unittest.TestCase):
	
	def test2_pVcfFilter (self):
		
		pVcfFilter.input        = [getfile('example.vcf')]
		pVcfFilter.forks        = 2
		
		pVcfFilter1              = pVcfFilter.copy()
		pVcfFilter2              = pVcfFilter.copy()
		pVcfFilter3              = pVcfFilter.copy()
		pVcfFilter4              = pVcfFilter.copy()
		pVcfFilter5              = pVcfFilter.copy()

		pVcfFilter1.args.filters = 'lambda record, samples: record.ID == "rs6054257"'
		
		pVcfFilter2.args.filters = 'lambda record, samples: record.ID == "rs6054257"'
		pVcfFilter2.args.keep    = False
		
		pVcfFilter3.args.keep     = False
		pVcfFilter3.args.filters = 'lambda record, samples: record.QUAL > 50'
		
		pVcfFilter4.args.keep     = False
		pVcfFilter4.args.filters = 'lambda record, samples: not record.FILTER or "PASS" in record.FILTER'
		
		pVcfFilter5.args.keep     = False
		pVcfFilter5.args.filters = 'lambda record, samples: all([sample.data.GQ >= 48 for sample in samples])'
		
		PyPPL().start([
			pVcfFilter1,
			pVcfFilter2,
			pVcfFilter3,
			pVcfFilter4,
			pVcfFilter5
		]).run()
		procOK(pVcfFilter1, 'example-rs6054257-keep.vcf', self)
		procOK(pVcfFilter2, 'example-rs6054257-nokeep.vcf', self)
		procOK(pVcfFilter3, 'example-q50.vcf', self)
		procOK(pVcfFilter4, 'example-pass.vcf', self)
		procOK(pVcfFilter5, 'example-gq48.vcf', self)

	def test2_pVcfAnno (self):
		pVcfAnno.input             = [getfile('example.vcf')]
		pVcfAnno.forks             = 2
		pVcfAnno1                  = pVcfAnno.copy()
		pVcfAnno2                  = pVcfAnno.copy()
		pVcfAnno3                  = pVcfAnno.copy()
		pVcfAnno1.args.tool        = 'snpeff'
		pVcfAnno2.args.tool        = 'vep'
		pVcfAnno2.args.params.port = 3337
		pVcfAnno2.args.genome      = 'GRCh37'
		pVcfAnno3.args.tool        = 'annovar'
		pVcfAnno3.args.gz          = True
		
		PyPPL().start([
			pVcfAnno1,
			pVcfAnno2,
			pVcfAnno3
		]).run()
		procOK(pVcfAnno1, 'example.snpeff.vcf', self)
		procOK(pVcfAnno2, 'example.vep.vcf', self)
		procOK(pVcfAnno3, 'example.annovar.vcf.gz', self)
		
	def test3_pCallRate (self):
		pFiles2Dir1 = pFiles2Dir.copy('vcfstats')
		pFiles2Dir1.input = [[getfile('example.vcf')] * 4]
		pCallRate1 = pCallRate.copy()
		pCallRate1.depends =  pFiles2Dir1
		
		PyPPL().start(pFiles2Dir1).run()
		procOK(pCallRate1, 'callrate', self)
		
	def test4_pPlotStats(self):

		pVcfAnno4                  = pVcfAnno.copy()
		pVcfAnno4.args.tool        = 'snpeff'
		pVcfAnno4.args.snpeffStats = True
		
		pFiles2Dir2         = pFiles2Dir.copy()
		pFiles2Dir2.depends = pVcfAnno4
		pFiles2Dir2.tag     = 'plotstats'
		pFiles2Dir2.input   = lambda ch: [ch.colAt(1).expand(0, '*.stats.csv').flatten()*4]
		
		pVcfStatsPlot2         = pVcfStatsPlot.copy()
		pVcfStatsPlot2.depends = pFiles2Dir2
		pVcfStatsPlot2.expect  = 'ls {{out.outdir}}/mats/*.mat.txt'
		PyPPL().start(pVcfAnno4).run()

	def test5_pVcfSplit(self):
		pTabix.input = (
			"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz", 
			getfile('region.txt')
		)

		pVcfSplit1 = pVcfSplit.copy()
		pVcfSplit1.depends = pTabix
		pVcfSplit1.expect  = "ls -l {{out.outdir}}/HG00098.vcf {{out.outdir}}/NA20828.vcf"

		pVcfSplit3 = pVcfSplit.copy()
		pVcfSplit3.depends = pTabix
		pVcfSplit3.input   = lambda ch: ch.cbind(','.join(['NA' + str(i) for i in range(19074, 19190)]))
		pVcfSplit3.expect  = "ls -l {{out.outdir}}/NA19189.vcf {{out.outdir}}/NA19119.vcf {{out.outdir}}/NA19074.vcf"
		pVcfSplit3.args.nthread = 20

		# dont have chr prefix on chromosomes!
		pVcfSplit4 = pVcfSplit.copy()
		# gatk is toooooo slow for test
		#pVcfSplit4.depends = pTabix 
		pVcfSplit4.input   = lambda ch: ch.cbind(','.join(['NA' + str(i) for i in range(19074, 19140)]))
		pVcfSplit4.expect  = "ls -l {{out.outdir}}/NA19189.vcf {{out.outdir}}/NA19119.vcf"
		pVcfSplit4.args.nthread = 20
		pVcfSplit4.args.tool  = 'gatk'

		pVcfSplit.tag = 'sample'
		pVcfSplit.depends = pTabix
		pVcfSplit.input   = lambda ch: ch.cbind("NA19660,NA19917")
		pVcfSplit.expect  = "ls -l {{out.outdir}}/NA19660.vcf {{out.outdir}}/NA19917.vcf"

		PyPPL().start(pTabix).run()

	def test6_pPlotStats_Real(self):
		
		pVcfAnno5 = pVcfAnno.copy()
		pVcfAnno5.input = [getfile('example.vcf')]
		pVcfAnno5.args.snpeffStats = True

		pFiles2Dir3 = pFiles2Dir.copy()
		pFiles2Dir3.depends = pVcfAnno5
		pFiles2Dir3.input = lambda ch: [ch.colAt(1).expand(0, '*.stats*.csv').flatten()]

		pVcfStatsPlot3 = pVcfStatsPlot.copy()
		pVcfStatsPlot3.depends = pFiles2Dir3
		pVcfStatsPlot3.expect  = 'ls {{out.outdir}}/mats/*.mat.txt'

		PyPPL().start(pVcfAnno5).run()

	def test7_pVcf2Maf(self):
		pFiles2Dir4 = pFiles2Dir.copy()

		pTabixIndex.input      = pVcfSplit.channel.outdir.expand(pattern = '*.vcf')
		pTabixIndex.forks      = 10
		pFiles2Dir4.depends    = pTabixIndex
		pFiles2Dir4.input      = lambda ch: [ch.flatten()]
		pVcfMerge.depends      = pFiles2Dir4

		pAddChr2         = Proc()
		pAddChr2.depends = pVcfMerge
		pAddChr2.input   = "infile:file"
		pAddChr2.output  = "outfile:file:{{in.infile|bn}}"
		pAddChr2.script  = """awk '$0 ~ /^#/ {print} $0 !~ /^#/ {print "chr"$0}' {{in.infile}} > {{out.outfile}}"""

		pVcf2Maf.depends       = pAddChr2
		PyPPL().start(pTabixIndex).run()
		procOK(pVcf2Maf, 'vcf2maf.maf', self)
		
if __name__ == '__main__':
	unittest.main(failfast=True)