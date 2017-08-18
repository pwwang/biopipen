import os, sys, unittest, addPath

from pyppl import pyppl, doct, channel, proc
from bioprocs.fastx import pFastqPESim
from bioprocs.vcf import pVcfFilter, pVcfAnno, pCallRate
from bioprocs.vcfnext import pStats2Matrix, pPlotStats
from bioprocs.web import pDownloadGet
from bioprocs.common import pFiles2Dir
from bioaggrs import params

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

class TestVcf (unittest.TestCase):
	
	data = doct()
	def test0_prepareData(self):
		
		pDownloadGet.input = {pDownloadGet.input: ["https://raw.githubusercontent.com/jamescasbon/PyVCF/master/vcf/test/example-4.1.vcf"]}
		pAddChr = proc ()
		pAddChr.input = "infile:file"
		pAddChr.output = "outfile:file:{{infile | bn}}"
		pAddChr.depends = pDownloadGet
		pAddChr.script = """
		awk '$0 ~ /^#/ && $0 !~ /^##contig/ {print} $0 ~ /##contig/ {print "##contig=<ID=chr20,length=63025520>"} $0 !~ /^#/ {print "chr"$0}' "{{infile}}" > "{{outfile}}"
		"""
		pAddChr.callback = lambda p: self.data.update({'vcfs': p.channel})
		pyppl().starts(pDownloadGet).run()
		
	def test2_pVcfFilter (self):
		
		pVcfFilter.input        = {pVcfFilter.input: self.data.vcfs}
		pVcfFilter.forks        = 2
		
		pVcfFilter1             = pVcfFilter.copy()
		pVcfFilter2             = pVcfFilter.copy()
		pVcfFilter3             = pVcfFilter.copy()
		pVcfFilter4             = pVcfFilter.copy()
		pVcfFilter1.args.tool   = 'gatk'
		pVcfFilter2.args.tool   = 'bcftools'
		pVcfFilter3.args.tool   = 'snpsift'
		pVcfFilter4.args.tool   = 'vcflib'
		pVcfFilter5             = pVcfFilter1.copy()
		pVcfFilter6             = pVcfFilter2.copy()
		pVcfFilter7             = pVcfFilter3.copy()
		pVcfFilter8             = pVcfFilter4.copy()
		
		pVcfFilter5.args.ref    = params.ref
		pVcfFilter1.args.selectors.type     = "indel"
		pVcfFilter1.expect = 'grep 1230237 {{outfile}} && grep 1234567 {{outfile}} && !(grep 14370 {{outfile}}) && !(grep 17330 {{outfile}}) && !(grep 1110696 {{outfile}})'
		
		pVcfFilter2.args.selectors.filter   = "PASS"
		pVcfFilter2.expect = 'grep 1230237 {{outfile}} && grep 1234567 {{outfile}} && grep 14370 {{outfile}} && !(grep 17330 {{outfile}}) && grep 1110696 {{outfile}}'
		
		pVcfFilter3.args.selectors.genotype = {0: '0/0'}
		pVcfFilter3.expect = 'grep 1230237 {{outfile}} && !(grep 1234567 {{outfile}}) && grep 14370 {{outfile}} && grep 17330 {{outfile}} && !(grep 1110696 {{outfile}})'
		
		pVcfFilter4.args.selectors.type     = "snp"
		pVcfFilter4.args.selectors.genotype = {2: '1/1'}
		pVcfFilter4.expect = '!(grep 1230237 {{outfile}}) && !(grep 1234567 {{outfile}}) && grep 14370 {{outfile}} && !(grep 17330 {{outfile}}) && !(grep 1110696 {{outfile}})'
		
		pVcfFilter5.args.selectors.type     = "snp"
		pVcfFilter5.args.filters.lowQUAL    = "QUAL < 10"
		pVcfFilter5.args.filters.lowDP      = "DP < 10"
		pVcfFilter5.args.keep               = False
		pVcfFilter5.expect = '!(grep 1230237 {{outfile}}) && !(grep 1234567 {{outfile}}) && grep 14370 {{outfile}} && !(grep 17330 {{outfile}}) && grep 1110696 {{outfile}}'
		
		pVcfFilter6.args.selectors.type     = "snp"
		pVcfFilter6.args.filters.lowQUAL    = "QUAL < 10 | DP < 10"
		pVcfFilter6.expect = '!(grep 1230237 {{outfile}}) && !(grep 1234567 {{outfile}}) && grep 14370 {{outfile}} && !(grep 17330 {{outfile}}) && grep 1110696 {{outfile}}'
		
		pVcfFilter7.args.selectors.type     = "snp"
		pVcfFilter7.args.filters.lowQUAL    = "(QUAL < 10) | (DP < 10)"
		pVcfFilter7.args.keep               = False
		pVcfFilter7.expect = '!(grep 1230237 {{outfile}}) && !(grep 1234567 {{outfile}}) && grep 14370 {{outfile}} && !(grep 17330 {{outfile}}) && grep 1110696 {{outfile}}'
		
		pVcfFilter8.args.selectors.type     = "snp"
		pVcfFilter8.args.filters.lowQUAL    = "QUAL < 10"
		pVcfFilter8.args.filters.lowDP      = "DP < 10"
		pVcfFilter8.args.gz                 = True
		pVcfFilter8.expect = '!(zcat {{outfile}} | grep 1230237) && !(zcat {{outfile}} | grep 1234567) && zcat {{outfile}} | grep 14370 && !(zcat {{outfile}} | grep 17330) && zcat {{outfile}} | grep 1110696'
		
		starts = [
			pVcfFilter1,
			pVcfFilter2,
			pVcfFilter3,
			pVcfFilter4,
			pVcfFilter5,
			pVcfFilter6,
			pVcfFilter7,
			pVcfFilter8,
		]
		pyppl().starts(*starts).run()
	
	def test2_pVcfAnno (self):
		pVcfAnno.input        = {pVcfAnno.input: self.data.vcfs}
		pVcfAnno.forks        = 2
		pVcfAnno1             = pVcfAnno.copy()
		pVcfAnno2             = pVcfAnno.copy()
		pVcfAnno3             = pVcfAnno.copy()
		pVcfAnno1.args.tool   = 'snpeff'
		pVcfAnno2.args.tool   = 'vep'
		pVcfAnno2.args.dbpath = '/data2/junwenwang/shared/reference/hg19/vep/cache/'
		pVcfAnno2.args.params = '--port 3337'
		pVcfAnno2.args.genome = 'GRCh37'
		pVcfAnno3.args.tool   = 'annovar'
		pVcfAnno3.args.dbpath = '/data2/junwenwang/shared/tools/annovar/humandb/'
		pVcfAnno3.args.gz     = True
		pVcfAnno1.callback    = lambda p: self.data.update({'vcfs': p.channel.colAt(0)})
		
		starts = [
			pVcfAnno1,
			pVcfAnno2,
			pVcfAnno3
		]
		pyppl().starts(*starts).run()
		
	def test3_pCallRate (self):
		pFiles2Dir1 = pFiles2Dir.copy('vcfstats')
		pFiles2Dir1.input = {pFiles2Dir1.input: [self.data.vcfs.toList() * 4]}
		pCallRate1 = pCallRate.copy()
		pCallRate1.depends =  pFiles2Dir1
		
		pyppl().starts(pFiles2Dir1).run()
		
	def test4_pPlotStats(self):
		pVcfAnno4             = pVcfAnno.copy()
		pVcfAnno4.args.tool   = 'snpeff'
		pVcfAnno4.args.snpeffStats = True
		
		pFiles2Dir2 = pFiles2Dir.copy()
		pFiles2Dir2.depends = pVcfAnno4
		pFiles2Dir2.tag = 'plotstats'
		pFiles2Dir2.input = {pFiles2Dir2.input: lambda ch: [ch.colAt(1).expand(0, '*.stats.csv').toList()*4]}
		
		pStats2Matrix.depends = pFiles2Dir2
		pStats2Matrix.expect  = 'ls {{outdir}}/*.mat.txt'
		pPlotStats.depends = pStats2Matrix
		pyppl().starts(pVcfAnno4).run()
		
		
if __name__ == '__main__':
	unittest.main(failfast=True)