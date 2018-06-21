import testly, helpers
from os import path
from pyppl import Channel, PyPPL, Proc
from bioaggrs.wxs import aBam2MutCnv
from bioprocs.sambam import pBamMarkdup, pSam2Bam
from bioprocs.fastx import pFastq2Sam
from bioprocs.common import pFiles2Dir

class TestWxs(testly.TestCase):
	
	# use the files from tests/procs/TestSambam
	testdir, indir, outdir      = helpers.testdirs('TestSambam')
	
	def dataProvider_test_aBam2MutCnv(self):
		pFastq2Sam.args.nthread = 2
		pFastq2Sam.forks        = 2
		pFastq2Sam.input        = Channel.fromPairs(path.join(self.indir, '*.fastq'))
		pFastq2Sam.args.ref     = path.join(self.indir, 'ref.fa')
		pSam2Bam.depends        = pFastq2Sam
		pSam2Bam.forks          = 2
		pSam2Bam.args.markdup   = True
		pFiles2Dir.depends      = pSam2Bam
		pFiles2Dir.input        = lambda ch: [ch.flatten()]
		pSampleInfo             = Proc(desc = 'Generate sample information file from directory', tag = 'test')
		pSampleInfo.depends     = pFiles2Dir
		pSampleInfo.input       = 'indir:file'
		pSampleInfo.output      = 'outfile:file:{{in.indir | fn}}.saminfo'
		pSampleInfo.lang        = 'python'
		pSampleInfo.script      = """#
		# PYPPL INDENT REMOVE 
		from os import path
		from glob import glob
		from bioprocs.utils.tsvio import TsvWriter, TsvRecord
		bams                    = glob(path.join({{in.indir | quote}}, '*.bam'))
		writer                  = TsvWriter({{out.outfile | quote}})
		writer.meta.add('Sample', 'Patient', 'Group')
		writer.writeHead()
		for i, bam in enumerate(bams):
			r                   = TsvRecord()
			r.Sample            = path.basename(bam)
			r.Patient           = 'patient{index}'.format(index = int(i/2))
			r.Group             = ['Normal', 'Tumor'][i%2]
			writer.write(r)
		"""
		PyPPL(helpers.config).start(pFastq2Sam).run()
		yield pFiles2Dir.channel.get(), pSampleInfo.channel.get(), '', '', ''
		
	
	def test_aBam2MutCnv(self, samdir, saminfo, gmutfile, smutfile, cnvfile):
		aBam2MutCnvTest          = aBam2MutCnv.copy()
		aBam2MutCnvTest.input    = [[samdir], [saminfo]]
		aBam2MutCnvTest.forks    = 2
		aBam2MutCnvTest.args.ref = path.join(self.indir, 'ref.fa')
		PyPPL(helpers.config).start(aBam2MutCnvTest).run()
		
	
if __name__ == '__main__':
	testly.main()