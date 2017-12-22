import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.genomeplot import pGenomePlot, pGeneTrack, pDataTrack, pAnnoTrack, pUcscTrack, pInteractionTrack

class TestGenomePlot (unittest.TestCase):
	
	def testpGenomePlot (self):
		pGenomePlot0       = pGenomePlot.copy()
		pGenomePlot0.input = [([], "chr10:21234556-23405060")]
		PyPPL(config).start(pGenomePlot0).run()

	def testExtraGeneTrack(self):
		pGeneTrack.input = [('extraGene', "chr10:21234556-23405060")]

		pGenomePlot1         = pGenomePlot.copy()
		pGenomePlot1.depends = pGeneTrack
		pGenomePlot1.input   = lambda ch: [([ch.get()], "chr10:21234556-23405060")]
		PyPPL(config).start(pGeneTrack).run()

	def testDataTrack(self):
		pDataTrack.input = ('bedgraph', getfile('test.bedGraph'), 'chr19')

		pGenomePlot2         = pGenomePlot.copy()
		pGenomePlot2.depends = pDataTrack
		pGenomePlot2.input   = lambda ch: [([ch.get()], "chr19:49231000-49315700")]
		PyPPL(config).start(pDataTrack).run()

	def testAnnoTrack(self):
		pAnnoTrack.input = ('Annotation', getfile('test.bam'), 'chr1')

		pGenomePlot3         = pGenomePlot.copy()
		pGenomePlot3.depends = pAnnoTrack
		pGenomePlot3.input   = lambda ch: [([ch.get()], "chr1:189891483-190087517")]
		PyPPL(config).start(pAnnoTrack).run()

	def testUcscTrack(self):
		pUcscTrack.input                = ('SNPs', 'snp147', 'AnnotationTrack', 'chr1:189891483-190087517')
		pUcscTrack.args.params.start    = 'chromStart'
		pUcscTrack.args.params.end      = 'chromEnd'
		pUcscTrack.args.params.id       = 'name'
		pUcscTrack.args.params.feature  = 'func'
		pUcscTrack.args.params.strand   = 'strand'
		pUcscTrack.args.params.shape    = 'box'
		pUcscTrack.args.params.stacking = 'dense'

		pGenomePlot4         = pGenomePlot.copy()
		pGenomePlot4.depends = pUcscTrack
		pGenomePlot4.input   = lambda ch: [([ch.get()], "chr1:188891483-190087517")]
		PyPPL(config).start(pUcscTrack).run()

	def testInteractionTrack(self):
		pInteractionTrack.input = ('HiC', getfile('chiapet.tool.txt'), "chr10")
		pInteractionTrack.args.intype = "chiapet.tool"
		pInteractionTrack.args.params = {
			'col.interactions'      : "red",
			'col.anchors.fill'      : "blue",
			'col.anchors.line'      : "black",
			'interaction.dimension' : "height",
			'interaction.measure'   : "counts",
			'plot.trans'            : False,
			'plot.outside'          : True,
			'col.outside'           : "lightblue",
			'anchor.height'         : 0.1
		}

		pGenomePlot5         = pGenomePlot.copy()
		pGenomePlot5.depends = pInteractionTrack
		pGenomePlot5.input   = lambda ch: [([ch.get()], "chr10:103449060-103918049")]
		pGenomePlot5.args.params.general.sizes = 'R:c(1, 1, 2, 4)'
		PyPPL(config).start(pInteractionTrack).run()

if __name__ == '__main__':
	unittest.main()
		