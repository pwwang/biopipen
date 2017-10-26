import unittest, addPath
from os import path
from pyppl import PyPPL, Channel, Proc
from bioprocs import params
from bioprocs.snp import pSnp2Bed, pSnp2Avinput

params = params.toDict()

snps = [
	'rs963932939', 'rs894694327', 'rs1014515781', 'rs558262067', 'rs970434764', 'rs988360965', 'rs913337047', 'rs547820522', 
	'rs573143395', 'rs6693167', 'rs905822467', 'rs530715173', 'rs957863901', 'rs552799991', 'rs989309163', 'rs574400275', 
	'rs1038648395', 'rs541590342', 'rs777070779', 'rs192054622', 'rs552022974', 'rs931861179', 'rs987418923', 'rs879195221', 
	'rs533021860', 'rs1007893826', 'rs940724512', 'rs1017978027', 'rs524312', 'rs770328307'
]

snpfile = path.join(params.tmpdir, 'snps.txt')
snpstr  = ''
if path.isfile(snpfile):
	with open(snpfile) as f: snpstr = f.read()

osnpstr  = 'to be skipped\n'
osnpstr += '# skip again\n'
for i, snp in enumerate(snps):
	osnpstr += '%s\t%s\n' % (i, snp)

if osnpstr != snpstr:
	with open(snpfile, 'w') as f:
		f.write(osnpstr)

class TestGene (unittest.TestCase):
	def testpSnp2Bed(self):
		pSnp2Bed.input     = [snpfile]
		pSnp2Bed.errhow    = 'terminate'
		pSnp2Bed.args.skip = 1
		pSnp2Bed.args.col  = 1
		PyPPL().start(pSnp2Bed).run()

	def testpSnp2Avinput(self):
		pSnp2Avinput.input = [snpfile]
		pSnp2Avinput.errhow    = 'terminate'
		pSnp2Avinput.args.skip = 1
		pSnp2Avinput.args.col  = 1
		PyPPL().start(pSnp2Avinput).run()

if __name__ == '__main__':
	unittest.main(failfast=True)