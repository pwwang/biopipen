import helpers, unittest, testly
from os import path, makedirs
from bioprocs.utils.gene import genenorm
from medoo import Medoo

class TestUtilsGene(helpers.TestCase):

	testdir, indir, outdir = helpers.testdirs('TestUtilsGene')

	def dataProvider_testGenenorm(self):
		infile = path.join(self.indir, 'genenorm.txt')
		yield infile, None, 'skip', 'symbol, alias', 'symbol', 'hg19', {'cnames': ['Gene', 'Index']}, {}, None, self.testdir, {
			'CDH13': {'alias': ['CDHH', 'P105'], 'symbol': 'CDH13', '_id': '1012', 'taxid': 9606},
			'CDKN2': {'alias': ['CDKN2', 'p33(CDK2)'], 'symbol': 'CDK2', '_id': '1017', 'taxid': 9606},
			'TAF1A': {'alias': ['MGC:17061', 'RAFI48', 'SL1', 'TAFI48'], 'symbol': 'TAF1A', '_id': '9015', 'taxid': 9606},
			'CRYG5': {'alias': ['CRY-g-A', 'CRYG1', 'CRYG5'], 'symbol': 'CRYGA', '_id': '1418', 'taxid': 9606},
			'nfs1': {'alias': ['HUSSY-08', 'IscS', 'NIFS'], 'symbol': 'NFS1', '_id': '9054', 'taxid': 9606},
			'CALM1P1': {'alias': ['CALML1'], 'symbol': 'CALM1P1', '_id': '802', 'taxid': 9606},
			'RCBTB2': {'alias': ['CHC1L', 'RLG'], 'symbol': 'RCBTB2', '_id': '1102', 'taxid': 9606},
			'SLC25A5P3': {'alias': [], 'symbol': 'SLC25A5P3', '_id': 'ENSG00000213673', 'taxid': 9606},
			'znf439': {'alias': [], 'symbol': 'ZNF439', '_id': '90594', 'taxid': 9606},
			'GLE1': {'alias': ['GLE1L', 'LCCS', 'LCCS1', 'hGLE1'], 'symbol': 'GLE1', '_id': '2733', 'taxid': 9606}
		}

	def testGenenorm(self, infile, outfile, notfound, frm, to, genome, inopts, outopts, genecol, cachedir, output):
		self.maxDiff = None

		ret = genenorm(infile, outfile, notfound, frm, to, genome, inopts, outopts, genecol, cachedir)
		for key, val in output.items():
			self.assertDictEqual(ret[key], val)

	def dataProvider_testGenenormCached(self):
		infile = path.join(self.indir, 'genenorm.txt')
		infile2 = path.join(self.indir, 'genenorm2.txt')
		cachedir = path.join(self.testdir, 'cached')
		makedirs(cachedir)
		# cache part of the query
		genenorm(infile2, None, 'skip', 'symbol, alias', 'alias', 'hg19', {'cnames': ['Gene', 'Index']}, {}, None, cachedir)

		yield infile, None, 'skip', 'symbol, alias', 'symbol', 'hg19', {'cnames': ['Gene', 'Index']}, {}, None, cachedir, {
			'CDH13': {'alias': ['CDHH', 'P105'], 'symbol': u'CDH13', '_id': '1012', 'taxid': 9606},
			'CDKN2': {'alias': ['CDKN2', 'p33(CDK2)'], 'symbol': u'CDK2', '_id': '1017', 'taxid': 9606},
			'TAF1A': {'alias': ['MGC:17061', 'RAFI48', 'SL1', 'TAFI48'], 'symbol': u'TAF1A', '_id': '9015', 'taxid': 9606},
			'CRYG5': {'alias': ['CRY-g-A', 'CRYG1', 'CRYG5'], 'symbol': u'CRYGA', '_id': '1418', 'taxid': 9606},
			'nfs1': {'alias': ['HUSSY-08', 'IscS', 'NIFS'], 'symbol': u'NFS1', '_id': '9054', 'taxid': 9606},
			'CALM1P1': {'alias': ['CALML1'], 'symbol': u'CALM1P1', '_id': u'802', 'taxid': 9606},
			'RCBTB2': {'alias': ['CHC1L', 'RLG'], 'symbol': u'RCBTB2', '_id': '1102', 'taxid': 9606},
			'SLC25A5P3': {'alias': ['ANTP1'], 'symbol': u'SLC25A5P3', '_id': u'222005', 'taxid': 9606},
			'znf439': {'alias': [], 'symbol': u'ZNF439', '_id': '90594', 'taxid': 9606},
			'GLE1': {'alias': ['GLE1L', 'LCCS', 'LCCS1', 'hGLE1'], 'symbol': u'GLE1', '_id': '2733', 'taxid': 9606}
		}

	def testGenenormCached(self, infile, outfile, notfound, frm, to, genome, inopts, outopts, genecol, cachedir, output):
		self.maxDiff = None
		ret = genenorm(infile, outfile, notfound, frm, to, genome, inopts, outopts, genecol, cachedir)
		for key, val in output.items():
			retalias = ret[key]['alias']
			outalias = val['alias']
			del ret[key]['alias']
			del val['alias']
			self.assertCountEqual(retalias, outalias)
			self.assertDictEqual(dict(ret[key]), val)

	def dataProvider_testGenenormOutfile(self):
		infile = path.join(self.indir, 'genenorm.txt')
		cachedir = path.join(self.testdir, 'outfile')
		makedirs(cachedir)
		# cache part of the query
		genenorm(infile, None, 'skip', 'symbol, alias', 'alias', 'hg19', {'cnames': ['Gene', 'Index']}, {}, None, cachedir)

		outfile = path.join(self.testdir, 'genenorm.out')
		exptfile = path.join(self.outdir, 'genenorm.out')

		yield infile, outfile, 'ignore', 'symbol, alias', 'symbol', 'hg19', {'cnames': ['Gene', 'Index']}, {'query': True}, None, cachedir, exptfile

	def testGenenormOutfile(self, infile, outfile, notfound, frm, to, genome, inopts, outopts, genecol, cachedir, exptfile):
		self.maxDiff = None
		ret = genenorm(infile, outfile, notfound, frm, to, genome, inopts, outopts, genecol, cachedir)
		self.assertFileEqual(outfile, exptfile)

if __name__ == '__main__':
	testly.main(verbosity = 2)
