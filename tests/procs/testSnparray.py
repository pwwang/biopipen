import unittest
from os import path
from subprocess import check_output, CalledProcessError
from pyppl import PyPPL, Proc
from helpers import getfile, procOK, config, fileOKIn
from bioprocs.snparray import pGistic
from bioprocs import params

def gisticdir():
	try:
		gisticbin = check_output(['which', params.gistic.value]).strip()
		return path.dirname(path.realpath(gisticbin))
	except CalledProcessError:
		return ''

class TestSnparray (unittest.TestCase):
	
	@unittest.skipIf(not path.isdir(gisticdir()), 'Gistic not installed, or gistic2 not in PATH.')
	def test2_pGistic (self):
		gdir             = gisticdir()
		segfile          = path.join(gdir, 'examplefiles', 'segmentationfile.txt')
		mkfile           = path.join(gdir, 'examplefiles', 'markersfile.txt')
		alfile           = path.join(gdir, 'examplefiles', 'arraylistfile.txt')
		cnvfile          = path.join(gdir, 'examplefiles', 'cnvfile.txt')
		pGistic.input    = (segfile, mkfile, alfile, cnvfile)
		pGistic.args.mcr = '/home/biotools/matlab/r2014a/'
		PyPPL(config).start(pGistic).run()

		
if __name__ == '__main__':
	unittest.main(failfast=True)