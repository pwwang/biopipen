import unittest
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.hic import pPartners

class TestHiC (unittest.TestCase):
	
	def testpPartners (self):
		pPartners.input = (getfile('regions.bed'), getfile('interactions.bedpe'))
		
		PyPPL(config).start(pPartners).run()
		

		
if __name__ == '__main__':
	unittest.main()