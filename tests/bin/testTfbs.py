import unittest
from helpers import getfile, getbin

class testTFBS (unittest.TestCase):
	def testPWMScanNoArgs(self):
		pwmscan = getbin('pwmscan')
		


if __name__ == '__main__':
	unittest.main(failfast = True)