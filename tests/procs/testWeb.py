import unittest, json
from pyppl import PyPPL, Channel
from helpers import getfile, procOK, config
from bioprocs.web import pDownloadForm, pDownloadGet, pDownloadPost

class testWeb (unittest.TestCase):
	
	def testDownloadForm(self):
		url        = "https://www.google.com"
		data1      = {'q': 'pwwang'}
		data2      = {'q': 'job zhou'}
		submit     = "//*[@name='btnG']"
		#submitbtn = '//*[@id="content"]/fieldset[1]/form/table/tbody/tr[3]/td[3]/div/input'
		#nextbtn   = '//*[@id="content"]/span/a[.="Next"]'

		pDownloadForm1 = pDownloadForm.copy()
		pDownloadForm1.input = [(url,data1,submit), (url,data2,submit)]
		pDownloadForm1.forks = 2

		PyPPL(config).start(pDownloadForm1).run()

		pDownloadForm2 = pDownloadForm.copy()
		pDownloadForm2.input = (url, {'q': 'pyppl+pwwang'}, submit, '//*[@id="nav"]//td[last()]//span[.="Next"]')
		PyPPL(config).start(pDownloadForm2).run()

	def testDownloadGet(self):
		pDownloadGet.input = ["http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt"]
		PyPPL(config).start(pDownloadGet).run()
		procOK(pDownloadGet, 'downloadget.txt', self)

	def testDownloadPost(self):
		pDownloadPost.input = ("https://www.w3schools.com/action_page.php", {'firstname': 'Hello', 'lastname': 'Kitty'})
		PyPPL(config).start(pDownloadPost).run()
		procOK(pDownloadPost, 'downloadpost.txt', self)



if __name__ == '__main__':
	unittest.main()
