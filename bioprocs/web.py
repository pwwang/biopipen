from pyppl import Proc

# Web utils
"""
@name:
	pDownloadPost
@description:
	Download results by submitting a form, supporting pagination.
@input:
	`url`: the URL contains the form
	`submitbtn`: the submit button to click to submit the form (use Xpath).
	`nextbtn`: the button for next page (use Xpath)
	`params`: the params used to fill the form (JSON string or transformed from dict by json.dumps).
@output:
	`outdir:file`: The directory saves the results
@args:
	`interval`: seconds to wait between fetching each page. Default: 1
@requires:
	[`Splinter`](https://splinter.readthedocs.io/en/latest/index.html)
	[`Phantomjs`](http://phantomjs.org/)
"""
pDownloadPost = Proc (desc="Download results by submitting to a form")
pDownloadPost.input  = "url, submitbtn, nextbtn, params"
pDownloadPost.args   = {'interval': 1}
pDownloadPost.output = "outdir:file:outdir-{{#}}"
pDownloadPost.script = """
#!/usr/bin/env python
from splinter import Browser
import json, os, sys, time
if not os.path.exists ("{{out.outdir}}"):
	os.makedirs ("{{out.outdir}}")

browser = Browser('phantomjs')
params = json.loads('''{{params}}''')

sys.stderr.write ("Visiting {{in.url}} ...\\n")
browser.visit("{{in.url}}")

sys.stderr.write ("Filling parameters ...\\n")
for key, val in params.iteritems():
	browser.find_by_xpath(key).fill(val)

i = 0
browser.find_by_xpath('''{{submitbtn}}''').click()
while browser.is_element_present_by_xpath ('''{{nextbtn}}'''):
	sys.stderr.write ("Fetching and saving page %s ...\\n" % (i+1))
	open ("{{out.outdir}}/out-%s.html" % i, 'w').write (browser.html)
	i += 1
	time.sleep({{args.interval}})
	browser.find_by_xpath('''{{nextbtn}}''').click()

open ("{{out.outdir}}/out-%s.html" % i, 'w').write (browser.html)
sys.stderr.write ("Qutting browser ...\\n")
browser.quit()
"""

"""
@name:
	pDownloadGet
@description:
	Download results by urls.
@input:
	`url`: the URLs to download
@output:
	`outfile:file`: The output file
"""
pDownloadGet = Proc (desc="Download from URLs")
pDownloadGet.input  = "url"
pDownloadGet.output = "outfile:file:{{in.url | bn | .replace('?', '__Q__').replace('&', '__N__')  }}"
pDownloadGet.script = """
#!/usr/bin/env python
import urllib
urllib.urlretrieve ('''{{in.url}}''', "{{out.outfile}}")
"""
