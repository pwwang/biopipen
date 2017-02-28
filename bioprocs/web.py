from pyppl import proc

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
	- [`Splinter`](https://splinter.readthedocs.io/en/latest/index.html)
	- [`Phantomjs`](http://phantomjs.org/)
"""
pDownloadPost = proc ()
pDownloadPost.input  = "url, submitbtn, nextbtn, params"
pDownloadPost.args   = {'interval': 1}
pDownloadPost.output = "outdir:file:outdir-{{#}}"
pDownloadPost.script = """
#!/usr/bin/env python
from splinter import Browser
import json, os, sys, time
if not os.path.exists ("{{outdir}}"):
	os.makedirs ("{{outdir}}")

browser = Browser('phantomjs')
params = json.loads('''{{params}}''')

sys.stderr.write ("Visiting {{url}} ...\\n")
browser.visit("{{url}}")

sys.stderr.write ("Filling parameters ...\\n")
for key, val in params.iteritems():
	browser.find_by_xpath(key).fill(val)

i = 0
browser.find_by_xpath('''{{submitbtn}}''').click()
while browser.is_element_present_by_xpath ('''{{nextbtn}}'''):
	sys.stderr.write ("Fetching and saving page %s ...\\n" % (i+1))
	open ("{{outdir}}/out-%s.html" % i, 'w').write (browser.html)
	i += 1
	time.sleep({{proc.args.interval}})
	browser.find_by_xpath('''{{nextbtn}}''').click()

open ("{{outdir}}/out-%s.html" % i, 'w').write (browser.html)
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
@args:
	`keepname`: bool, whether to keep the basename, otherwise use {{#}}.<ext>, default: True
@output:
	`outdir:file`: The directory saves the results
"""
pDownloadGet = proc ()
pDownloadGet.input  = "url"
pDownloadGet.args   = {'keepname': True}
pDownloadGet.output = "outdir:file:outdir"
pDownloadGet.script = """
#!/usr/bin/env python
import time, urllib, os
from urlparse import urlparse
url = urlparse('''{{url}}''')
ext = os.path.splitext (url.path)[1]
if {{proc.args.keepname}}:
	name = os.path.basename (url.path)
else:
	name = "{{#}}%s" % ext
if not os.path.exists("{{outdir}}"):
	os.makedirs("{{outdir}}")
outfile = os.path.join ("{{outdir}}", name)

urllib.urlretrieve (url.geturl(), outfile)
"""
