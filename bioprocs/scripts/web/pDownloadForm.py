
from splinter import Browser
import json, os, sys, time


browser = Browser('phantomjs')

try:
	data    = {{in.data}}

	sys.stderr.write ('Visiting {{in.url | quote}} ...\n')
	browser.visit({{in.url | quote}})

	sys.stderr.write ("Filling parameters ...\n")
	for key, val in data.iteritems():
		print key
		if key.startswith('/'):
			element = browser.find_by_xpath(key)
			if not element:
				raise KeyError('Cannot find element: %s' % key)
			element.fill(val)
		else:
			browser.fill(key, val)

	i = 0
	browser.find_by_xpath({{in.submit | quote}}).first.click()

	while {{in.next | lambda x: bool(x)}} and \
		browser.is_element_present_by_xpath ({{in.next | quote}}):
		sys.stderr.write ("Fetching and saving page %s ...\n" % (i+1))
		with open ("{{out.outdir}}/out-%s.html" % i, 'w') as f:
			f.write (browser.html.encode('utf-8'))
		i += 1
		time.sleep({{args.interval}})
		browser.find_by_xpath({{in.next | quote}}).first.click()

	sys.stderr.write ("Fetching and saving page %s ...\n" % (i+1))
	with open ("{{out.outdir}}/out-%s.html" % i, 'w') as f:
		f.write (browser.html.encode('utf-8'))
	sys.stderr.write ("Qutting browser ...\n")
except Exception:
	raise
finally:
	browser.quit()
