"""Download/Get data from Websites instead of APIs"""
from pyppl import Proc
from . import params
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

# Web utils
@procfactory
def _pDownloadForm():
	"""
	@name:
		pDownloadForm
	@description:
		Download results by submitting a form, supporting pagination.
	@input:
		`url`   : the URL contains the form
		`data`  : the data used to fill the form (JSON string or transformed from dict by json.dumps).
		`submit`: the submit button to submit the form (use Xpath).
		`next`  : the button for next page (use Xpath)
	@output:
		`outdir:file`: The directory saves the results
	@args:
		`interval`: seconds to wait between fetching each page. Default: 1
	@requires:
		[`Splinter`](https://splinter.readthedocs.io/en/latest/index.html)
		[`Phantomjs`](http://phantomjs.org/)
	"""
	pDownloadForm               = Proc (desc="Download results by submitting to a form")
	pDownloadForm.input         = "url, data, submit, next"
	pDownloadForm.output        = "outdir:dir:{{i.url | bn | lambda x: x if x else 'outdir'}}"
	pDownloadForm.args.interval = 1
	pDownloadForm.lang          = params.python.value
	pDownloadForm.script        = "file:scripts/web/pDownloadForm.py"
	return pDownloadForm

@procfactory
def _pDownloadGet(alias = 'pDownload'):
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
	pDownloadGet.output = "outfile:file:{{i.url | bn | lambda x: x if x else 'outfile' | .replace('?', '__Q__').replace('&', '__N__')}}"
	pDownloadGet.lang   = params.python.value
	pDownloadGet.script = "file:scripts/web/pDownloadGet.py"
	return pDownloadGet

@procfactory
def _pDownloadPost():
	"""
	@name:
		pDownloadPost
	@description:
		Download results by POST.
	@input:
		`url` : the URLs to download
		`data`: the POST data.
	@output:
		`outfile:file`: The output file
	"""
	pDownloadPost = Proc (desc="Download from URLs")
	pDownloadPost.input  = "url, data"
	pDownloadPost.output = "outfile:file:{{i.url | bn | lambda x: x if x else 'outfile'}}"
	pDownloadPost.lang   = params.python.value
	pDownloadPost.script = "file:scripts/web/pDownloadPost.py"
	return pDownloadPost

