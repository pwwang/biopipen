from pyppl import proc
from .utils import txt, download

"""
@name:
	pTxt
@description:
	Download CSV format files.
@input:
	`in`: The name of the resource
@output:
	`outfile:file`: The output file
@args:
	`cols`:      Select the columns to keep. Default: '' (all cols)
	`rowfilter`: Filter rows. For example, to filter out rows not start with 'Chr':
		- `"lambda x: not x[0].startswith('Chr')"`
		- Note that rowfilter applied before cols filter.
	`urls`:      Available resources and their urls.
	`gz`:        Whether to gzip the output file.
@requires:
	[`curl`](https://en.wikipedia.org/wiki/CURL)
"""
pTxt                    = proc(desc = 'Download CSV format files.')
pTxt.input              = "in"
pTxt.output             = "outfile:file:{{in}}.txt{{args.gz | lambda x: '.gz' if x else ''}}"
pTxt.args.gz            = False
pTxt.args.delimit       = "\t"
pTxt.args.cols          = ''
pTxt.args.header        = True
pTxt.args.rowfilter     = ''
pTxt.args.username      = ''
pTxt.args.password      = ''
pTxt.args.curl          = 'curl'
pTxt.args._txtFilter    = txt.filter.python
pTxt.args._downloadCurl = download.curl.python
pTxt.args.urls          = {
	'drugbank-target-all': 'https://www.drugbank.ca/releases/5-0-7/downloads/target-all-uniprot-links',
	'drugbank-target-approved': 'https://www.drugbank.ca/releases/5-0-7/downloads/target-approved-uniprot-links'
}
pTxt.lang   = 'python'
pTxt.script = """
import os, shutil
from subprocess import check_output

{{args._downloadCurl}}

url = {{args.urls | json}}["{{in}}"]
tmpdir = "{{job.outdir}}/tmp"
if not os.path.exists(tmpdir):
	os.makedirs(tmpdir)
	
downfile = os.path.join(tmpdir, 'downloaded')
downloadCurl(url, downfile, {{args.username | quote}}, {{args.password | quote}}, {{args.curl | quote}})

# determin file types
output = check_output(['file', downfile])
if 'gzip' in output:
	ugfile = downfile + '.ungz'
	with open(ugfile, 'w') as f:
		f.write(check_output(['gunzip', downfile, '-c']))
	downfile = ugfile
elif 'Zip' in output:
	zipdir = os.path.join(tmpdir, '_unzipped')
	import zipfile, glob
	zipref = zipfile.ZipFile(downfile, 'r')
	zipref.extractall(zipdir)
	zipref.close()
	downfile = glob.glob(os.path.join(zipdir, '*'))[0]

{{args._txtFilter}}
txtFilter(downfile, "{{outfile}}", {{args.cols | lambda x: __import__('json').dumps(x) if isinstance(x, list) else '"'+ x +'"'}}, {{args.rowfilter | lambda x: 'False' if not x else x}}, {{args.header}}, {{args.delimit | quote}})

#shutil.rmtree(tmpdir)
"""