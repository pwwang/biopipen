import re
from os import path
from pyppl import Box, utils

def _getsource(*prepends):
	fn  = utils.varname(1, incldot = True)
	cwd = path.dirname(path.realpath(__file__))
	ret = ''
	for prepend in prepends:
		ret += prepend + '\n'
	with open(path.join(cwd, 'scripts', 'utils', fn)) as f:
		ret += f.read()
	return ret

def dirnameFiles(files):
	file0 = sorted(files)[0]
	file0 = path.basename(file0).split('.')[0]
	return re.sub('[^A-Za-z0-9]*\([^)]+?\)|[^A-Za-z0-9]+$', '', file0) + '_etc'

def dirnamePattern(dir, pattern = '*'):
	from glob import glob
	files = glob(path.join(dir, pattern))
	return dirnameFiles(files)

"""
Convert memories
"""
mem2    = Box()
mem2.py = _getsource()
mem2.r  = _getsource()

"""
Run command
"""
runcmd    = Box()
runcmd.py = _getsource()
runcmd.r  = _getsource()

"""
Polling helpers
"""
polling = Box(
	non1st = {},
	first  = {},
	last   = {},
	all    = {}
)
polling.non1st.py = _getsource(runcmd.py)
polling.non1st.r  = _getsource(runcmd.r)
polling.first.py  = _getsource(runcmd.py)
polling.last.r    = _getsource(runcmd.r)
polling.all.py    = _getsource(runcmd.py)
polling.all.r     = _getsource(runcmd.r)

"""
Build reference indices
"""
buildref = Box({
	'fai'  : {},
	'dict' : {},
	'index': {} # general
})
# args.ref: args.samtools needed
buildref.fai.bash  = _getsource()
# args.ref: args.picard needed
buildref.dict.bash = _getsource()
buildref.index.py  = _getsource(polling.non1st.py)

"""
Check reference existence
"""
checkref = Box({
	'fa'  : {},
	'gene': {}
})
# args.ref
checkref.fa.bash   = _getsource()
# args.refgene
checkref.gene.bash = _getsource()

"""
Helpers
"""
helpers = Box({
	'cbindfill'     : {},
	'rbindfill'     : {},
	'params2CmdArgs': {}
})
helpers.cbindfill.r = _getsource()
helpers.rbindfill.r = _getsource()
helpers.params2CmdArgs.py = _getsource()
helpers.params2CmdArgs.r  = _getsource()

plot = Box({
	'hist'   : {},
	'boxplot': {},
	'heatmap': {},
	'maplot' : {},
	'volplot': {}, # volcano plot
	'symbols': {},
	'text'   : {},
	'venn'   : {},
	'upset'  : {},
	'pie'    : {}
})
plot.hist.r    = _getsource()
plot.boxplot.r = _getsource()
plot.maplot.r  = _getsource()
plot.volplot.r = _getsource()
plot.heatmap.r = _getsource()

plot.symbols.r = """
if (!exists('plotSymbols')) {
	plotSymbols = function (x, y, params) {
		do.call(symbols, c(list(x=x, y=y), params))
	}
}
"""

plot.text.r = """
if (!exists('plotText')) {
	plotText = function (x, y, params) {
		do.call(text, c(list(x=x, y=y), params))
	}
}
"""

plot.venn.r  = _getsource()
plot.upset.r = _getsource()
plot.pie.r   = _getsource()

txt = Box({
	'filter'    : {},
	'sampleinfo': {},
	'transform' : {},
})
txt.filter.py    = _getsource()
txt.transform.py = _getsource()

# read sample info file
#
# [Sample	]Patient	Group	Batch
# sample1	patient1	group1	batch1
# sample2	patient1	group2	batch1
# sample3	patient2	group1	batch2
# sample4	patient2	group2	batch2
# Remember NORMAL/CONTROL samples come last if you wanna do comparison between groups
txt.sampleinfo.py = _getsource()
txt.sampleinfo.r  = _getsource()

download = Box({
	'curl': {},
})
download.curl.py = _getsource()

parallel    = Box()
parallel.py = _getsource(runcmd.py)

sql = Box(
	dsnparse    = {},
	schemaparse = {}
)
sql.dsnparse.py    = _getsource()
sql.schemaparse.py = _getsource()

read = Box(
	record = Box(),
	meta   = Box(),
	base   = Box(),
	head   = Box(),
	bed12  = Box(),
	bedpe  = Box(),
	bedx   = Box(),
	bed    = Box(),
)
read.meta.py   = _getsource()
read.record.py = _getsource()
read.base.py   = _getsource(read.meta.py, read.record.py)
read.head.py   = _getsource(read.base.py)
read.bed12.py  = _getsource(read.base.py)
read.bedpe.py  = _getsource(read.base.py)
read.bedx.py   = _getsource(read.base.py)
read.bed.py    = _getsource(read.base.py)

write = Box(
	base = Box(),
	bed  = Box(),
	bedx = Box()
)
write.base.py = _getsource(read.meta.py)
write.bed.py  = _getsource(read.bed.py, write.base.py)
write.bedx.py = _getsource(read.bedx.py, write.base.py)


genenorm    = Box()
genenorm.py = _getsource(read.head.py, read.bed12.py, read.bedpe.py, read.bedx.py, read.bed.py, write.base.py)
