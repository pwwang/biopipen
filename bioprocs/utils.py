from os import path
from pyppl import Box, utils

def _getsource(prepend = ''):
	fn  = utils.varname(1, incldot = True)
	cwd = path.dirname(path.realpath(__file__))
	ret = prepend
	with open(path.join(cwd, 'scripts', 'utils', fn)) as f:
		ret += ''.join(f.readlines()) + '\n'
	return ret

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
polling = Box({
	'non1st': {},
	'first' : {},
	'all'   : {}
})
polling.non1st.py = _getsource(runcmd.py)
polling.non1st.r  = _getsource(runcmd.r)
polling.first.py  = _getsource(runcmd.py)
polling.first.r   = _getsource(runcmd.r)
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
	'params2CmdArgs': {}
})
helpers.cbindfill.r = _getsource()
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

plot.venn.r = _getsource()

plot.upset.r = """
if (!exists('plotUpset')) {
	plotUpset = function(mat, filename, params) {
		library(UpSetR)
		png (filename, res=300, width=2000, height=2000)
		if (! "nintersects" %in% names(params)) {
			params$nintersects = NA
		}
		v = do.call (upset, c(list(data = mat, sets = colnames(mat)), params))
		print (v)
		dev.off()
	}
}
"""

txt = Box({
	'filter'    : {},
	'sampleinfo': {},
	'transform' : {},
})
txt.filter.py    = _getsource()
txt.transform.py = _getsource()

# read sample info file
#
# Sample	Patient	Group	Batch
# sample1	patient1	group1	batch1
# sample2	patient1	group2	batch1
# sample3	patient2	group1	batch2
# sample4	patient2	group2	batch2
# Remember NORMAL/CONTROL samples come last if you wanna do comparison between groups
txt.sampleinfo.py = _getsource()
txt.sampleinfo.r      = _getsource()

download = Box({
	'curl': {},
})
download.curl.py = _getsource()

