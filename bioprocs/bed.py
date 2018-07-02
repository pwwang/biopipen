from pyppl import Proc, Box
from .bedtools import *
from . import params

"""
@name:
	pBedSort
@description:
	Sort bed files
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file
@args:
	`tool`:         The tool used to sort the file. Default: sort (bedtools, bedops)
	`bedtools`:     The path to bedtools. Default: bedtools
	`bedops_sort`:  The path to bedops' sort-bed. Default: sort-bed
	`mem`:          The memory to use. Default: 8G
	`tmpdir`:       The tmpdir to use. Default: `$TMPDIR`
	`unique`:       Remove the dupliated records? Default: True
	`params`:       Other params for `tool`. Default: {}
@requires:
	[`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)
	[`bedops`](https://github.com/bedops/bedops)
"""
pBedSort                     = Proc(desc = 'Sort bed files.')
pBedSort.input               = "infile:file"
pBedSort.output              = "outfile:file:{{in.infile | bn}}"
pBedSort.args.tool           = 'sort'
pBedSort.args.bedtools       = params.bedtools.value
pBedSort.args.bedops         = params.bedops_sort.value
pBedSort.args.mem            = '8G'
pBedSort.args.unique         = True
pBedSort.args.params         = Box()
pBedSort.args.tmpdir         = params.tmpdir.value
pBedSort.lang                = params.python.value
pBedSort.script              = "file:scripts/bed/pBedSort.py"

pBedLiftover               = Proc(desc = 'Lift over bed files.')
pBedLiftover.input         = 'infile:file'
pBedLiftover.output        = 'outfile:file:{{in.infile | bn}}, umfile:file:{{in.infile | fn}}.unmapped{{in.infile | ext}}'
pBedLiftover.args.liftover = params.liftover.value
pBedLiftover.args.lochain  = params.lochain.value
pBedLiftover.args.params   = Box()
pBedLiftover.lang          = params.python.value
pBedLiftover.script        = "file:scripts/bed/pBedLiftover.py"