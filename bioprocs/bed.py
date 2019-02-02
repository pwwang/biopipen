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
	`by`:           Sort by coordinates("coord", default) or name("name")
		- Only available when use tool `sort`
	`tmpdir`:       The tmpdir to use. Default: `$TMPDIR`
	`unique`:       Remove the dupliated records? Default: True
	`params`:       Other params for `tool`. Default: {}
@requires:
	[`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)
	[`bedops`](https://github.com/bedops/bedops)
"""
pBedSort               = Proc(desc = 'Sort bed files.')
pBedSort.input         = "infile:file"
pBedSort.output        = "outfile:file:{{i.infile | bn}}"
pBedSort.args.tool     = 'sort'
pBedSort.args.bedtools = params.bedtools.value
pBedSort.args.bedops   = params.bedops_sort.value
pBedSort.args.mem      = '8G'
pBedSort.args.by       = 'coord'
pBedSort.args.unique   = True
pBedSort.args.params   = Box()
pBedSort.args.tmpdir   = params.tmpdir.value
pBedSort.lang          = params.python.value
pBedSort.script        = "file:scripts/bed/pBedSort.py"

"""
@name:
	pBedLiftover
@description:
	Lift over bed files.
@input:
	`infile:file`: The input bed file
@output:
	`outfile:file`: The output file
	`umfile:file` : The unmapped file
@args:
	`liftover`: The liftover program
	`lochain` : the liftover chain file
@require:
	`liftover` from UCSC
"""
pBedLiftover               = Proc(desc = 'Lift over bed files.')
pBedLiftover.input         = 'infile:file'
pBedLiftover.output        = 'outfile:file:{{i.infile | bn}}, umfile:file:{{i.infile | fn}}.unmapped{{i.infile | ext}}'
pBedLiftover.args.liftover = params.liftover.value
pBedLiftover.args.lochain  = params.lochai.value
pBedLiftover.args.params   = Box()
pBedLiftover.lang          = params.python.value
pBedLiftover.script        = "file:scripts/bed/pBedLiftover.py"

"""
@name:
	pGff2Bed
@description:
	Convert GTF/GFF file to BED file
@input:
	`infile:file`: The input gtf/gff file
@output:
	`outfile:file`: The converted bed file
@args:
	`attr2name`: The function used to convert attributes from GTF/GFF file to BED field 'name'
"""
pGff2Bed                = Proc(desc = 'Convert GTF/GFF file to BED file')
pGff2Bed.input          = 'infile:file'
pGff2Bed.output         = 'outfile:file:{{i.infile | fn}}.bed'
pGff2Bed.args.attr2name = None
pGff2Bed.args.keepinfo  = True
pGff2Bed.lang           = params.python.value
pGff2Bed.script         = "file:scripts/bed/pGff2Bed.py"