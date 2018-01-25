from pyppl import Proc, Box
from .bedtools import pBedGetfasta, pBedRandom, pBedFlank, pBedIntersect, pBedIntersect2, pBedMerge, pBedsMerge
from . import params
from .utils import helpers, runcmd

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
pBedSort.envs.runcmd         = runcmd.py
pBedSort.envs.params2CmdArgs = helpers.params2CmdArgs.py
pBedSort.script              = "file:scripts/bed/pBedSort.py"

"""
@name:
	pBedCluster
@description:
	Assign cluster id to each record
@input:
	`infile:file`: The input bed file
@output:
	`outfile:file`: The output file
@args:
	`tool`:         The tool used to sort the file. Default: bedtools
	`bedtools`:     The path to bedtools. Default: bedtools
	`params`:       Other params for `tool`. Default: ''
@requires:
	[`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)
"""
pBedCluster               = Proc(desc = 'Assign cluster id to each record.')
pBedCluster.input         = "infile:file"
pBedCluster.output        = "outfile:file:{{infile | fn}}.cluster.bed"
pBedCluster.args.tool     = 'bedtools'
pBedCluster.args.bedtools = 'bedtools'
pBedCluster.args.params   = ''
pBedCluster.script        = """
case {{args.tool | quote}} in 
	bedtools)
		{{args.bedtools}} cluster -i "{{infile}}" {{args.params}} > "{{outfile}}"
esac
"""