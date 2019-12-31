"""Processes for BED files"""
from diot import Diot
from pyppl import Proc
from . import params, proc_factory
from .bedtools import *

pBedSort = proc_factory(
	desc = 'Sort bed files.',
	config = Diot(annotate = """
	@name:
		pBedSort
	@description:
		Sort bed files
	@input:
		`infile:file`: The input file
	@output:
		`outfile:file`: The output file
	@args:
		`tool`:     The tool used to sort the file. Default: sort (bedtools, bedops)
		`bedtools`: The path to bedtools. Default: bedtools
		`bedops`:   The path to bedops' sort-bed. Default: sort-bed
		`mem`:      The memory to use. Default: 8G
		`by`:       Sort by coordinates("coord", default) or name("name")
			- Only available when use tool `sort`
		`tmpdir`:   The tmpdir to use. Default: `$TMPDIR`
		`unique`:   Remove the dupliated records? Default: True
		`params`:   Other params for `tool`. Default: {}
		`chrorder`: The chromosome order used to sort. Default: `None`
			- `None`: Sort by natural order (chr1 followed by chr10, instead of chr2)
			- Only available when using `sort` (`args.tool = 'sort'`)
	@requires:
		[`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)
		[`bedops`](https://github.com/bedops/bedops)
	"""))
pBedSort.input         = "infile:file"
pBedSort.output        = "outfile:file:{{i.infile | bn}}"
pBedSort.args.tool     = 'sort'
pBedSort.args.bedtools = params.bedtools.value
pBedSort.args.bedops   = params.bedops_sort.value
pBedSort.args.mem      = '8G'
pBedSort.args.by       = 'coord'
pBedSort.args.unique   = True
pBedSort.args.params   = Diot()
pBedSort.args.chrorder = None
pBedSort.args.tmpdir   = params.tmpdir.value
pBedSort.lang          = params.python.value

pBedLiftover = proc_factory(
	desc = 'Lift over bed files.',
	config = Diot(annotate = """
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
	"""))
pBedLiftover.input         = 'infile:file'
pBedLiftover.output        = 'outfile:file:{{i.infile | bn}}, umfile:file:{{i.infile | fn}}.unmapped{{i.infile | ext}}'
pBedLiftover.args.liftover = params.liftover.value
pBedLiftover.args.lochain  = params.lochain.value
pBedLiftover.args.params   = Diot()
pBedLiftover.lang          = params.python.value

pGff2Bed = proc_factory(
	desc   = 'Convert GTF/GFF file to BED file',
	config = Diot(annotate = """
	@input:
		infile: The input gtf/gff file
	@output:
		outfile: The converted bed file
	@args:
		bedcols:  Strings of python functions used to convert GTF/GFF records to BED fields.
			- You can define the NAME column here, and extra columns after the 6th column.
			- For example: `args.bedcols = {"NAME": "lambda attrs: rec.CHR + ':' + rec.START"}`
				- `attrs` are the attributes of GFF records, plus CHR, START, END, SCORE and STRAND.
				- See: https://github.com/pwwang/pygff
				- By default, NAME will use `id` in attributes, and then `name`. Otherwise `CHR:START-END` will  be used.
			- You can also add extra columns starting from 7th column of BED file, for example:
				- `args.bedcols = {"CN": "lambda attrs: attrs['CopyNumber']"}`
		keepattrs: Keep the original attributes at last column of the output BED file.
		outhead: Put head to output file or not.
			- Could be prefix to the head.
	"""))
pGff2Bed.lang   = params.python.value,
pGff2Bed.input  = 'infile:file',
pGff2Bed.output = 'outfile:file:{{i.infile | stem}}.bed',
pGff2Bed.args   = Diot(bedcols = Diot(), keepattrs = True, outhead = '#')

pBedFromGff = pGff2Bed.copy()
