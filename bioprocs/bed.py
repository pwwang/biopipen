from pyppl import proc
"""
Utilities for bed files.
"""

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
	`sort`:         The path to linux's sort. Default: sort
	`mem`:          The memory to use. Default: 8G
	`tmpdir`:       The tmpdir to use. Default: `$TMPDIR`
	`unique`:       Remove the dupliated records? Default: True
	`params`:       Other params for `tool`. Default: ''
@requires:
	[`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)
	[`bedops`](https://github.com/bedops/bedops)
"""
pBedSort                   = proc(desc = 'Sort bed files.')
pBedSort.input             = "infile:file"
pBedSort.output            = "outfile:file:{{infile | bn}}"
pBedSort.args.tool         = 'sort'
pBedSort.args.bedtools     = 'bedtools'
pBedSort.args.bedops_sort  = 'bedops'
pBedSort.args.sort         = 'sort'
pBedSort.args.mem          = '8G'
pBedSort.args.unique       = True
pBedSort.args.params       = ''
pBedSort.args.tmpdir       = __import__('tempfile').gettempdir()

pBedSort.script            = """
export TMPDIR="{{args.tmpdir}}"
uq=""

case {{args.tool | quote}} in 
	sort)
		cmd="{{args.sort}} -k1,1 -k2,2n {{args.params}} \"{{infile}}\""
		;;
	bedtools)
		cmd="{{args.bedtools}} sort -i \"{{infile}}\" {{args.params}}"
		;;
	bedops)
		cmd="{{args.bedops_sort}} --max-mem {{args.mem}} --tmpdir \"{{args.tmpdir}}\" {{args.params}} \"{{infile}}\""
		;;
esac
if [[ "{{args.unique | Rbool}}" == "TRUE" ]]; then
	exec $cmd | uniq > "{{outfile}}"
else
	exec $cmd > "{{outfile}}"
fi
"""

"""
@name:
	pBedIntersect
@description:
	Find intersections of two bed files.
	Input files must be sorted.
@input:
	`infile1:file`: The 1st input bed file
	`infile2:file`: The 2nd input bed file
@output:
	`outfile:file`: The output file
@args:
	`tool`:         The tool used to sort the file. Default: bedtools (bedops)
	`bedtools`:     The path to bedtools. Default: bedtools
	`bedops`:  The path to bedops. Default: bedops
	`params`:       Other params for `tool`. Default: ''
@requires:
	[`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)
	[`bedops`](https://github.com/bedops/bedops)
"""
pBedIntersect               = proc(desc = 'Find intersections of two bed files.')
pBedIntersect.input         = "infile1:file, infile2:file"
pBedIntersect.output        = "outfile:file:{{infile1 | fn | fn}}-{{infile2 | fn | fn}}.intersect.bed"
pBedIntersect.args.tool     = 'bedtools'
pBedIntersect.args.bedtools = 'bedtools'
pBedIntersect.args.bedops   = 'bedops'
pBedIntersect.args.params   = ''
pBedIntersect.script        = """
case {{args.tool | quote}} in 
	bedtools)
		{{args.bedtools}} intersect -a "{{infile1}}" -b "{{infile2}}" {{args.params}} > "{{outfile}}"
		;;
	bedops)
		{{args.bedops}} --intersect "{{infile1}}" "{{infile2}}" {{args.params}} > "{{outfile}}"
		;;
esac
"""

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
pBedCluster               = proc(desc = 'Assign cluster id to each record.')
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