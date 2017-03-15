from pyppl import proc

"""
@name:
	pSort
@description:
	Sort file using linux command `sort`
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file
@args:
	`sort-args`: The arguments used by `sort`
"""
pSort = proc ()
pSort.input  = "infile:file"
pSort.output = "outfile:file:{{infile.bn}}.sorted"
pSort.args   = {"sort-args": ""}
pSort.script = """
#!/usr/bin/env bash
sort {{proc.args.sort-args}} "{{infile}}" > "{{outfile}}"
"""