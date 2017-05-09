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
	`skip`:   To skip first N lines. Default: 0
	`params`: The arguments used by `sort`
"""
pSort = proc ()
pSort.input  = "infile:file"
pSort.output = "outfile:file:{{infile.bn}}.sorted"
pSort.args   = {"params": "", "skip": 0}
pSort.script = """
#!/usr/bin/env bash
if [[ "{{proc.args.skip}}" == "0" ]]; then
	sort {{proc.args.params}} "{{infile}}" > "{{outfile}}"
else
	nhead={{proc.args.skip}}
	ntail=$((nhead+1))
	( head -n $nhead "{{infile}}" && tail -n +$ntail "{{infile}}" | sort {{proc.args.params}} ) > "{{outfile}}"
fi
"""