from bioprocs.utils import runcmd, cmdargs

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
header  = {{args.header | repr}}
by      = {{args.by | quote}}
tool    = {{args.tool | quote}}

if header:
	# write the header to outfile
	params = cmdargs({
			'e': '^#',
		},
		dash  = '-',
		equal = ' '
	)
	if infile.endswith('.gz'):
		cmd = 'zcat {infile} | grep {cmdargs} > {outfile}'
	else:
		cmd = 'grep {cmdargs} {infile} > {outfile}'
	runcmd(cmd.format(
		cmdargs = params,
		infile  = infile,
		outfile = outfile
	))

if tool == 'sort':
	if infile.endswith('.gz'):
		cmd = 'zcat {infile} | grep "^#" | sort {cmdargs} >> {outfile}'
	else:
		cmd = 'grep -v "^#" {infile} | sort {cmdargs} >> {outfile}'

	if by.lower().startswith('coord'):
		params = cmdargs({
			'k#1': '1,1',
			'k#2': '2,2n'
		}, dash = '-', equal = ' ')
	else:
		params = cmdargs({
			'k#1': '3,3',
			'k#2': '1,1',
			'k#3': '2,2n'
		}, dash = '-', equal = ' ')
	runcmd(cmd.format(
		cmdargs = params,
		infile  = infile,
		outfile = outfile
	))
