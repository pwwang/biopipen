from bioprocs.utils import runcmd, cmdargs

infile1  = {{in.infile1 | quote}}
infile2  = {{in.infile2 | quote}}
outfile  = {{out.outfile | quote}}
rmany    = {{args.any | repr}}
tool     = {{args.tool | quote}}
header   = {{args.header | repr}}
bedtools = {{args.bedtools | quote}}

if tool == 'bedtools':
	if header:
		# write the header to outfile
		params = cmdargs({
				'e': '^#',
			},
			dash  = '-',
			equal = ' '
		)
		if infile1.endswith('.gz'):
			cmd = 'zcat {infile} | grep {cmdargs} > {outfile}'
		else:
			cmd = 'grep {cmdargs} {infile} > {outfile}'
		runcmd(cmd.format(
			cmdargs = params,
			infile  = infile1,
			outfile = outfile
		))

	cmd = '{bedtools} subtract {params} >> {outfile}'
	params = {
		'a': infile1,
		'b': infile2
	}
	if rmany:
		params['A'] = True
	
	params = cmdargs(params, dash = '-', equal = ' ')
	runcmd(cmd.format(
		bedtools = bedtools,
		params   = params,
		outfile  = outfile
	))