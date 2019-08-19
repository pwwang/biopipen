import sys
from os import path
from pyppl import Box
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
statfile = {{o.statfile | quote}}
bcftools = {{args.bcftools | quote}}
params   = {{args.params | repr}}
nthread  = {{args.nthread | repr}}
keep     = {{args.keep | repr}}
gz       = {{args.gz | repr}}
stat     = {{args.stat | repr}}
include  = {{args.include | repr}}
exclude  = {{args.exclude | repr}}

shell.load_config(bcftools = bcftools)

params.O = 'z' if gz else 'v'

def runBcftools(command, *args, **kwargs):
	cmd = getattr(shell.bcftools, command)(*args, **kwargs).cmd
	sys.stderr.write("\n%s RUNNING %s\n%s\n%s\n\n" % ("-" * 40, "-" * 40, cmd, "-" *  89))

def normExpr(expr):
	if not expr:
		return {}
	if isinstance(expr, list):
		return {("Filter%s" % i+1): ex for i, ex in enumerate(expr)}
	if isinstance(expr, dict):
		return expr
	return {"Filter1": expr}

# get the normalized filters
# if not include and not exclude:
# 	raise ValueError('At least one expression for `args.include` or `args.exclude` is required.')

if include and exclude:
	raise ValueError('We can only handle `include` or `exclude` expression at one time, please rephrase one of them into the same type.')

includes = normExpr(include)
excludes = normExpr(exclude)

flag = 'include' if includes else 'exclude'
exprs = includes or excludes
for filter_name, filter_expr in exprs.items():
	# run each filter at one time, because we want to give each filter a different name
	sys.stderr.write("Running filter: %s ...\n" % filter_name)
	tmpfile = outfile + '.' + filter_name
	params_one = params.copy()
	params_one[flag] = filter_expr
	params_one.s = filter_name
	params_one._ = infile
	params_one.o = tmpfile
	params_one.threads = nthread
	runBcftools('filter', **params_one)
	infile = tmpfile

# stat file
if stat:
	annfile = statfile + '.ann'
	with open(infile) as fin, open(annfile, 'w') as fann, open(statfile, 'w') as fstat:
		filters = {'PASS': 0}
		fann.write('"": "Variants may share multiple filters. Variants with PASS are desired ones and exclusive with other filters."\n')
		for line in fin:
			# ##FILTER=<ID=SIMPLEREPEAT,Description="Set if true: INFO/SIMPLEREPEAT">
			if line.startswith('##FILTER=<ID='):
				# SIMPLEREPEAT: "Set if true: INFO/SIMPLEREPEAT"
				line = line.strip()
				line = line[13:-1].replace(',Description=', ': ')
				fann.write(line + '\n')
			elif not line.startswith('#'):
				filts = line.split('\t')[6].split(';')
				for filt in filts:
					if filt == '.':
						filt = 'PASS'
					if filt not in filters:
						filters[filt] = 0
					filters[filt] += 1
		fstat.write(
			'\t'.join(
				[''] + 
				[f for f in filters.keys() if f != 'PASS'] + 
				['PASS']
			) + '\n')
		fstat.write(
			'\t'.join(
				[path.splitext(path.basename(infile))[0]] + 
				[str(count) for f, count in filters.items() if f != 'PASS'] +
				[str(filters['PASS'])]
			) + '\n')
else:
	open(statfile, 'w').close()

if not keep:
	runBcftools('view', f = '.,PASS', o = outfile, O = params.O, threads = nthread, _ = infile)
else:
	# don't do just simple move, in case of inconsistency of gz options of input and output
	runBcftools('convert', o = outfile, O = params.O, threads = nthread, _ = infile)
