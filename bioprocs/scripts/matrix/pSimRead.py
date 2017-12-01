from os import path
from simread import SimRead
from gzip import open as gzopen

fout = open({{out.outfile | quote}}, 'w')
files   = {{in.infiles}}
skip    = {{args.skip}}
gzip    = {{args.gzip | lambda x: '"auto"' if x == 'auto' else x}}
# get header before run
{% if args.usehead | lambda x: isinstance(x, int) %}
if not isinstance(skip, list):
	skip = [skip] * len(files)
if not isinstance(gzip, list):
	gzip = [gzip] * len(files)
uhskip   = 0 if len(skip) < {{args.usehead}} + 1 else skip[{{args.usehead}}]
uhfile   = files[{{args.usehead}}]
uhopen   = open 
if len(gzip) > {{args.usehead}} \
	and ( \
		gzip[{{args.usehead}}] == True \
		or (gzip[{{args.usehead}}] == 'auto' and uhfile.endswith('.gz')) \
	):
	uhopen = gzopen

with uhopen(uhfile) as f:
	for _ in range(uhskip):
		fout.write(f.readline())
{% elif args.usehead %}
usehead = {% if args.usehead | lambda x: isinstance(x, list) %}{{args.usehead}}{% else %}{{args.usehead | .split(',') | lambda x: [s.strip() for s in x]}}{% endif %}
fout.write('\t'.join(usehead) + '\n')
{% endif %}

# used for do func
def writeline(line, end = '\n'):
	fout.write(line + end)

def writelist(parts, delimit = '\t', end = '\n'):
	fout.write(delimit.join(parts) + end)

# used for match func
def compare(part1, part2, reverse = False):
	if not reverse:
		return 0 if part1 < part2 else 1 if part1 > part2 else -1
	else:
		return 0 if part1 > part2 else 1 if part1 < part2 else -1

r       = SimRead(*files, skip = skip, delimit = {% if args.delimit | lambda x: isinstance(x, list) %}{{args.delimit}}{% else %}{{ args.delimit | quote }}{% endif %}, gzip = gzip)
r.do    = {{args.do}}
{% if args.match %}
r.match = {{args.match}}
{% endif %}

r.run()
fout.close()