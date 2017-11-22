from os import path
from simread import SimRead
from gzip import open as gzopen

fout = open({{out.outfile | quote}}, 'w')
files   = {{in.infiles}}
skip    = {{args.skip}}
gzip    = {{args.gzip | lambda x: '"auto"' if x == 'auto' else x}}
# get header before run
{% if args.usehead | lambda x: x is not None %}
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
{% endif %}

r       = SimRead(*files, skip = skip, delimit = {{args.delimit}}, gzip = gzip)
r.do    = {{args.do}}
{% if args.match %}
r.match = {{args.match}}
{% endif %}

r.run()
fout.close()