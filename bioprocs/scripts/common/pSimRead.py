from os import path
from simread import SimRead

fout = open({{out.outfile | quote}}, 'w')
files   = {{in.infiles}}
r       = SimRead(*files, skip = {{args.skip}}, delimit = {{args.delimit}}, gzip = {{args.gzip}})
r.do    = {{args.do}}
{% if args.match %}
r.match = {{args.match}}
{% endif %}
r.run()
fout.close()