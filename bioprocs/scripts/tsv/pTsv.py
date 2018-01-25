{# Load the classes #}
{% for k, v in read.items() %}
{{v.py}}
{% endfor %}
{% for k, v in write.items() %}
{{v.py}}
{% endfor %}

inopts  = {{args.inopts}}
outopts = {{args.outopts}}

{# get the reader class #}
{% if args.inmeta | lambda x: isinstance(x, list) %}
metadata = {{args.inmeta}}
reader = readBase({{in.infile | quote}}, **inopts)
reader.meta.add(*metadata)
{% elif args.inmeta | lambda x: isinstance(x, dict) %}
metadata = {{args.inmeta}}.items()
reader = readBase({{in.infile | quote}}, **inopts)
reader.meta.add(*metadata)
{% elif args.inmeta %}
inmeta = {{args.inmeta | quote}}
inmeta = inmeta[0].upper() + inmeta[1:]
reader = locals()['read' + inmeta](
	infile  = {{in.infile | quote}},
	**inopts
)
{% else %}
inmeta = {{in.infile | ext | [1:] | quote}}
inmeta = inmeta[0].upper() + inmeta[1:]
reader = locals()['read' + inmeta](
	infile  = {{in.infile | quote}},
	**inopts
)
{% endif %}

{# get the writer class #}
writemeta  = True
writehead  = True
metaprefix = '##META/'
headprefix = '#'
if 'meta' in outopts:
	writemeta = outopts['meta']
	del outopts['meta']
if 'head' in outopts:
	writemeta = outopts['head']
	del outopts['head']
if metaprefix in outopts:
	metaprefix = outopts['metaprefix']
	del outopts['metaprefix']
if headprefix in outopts:
	headprefix = outopts['headprefix']
	del outopts['headprefix']
	
writer = writeBase({{out.outfile | quote}})
{% if args.outmeta | lambda x: not x %}
writer.meta.borrow(reader.meta)
{% elif args.outmeta | lambda x: isinstance(x, list) %}
writer.meta.add(*{{args.outmeta}})
{% elif args.outmeta | lambda x: isinstance(x, dict) %}
writer.meta.add(*{{args.outmeta}}.items())
{% else %}
outmeta = {{args.outmeta | lambda x: x[0].upper() + x[1:] | quote}}
writer  = locals()['write' + outmeta]({{out.outfile | quote}})
xcols   = [] if 'xcols' not in outopts else outopts['xcols']
if xcols: writer.meta.add(*xcols)
{% endif %}

opshelper = [line for line in {{args.opshelper | quote}}.splitlines() if line]
while opshelper and (all([line and line[0] == ' ' for line in opshelper]) or all([line and line[0] == '\t' for line in opshelper])):
	opshelper = [line[1:] for line in opshelper]
exec('\n'.join(opshelper))

{# write meta #}
if writemeta:
	writer.writeMeta(prefix = metaprefix)

if writehead:
	writer.writeHead(prefix = headprefix, delimit = outopts['delimit'])

ops = {{args.ops}}
for r in reader:
	if not ops:
		writer.write(r, delimit = outopts['delimit'])
	else:
		r = ops(r)
		if r:
			writer.write(r, delimit = outopts['delimit'])
