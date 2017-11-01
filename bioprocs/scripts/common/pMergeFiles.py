from os import rename, remove
files = {{in.infiles}}
lfile = len(files)

{% if args.skip | lambda x: not x %}
skip = [0] * lfile
{% elif args.skip | lambda x: not isinstance(x, list) %}
skip = [{{args.skip}}] * lfile
{% else %}
skip = [0] * lfile
skip[:len({{args.skip}})] = {{args.skip}}
{% endif %}

{% if args.comment | lambda x: not x %}
comment = ['#'] * lfile
{% elif args.comment | lambda x: not isinstance(x, list) %}
comment = [{{args.comment | quote}}] * lfile
{% else %}
comment = ['#'] * lfile
comment[:len({{args.comment | quote}})] = {{args.comment | quote}}
{% endif %}

{% if args.header | lambda x: x is None %}
header = [False] * lfile
{% elif args.header | lambda x: not isinstance(x, list) %}
header = [{{args.header}}] * lfile
{% else %}
header = [False] * lfile
header[:len({{args.header}})] = {{args.header}}
{% endif %}

headerRow = None
with open("{{out.outfile}}.nohead", 'w') as fout:
	for i, fn in enumerate(files):
		with open(fn) as f:
			for line in f:
				for _ in range(skip[i]): continue
				if header[i] and not headerRow:
					headerRow = line
					continue
				if line.startswith(comment[i]):
					continue
				fout.write(line)
if not headerRow:
	rename("{{out.outfile}}.nohead", {{out.outfile | quote}})
else:
	with open("{{out.outfile}}.nohead") as f, open({{out.outfile | quote}}, 'w') as fout:
		fout.write(headerRow)
		fout.write(f.read())
	remove("{{out.outfile}}.nohead")
		