import os #, shutil
from subprocess import check_output

{{downloadCurl}}
{{runcmd}}
{{txtFilter}}
{{txtTransform}}

url = {{args.urls | json}}["{{i.in}}"]
tmpdir = "{{job.outdir}}/tmp"
if not os.path.exists(tmpdir):
	os.makedirs(tmpdir)
	
downfile = os.path.join(tmpdir, 'downloaded')
downloadCurl(url, downfile, {{args.username | quote}}, {{args.password | quote}}, {{args.curl | quote}})

# determin file types
output = check_output(['file', downfile])
if 'gzip' in output:
	ugfile = downfile + '.ungz'
	with open(ugfile, 'w') as f:
		f.write(check_output(['gunzip', downfile, '-c']))
	downfile = ugfile
elif 'Zip' in output:
	zipdir = os.path.join(tmpdir, '_unzipped')
	import zipfile, glob
	zipref = zipfile.ZipFile(downfile, 'r')
	zipref.extractall(zipdir)
	zipref.close()
	downfile = glob.glob(os.path.join(zipdir, '*'))[0]

outfile = "{{o.outfile}}"[:-3] if {{args.gz | lambda x: bool(x)}} else "{{o.outfile}}"
txtFilter(downfile, outfile + '.filtered', {{args.cols | lambda x: [] if not x else __import__('json').dumps(x) if isinstance(x, list) else '"'+ x +'"'}}, {{args.rowfilter | lambda x: 'False' if not x else x}}, {{args.header}}, {{args.skip}}, {{args.delimit | quote}})

{% if i['in'] == 'KEGG_2016_gmt' %}
def transformer(parts):
	return [ p.replace('/', '-') if i == 0 else p[:-4] if p.endswith(',1.0') else p for i,p in enumerate(parts) ]

tmpfile = outfile + '.kegg'
txtTransform(outfile + '.filtered', tmpfile, transform = transformer, header = {{args.header}}, delimit = {{args.delimit | quote}})
os.rename(tmpfile, outfile + '.filtered')
{% endif %}

{% if args.transform %}
txtTransform(outfile + '.filtered', outfile, transform = {{args.transform}}, header = {{args.header}}, delimit = {{args.delimit | quote}})
{% else %}
os.rename(outfile + '.filtered', outfile)
{% endif %}

if {{args.gz | bool}}:
	runcmd ('gz "%s"' % outfile)

#shutil.rmtree(tmpdir)