import requests

data = {{in.data}}
r    = requests.post({{in.url | quote}}, data)

rc   = int(r.status_code)
if rc != 200:
	raise Exception('Failed to POST to {{in.url}}, status code is %s' % rc)

with open({{out.outfile | quote}}, 'w') as fout:
	fout.write(r.text)