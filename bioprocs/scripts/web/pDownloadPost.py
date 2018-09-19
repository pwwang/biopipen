import requests

data = {{i.data}}
r    = requests.post({{i.url | quote}}, data)

rc   = int(r.status_code)
if rc != 200:
	raise Exception('Failed to POST to {{i.url}}, status code is %s' % rc)

with open({{o.outfile | quote}}, 'w') as fout:
	fout.write(r.text)