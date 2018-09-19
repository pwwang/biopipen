from eutils import Client
from bioprocs.utils.tsvio import TsvWriter, TsvRecord
from time import sleep

term     = {{i.term | repr}}
outfile  = {{o.outfile | quote}}
prog     = {{args.prog | quote}}
apikey   = {{args.apikey | quote}} or None
db       = {{args.db | quote}}
joiner   = {{args.joiner | quote}}
record   = {{args.record}}
writer   = TsvWriter(outfile)
sleepsec = {{args.sleep | repr}}

def writerResults(rets):
	keys = []
	if rets:
		ret  = rets[0]
		keys = [key for key in dir(ret) if not key.startswith('_')]
		writer.meta.add(*keys)
		writer.writeHead()
	for ret in rets:
		r = TsvRecord()
		for key in keys:
			try:
				r[key] = getattr(ret, key)
			except:
				r[key] = ''
			if isinstance(r[key], (list, tuple)):
				try:
					r[key] = joiner.join([str(_) for _ in r[key]])
				except UnicodeEncodeError:
					r[key] = joiner.join([_.encode('utf-8') for _ in r[key]])
			else:
				try:
					r[key] = str(r[key])
				except UnicodeEncodeError:
					r[key] = r[key].encode('utf-8')
					
		if callable(record):
			r = record(r)
		elif record is not None:
			raise ValueError('Unknown record transform function (args.record).')
		if r:
			writer.write(r)

client = Client(api_key = apikey)
if prog == 'esearch':
	sret  = client.esearch(db = db, term = term)
	try:
		error = list(sret._xml_root.find('ErrorList').iterchildren())
	except:
		error = None
	
	print sret.count if not error else 0
	
	if not sret.ids:
		rets = []
	else:
		rets = client.efetch(db = db, id = sret.ids)
		rets = list(iter(rets))
	writerResults(rets)

else:
	fetches = client.efetch(db = db, id = term)
	key     = [key for key in dir(fetches) if not key.startswith('_')][0]
	rets    = getattr(fetches, key)
	writerResults(rets)

sleep (sleepsec)