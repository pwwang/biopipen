from bioprocs.utils.tsvio import TsvReader, TsvWriter

files   = {{in.infiles}}
lenfs   = len(files)
inopts  = {{args.inopts}}
outfile = {{out.outfile | quote}}
inopts_each = []
maxopen = {{args.maxopen}}
for i in range(lenfs):
	skip = 0
	if 'skip' in inopts:
		if not isinstance(inopts['skip'], list):
			skip = inopts['skip']
		elif i < len(inopts['skip']):
			skip = inopts['skip'][i]
	delimit = '\t'
	if 'delimit' in inopts:
		if not isinstance(inopts['delimit'], list):
			delimit = inopts['delimit']
		elif i < len(inopts['delimit']):
			delimit = inopts['delimit'][i]
	comment = '#'
	if 'comment' in inopts:
		if not isinstance(inopts['comment'], list):
			comment = inopts['comment']
		elif i < len(inopts['comment']):
			comment = inopts['comment'][i]
	ftype = ''
	if 'ftype' in inopts:
		if not isinstance(inopts['ftype'], list):
			ftype = inopts['ftype']
		elif i < len(inopts['ftype']):
			ftype = inopts['ftype'][i]
	cnames = ''
	if 'cnames' in inopts:
		if not isinstance(inopts['cnames'], list):
			cnames = inopts['cnames']
		elif i < len(inopts['cnames']):
			cnames = inopts['cnames'][i]
	head = ''
	if 'head' in inopts:
		if not isinstance(inopts['head'], list):
			head = inopts['head']
		elif i < len(inopts['head']):
			head = inopts['head'][i]
	inopts_each.append({'skip': skip, 'delimit': delimit, 'comment': comment, 'ftype': ftype, 'cnames': cnames, 'head': head})

def getReaders(fs, baseidx):
	readers = []
	for i,infile in enumerate(fs):
		reader = TsvReader(infile, **inopts_each[i + baseidx])
		if not reader.meta:
			reader.autoMeta()
		readers.append(reader)
	return readers
	
outopts = {
	'ftype':         '',
	'head'         : False, 
	'headPrefix'   : '', 
	'headDelimit'  : '\t', 
	'headTransform': None, 
	'delimit'      : '\t'
}
outopts.update({{args.outopts}})
writer = TsvWriter(outfile, ftype = outopts['ftype'], delimit = outopts['delimit'])

metaUpdated = False
for i in xrange(0, len(files), maxopen):
	fs = files[i:i + maxopen]
	readers = getReaders(fs, i)
	
	# empty file has no meta
	reader = [reader for reader in readers if reader.meta and not (len(reader.meta) == 1 and reader.meta.items() == [('', None)])]
	if not reader: continue
	reader = reader[0]
	if not metaUpdated:
		metaUpdated = True
		writer.meta.update(reader.meta)
		if outopts['head']:
			writer.writeHead(prefix = outopts['headPrefix'], delimit = outopts['headDelimit'], transform = outopts['headTransform'])
		
	for reader in readers:
		for r in reader:
			writer.write(r)
