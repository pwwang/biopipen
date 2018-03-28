from bioprocs.utils.tsvio import TsvReader, TsvWriter

files   = {{in.infiles}}
lenfs   = len(files)
inopts  = {{args.inopts}}
outfile = {{out.outfile | quote}}
inopts_each = []
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
	inopts_each.append({'skip': skip, 'delimit': delimit, 'comment': comment, 'ftype': ftype})
	
readers = []
for i,infile in enumerate(files):
	reader = TsvReader(infile, **inopts_each[i])
	if not reader.meta:
		reader.autoMeta()
	readers.append(reader)

outopts = {
	'head'         : False, 
	'headPrefix'   : '', 
	'headDelimit'  : '\t', 
	'headTransform': None, 
	'delimit'      : '\t'
}
outopts.update({{args.outopts}})
writer = TsvWriter(outfile, delimit = outopts['delimit'])
writer.meta.update(readers[0].meta)
if outopts['head']:
	writer.writeHead(prefix = outopts['headPrefix'], delimit = outopts['headDelimit'], transform = outopts['headTransform'])
	
for reader in readers:
	for r in reader:
		writer.write(r)
