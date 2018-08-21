bcode   = {{args.bcode | repr}}
acode   = {{args.acode | repr}}
infile  = {{in.infile | quote}}
outfile = {{out.outfile | quote}}
scale   = {{args.scale | dict}}

with open(outfile, 'w') as f:
	if not isinstance(bcode, list):
		bcode = [bcode]
	f.write('# Before image inserted:\n')
	for bc in bcode:
		f.write(bc + '\n')
	f.write('\n')
	
	f.write("# Inserting image: \n")
	scaleargs = ''
	for k,v in scale.items():
		scaleargs = ', {k} = {v}'.format(k = k, v = v)
	f.write("doc.add_picture({pic}{scale})\n".format(
		pic = repr(infile),
		scale = scaleargs
	))
	f.write("\n")

	f.write("# After image inserted:\n")
	if not isinstance(acode, list):
		acode = [acode]
	for ac in acode:
		f.write(ac + '\n')


