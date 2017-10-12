if 'txtFilter' not in vars() or not callable (txtFilter):
	# row filter comes first
	def txtFilter(infile, outfile, cols = [], rfilter = None, header = True, skip = 0, delimit = "\t"):
		if not isinstance(cols, list):
			cols = map(lambda x: x.strip(), cols.split(','))
			cols = [int(c) if c.isdigit() else c for c in cols]
		
		colnames = []
		colindex = cols
		if header:
			with open(infile) as f:
				colnames = f.readline().strip().split(delimit)
				nextline = f.readline().strip("\n").split(delimit)
				if len(nextline) == len(colnames) + 1:
					colnames.insert(0, '')
			for i, col in enumerate(cols):
				if isinstance(col, int): 
					colindex[i] = col
					continue
				if not col in colnames:
					raise ValueError('Unknown column names %s.' % col)
				colindex[i] = colnames.index(col)
		else:
			if not all([isinstance(c, int) for c in cols]):
				raise ValueError('Columns have to be indices if input file has no header: %s.' % cols)
	
		if not colindex and not rfilter and not skip:
			import os
			if os.path.exists(outfile):
				os.remove(outfile)
			os.symlink(infile, outfile)
			return

		import sys
		import csv
		csv.field_size_limit(sys.maxsize)

		with open (infile, 'r') as f, open(outfile, 'w') as fout:
			fcsv = csv.reader(f, delimiter = delimit)
			try:
				for x in range(skip):
					next(fcsv)
				if header: 
					fout.write("\t".join(colnames if not colindex else [colnames[i] for i in colindex]) + "\n")
					next(fcsv)
			except StopIteration:
				pass
			
			for parts in fcsv:
				# do row filter
				if callable(rfilter) and not rfilter(parts): continue
				outs = parts if not colindex else [parts[i] for i in colindex]
				fout.write("\t".join(outs) + "\n")