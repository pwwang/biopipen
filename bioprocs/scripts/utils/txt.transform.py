if 'txtTransform' not in vars() or not callable (txtTransform):
	def txtTransform(infile, outfile, cols = [], transform = None, header = True, skip = 0, delimit = "\t"):

		if not isinstance(cols, list):
			cols.strip().split(',')
		cols = [c if isinstance(c, int) else int(c) if c.isdigit() else c for c in cols]
		
		import sys, csv
		csv.field_size_limit(sys.maxsize)
		with open (infile, 'r') as f, open(outfile, 'w') as fout:
			for _ in range(skip):
				fout.write(f.readline())

			fcsv = csv.reader(f, delimiter = delimit)
			if header:
				headers = fcsv.next()
				cols    = [c if isinstance(c, int) else headers.index(c) for c in cols] if cols else range(len(headers))
				fout.write(delimit.join([headers[c] for c in cols]) + "\n")
			for row in fcsv:
				row = [row[c] for c in cols]
				row = transform(row) if transform else row
				fout.write(delimit.join(row) + "\n")