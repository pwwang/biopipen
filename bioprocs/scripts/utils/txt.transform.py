if 'txtTransform' not in vars() or not callable (txtTransform):
	def txtTransform(infile, outfile, transformer = None, header = True, delimit = "\t"):

		import sys
		import csv
		csv.field_size_limit(sys.maxsize)

		with open (infile, 'r') as f, open(outfile, 'w') as fout:
			fcsv = csv.reader(f, delimiter = delimit)
			try:
				if header: 
					cnames = next(fcsv)
					fout.write("\t".join(cnames) + "\n")
			except StopIteration:
				pass
			
			for parts in fcsv:
				# do row filter
				outs = parts if not callable(transformer) else transformer(parts)
				fout.write("\t".join(outs) + "\n")