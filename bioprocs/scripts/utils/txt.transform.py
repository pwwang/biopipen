if 'txtTransform' not in vars() or not callable (txtTransform):
	def txtTransform(infile, outfile, transformer = None, header = True, delimit = "\t", comment = "#"):

		import sys
		import csv
		csv.field_size_limit(sys.maxsize)

		with open (infile, 'r') as f, open(outfile, 'w') as fout:
			fcsv = csv.reader(f, delimiter = delimit)
			
			i = 0
			for parts in fcsv:
				# do row filter
				if i == 0 and header: continue
				if parts[0].startswith(comment): continue
				outs = parts if not callable(transformer) else transformer(parts)
				outs = [str(o) for o in outs]
				fout.write("\t".join(outs) + "\n")
				i += 1