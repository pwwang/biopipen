if 'txtSampleinfo' not in vars() or not callable (txtSampleinfo):
	from collections import OrderedDict
	def txtSampleinfo(sfile):
		info = OrderedDict()
		cols = {}
		firstLine = True
		headers = []
		with open(sfile) as f:
			for line in f:
				line  = line.strip()
				if not line: continue
				parts = line.split("\t")
				if firstLine: 
					headers = parts
					if headers[0] != 'Sample':
						headers.insert(0, 'Sample')
					if not set(headers) & set(['Patient', 'Group', 'Batch']):
						raise ValueError("Headers should be a subset of ['Patient', 'Group', 'Batch']")

					for header in headers: cols[header] = []
					firstLine = False
					continue

				info[parts[0]] = {}
				for i, part in enumerate(parts):
					cols[headers[i]].append(part)
					if i > 0:
						info[parts[0]][headers[i]] = part
				
		return info, cols