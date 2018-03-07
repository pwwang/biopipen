# PYPPL REPEAT START: readBedx
if 'readBedx' not in vars() or not callable(readBedx):
	class readBedx(readBase):
		META = [
			('CHR'  , 'Chromosome'),
			('START'  , 'Start position'),
			('END'  , 'End position'),
			('NAME'  , 'Name of the region pair'),
			('SCORE'  , 'The score of the region'),
			('STRAND'  , 'The strand of the region')
		]

		def __init__(self, infile, metaprefix = '##META/', skip = 0, comment = '#', delimit = '\t', xcols = None, headprefix = '#'):
			super(readBedx, self).__init__(infile, skip = skip, comment = comment, delimit = delimit)
			self.meta = readMeta(*readBedx.META)
			if xcols and isinstance(xcols, list):
				self.meta.add(*xcols)
			else:
				for line in self.opener:
					if line.startswith(metaprefix):
						xcol = line[len(metaprefix):].split(':', 1)
						if len(xcol) == 1: xcol.append('')
						mname, mdesc = xcol
						self.meta.add((mname.strip(), mdesc.strip()))
					elif line.startswith(headprefix):
						xcols = line.strip('#\n ').split(delimit)
						self.meta.add(*xcols)
				if len(self.meta.keys()) == 6:
					raise ValueError('You may have to specify extra columns before read the file.\n')
				self.opener.seek(0)
				if skip > 0:
					for _ in range(skip):
						self.opener.readline()
# PYPPL REPEAT END: readBedx