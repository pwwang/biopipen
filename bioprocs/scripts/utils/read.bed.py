if 'readBed' not in vars() or not callable(readBed):
	class readBed(readBase):
		META = [
			('CHR', 'Chromosome'),
			('START', 'Start position'),
			('END', 'End position'),
			('NAME', 'Name of the region pair'),
			('SCORE', 'The score of the region'),
			('STRAND', 'The strand of the region')
		]

		def __init__(self, infile, skip = 0, comment = '#', delimit = '\t'):
			super(readBed, self).__init__(infile, skip = skip, comment = comment, delimit = delimit)
			self.meta = readMeta(
				*readBed.META
			)
			self.index = 1

		def _parse(self, line):
			r = super(readBed, self)._parse(line)
			if not r.NAME: 
				r.NAME = 'BED' + str(self.index)
				self.index += 1
			if not r.SCORE: r.SCORE   = '0'
			if not r.STRAND: r.STRAND = '+'
			return r