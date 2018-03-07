# PYPPL REPEAT START: readBedpe
if 'readBedpe' not in vars() or not callable(readBedpe):
	class readBedpe(readBase):
		META = [
			('CHR1'    , 'Chromosome 1'),
			('START1'  , 'Start position 1'),
			('END1'    , 'End position 1'),
			('CHR2'    , 'Chromosome 2'),
			('START2'  , 'Start position 2'),
			('END2'    , 'End position 2'),
			('NAME'    , 'Name of the region pair'),
			('SCORE'   , 'The score of the region pair'),
			('STRAND1' , 'The strand of the first pair'),
			('STRAND2' , 'The strand of the second pair')
		]

		def __init__(self, infile, skip = 0, comment = '#', delimit = '\t'):
			super(readBedpe, self).__init__(infile, skip = skip, comment = comment, delimit = delimit)
			self.meta = readMeta(
				*readBedpe.META
			)
# PYPPL REPEAT END: readBedpe