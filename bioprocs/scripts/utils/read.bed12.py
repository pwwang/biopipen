# PYPPL REPEAT START: readBed12
if 'readBed12' not in vars() or not callable(readBed12):
	class readBed12(readBase):
		META = [
			('CHR'         , 'Chromosome'),
			('START'       , 'Start position'),
			('END'         , 'End position'),
			('NAME'        , 'Name of the region pair'),
			('SCORE'       , 'The score of the region'),
			('STRAND'      , 'The strand of the region'),
			('THICKSTART'  , 'The starting position at which the feature is drawn thickly'),
			('THICKEND'    , 'The ending position at which the feature is drawn thickly'),
			('ITEMRGB'     , 'An RGB value of the form R,G,B'),
			('BLOCKCOUNT'  , 'The number of blocks (exons) in the BED line'),
			('BLOCKSIZES'  , 'A comma-separated list of the block sizes'),
			('BLOCKSTARTS' , 'A comma-separated list of block starts')
		]

		def __init__(self, infile, skip = 0, comment = '#', delimit = '\t'):
			super(readBed12, self).__init__(infile, skip = skip, comment = comment, delimit = delimit)
			self.meta = readMeta(
				*readBed12.META
			)
# PYPPL REPEAT END: readBed12