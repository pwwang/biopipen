if 'writeBedx' not in vars() or not callable(writeBedx):
	class writeBedx(writeBase):
			
		def __init__(self, outfile):
			super(writeBedx, self).__init__(outfile)
			self.meta = readMeta(*readBedx.META)
