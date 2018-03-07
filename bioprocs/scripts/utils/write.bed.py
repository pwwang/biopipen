# PYPPL REPEAT START: writeBed
if 'writeBed' not in vars() or not callable(writeBed):
	class writeBed(writeBase):
			
		def __init__(self, outfile):
			super(writeBed, self).__init__(outfile)
			self.meta = readMeta(*readBed.META)
# PYPPL REPEAT END: writeBed