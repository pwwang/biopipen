# PYPPL REPEAT START: readHead
if 'readHead' not in vars() or not callable(readHead):
	
	class readHead(readBase):

		def __init__(self, infile, comment = '#', delimit = '\t', skip = 0):
			super(readHead, self).__init__(infile, skip = skip, comment = comment, delimit = delimit)
			self.meta = readMeta()
			
			header = self.opener.readline().strip('#\n ').split(delimit)
			self.tell = self.opener.tell()
			row1   = self.opener.readline().strip().split(delimit)
			if len(row1) == len(header) + 1:
				header.insert(0, 'ROWNAMES')
			self.meta.add(*header)
			
			self.rewind()
# PYPPL REPEAT END: readHead