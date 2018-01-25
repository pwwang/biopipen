if 'writeBase' not in vars() or not callable(writeBase):
	class writeBase(object):
		def __init__(self, outfile):
			self.meta    = readMeta()
			self.opener  = open(outfile, 'w')

		def writeMeta(self, prefix = '##META/'):
			for key in self.meta.keys():
				self.opener.write(prefix + key + ': %s' % getattr(self.meta, key) + '\n')
		
		def writeHead(self, prefix = '#', delimit = '\t'):
			self.opener.write(prefix + delimit.join(self.meta.keys()) + '\n')

		def write(self, x, delimit = '\t'):
			outs = []
			for key in self.meta.keys():
				outs.append(str(getattr(x, key)))
			self.opener.write(delimit.join(outs) + '\n')

		def __del__(self):
			if self.opener:
				self.opener.close()