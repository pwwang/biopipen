# PYPPL REPEAT START: readBase
if 'readBase' not in vars() or not callable(readBase):
	class readBase(object):
		def __init__(self, infile, delimit = '\t', comment = '#', skip = 0):
			if infile.endswith('.gz'):
				import gzip

			self.__dict__ = {
				'meta'   : readMeta(),
				'opener' : gzip.open(infile) if infile.endswith('.gz') else open(infile),
				'delimit': delimit,
				'comment': comment,
				'tell'   : 0
			}
			if skip > 0:
				for _ in range(skip):
					self.opener.readline()
			self.tell = self.opener.tell()

		def _parse(self, line):
			ret = {}
			for i, key in enumerate(self.meta.keys()):
				ret[key] = line[i] if i < len(line) else ''
			return readRecord(**ret)

		def next(self):
			line = self.opener.readline().rstrip('\n')
			while line.startswith(self.comment):
				line = self.opener.readline().rstrip('\n')
			if not line: raise StopIteration()
			return self._parse(line.split(self.delimit))

		def dump(self):
			return [r for r in self]

		def rewind(self):
			self.opener.seek(self.tell)

		def __iter__(self):
			return self

		def __del__(self):
			if self.opener:
				self.opener.close()
# PYPPL REPEAT END: readBase