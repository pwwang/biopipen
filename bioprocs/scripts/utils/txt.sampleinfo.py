if 'txtSampleinfo' not in vars() or not callable (txtSampleinfo):
	from collections import OrderedDict
	class txtSampleinfo(object):
		def __init__(self, sfile):
			if not isinstance(sfile, dict):
				with open(sfile) as f:
					self.data = [line.split('\t') for line in list(filter(None, [l.strip() for l in f.read().splitlines()]))]
				
				self.nrow = len(self.data) - 1
				if self.nrow < 1:
					raise ValueError('No data found in sample info file: %s' % sfile)
				
				self.colnames = self.data.pop(0)
				self.rownames = [row.pop(0) for row in self.data]
				self.ncol     = len(self.data[0]) 
				lencnames     = len(self.colnames)
				if lencnames - self.ncol > 1 or lencnames - self.ncol < 0:
					raise ValueError('Number of columns inconsistent between header and row 1 in sample info file: %s.' % sfile)
				
				if self.ncol == lencnames - 1:
					self.colnames.pop(0)

				if any([len(row)!=self.ncol for row in self.data]):
					raise ValueError('Number of columns inconsistent among rows in sample info file: %s.' % sfile)

				if not set(['Patient', 'Group', 'Batch']) & set(self.colnames):
					raise ValueError('Samaple info file must have one or more columns of Patient, Group, Batch.')
				
			else:
				self.data     = sfile['data']
				self.nrow     = sfile['nrow']
				self.ncol     = sfile['ncol']
				self.rownames = sfile['rownames']
				self.colnames = sfile['colnames']
		
		@staticmethod
		def _alwaysList(ele):
			if not isinstance(ele, tuple) and not isinstance(ele, list):
				if not ele: return []
				return [ele]
			return list(ele)

		def __getitem__(self, key):
			key = self._alwaysList(key)
			if len(key) == 0:
				raise KeyError('Need at least one subscript to extract data.')
			if len(key) > 2:
				raise KeyError('Cannot have >2 subscripts for a 2-dim matrix.')

			rkeystmp = self._alwaysList(key[0])
			rowkeys  = []
			for rowkey in rkeystmp:
				if not isinstance(rowkey, int):
					if not rowkey in self.rownames:
						raise KeyError('Unknown rowname: %s' % rowkey) 
					rowkey = self.rownames.index(rowkey)
				if rowkey < 0 or rowkey >= self.nrow:
					raise KeyError('Index out of range: %s.' % rowkey)
				rowkeys.append(rowkey)
			
			colkeys  = []
			if len(key) == 2:
				ckeystmp = self._alwaysList(key[1])
				for colkey in ckeystmp:
					if not isinstance(colkey, int):
						if not colkey in self.colnames:
							raise KeyError('Unknown colname: %s' % colkey) 
						colkey = self.colnames.index(colkey)
					if colkey < 0 or colkey >= self.ncol:
						raise KeyError('Index out of range: %s.' % colkey)
					colkeys.append(colkey)
			
			if not rowkeys: rowkeys = list(range(self.nrow))
			if not colkeys: colkeys = list(range(self.ncol))

			ret = {}
			ret['data'] = [[r for j, r in enumerate(row) if j in colkeys] for i, row in enumerate(self.data) if i in rowkeys]
			ret['nrow'] = len(rowkeys)
			ret['ncol'] = len(colkeys)
			ret['rownames'] = [row for i, row in enumerate(self.rownames) if i in rowkeys]
			ret['colnames'] = [col for i, col in enumerate(self.colnames) if i in colkeys]
			return txtSampleinfo(ret)

