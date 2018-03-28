from bioprocs.utils import alwaysList
from bioprocs.utils.tsvio import TsvReader

class SampleInfoException(Exception):
	pass	

class SampleInfo (object):
	def __init__(self, sfile):
		reader    = TsvReader(sfile, ftype = 'head')
		self.data = reader.dump()
		self.nrow = len(self.data)
		self.ncol = len(reader.meta) 
		self.colnames = reader.meta.keys()
		self.samcol   = self.colnames[0]
		self.rownames = [row[self.samcol] for row in self.data]
		
		expectColnames = ['ROWNAMES', 'Sample', 'Patient', 'Group', 'Batch']
		if not set(expectColnames) & set(self.colnames):
			raise SampleInfoException('Unexpected column names: %s.' % str(self.colnames))
	
	def select(self, sample = None, patient = None, group = None, batch = None, get = None):
		return [
			row if not get else row[get if get != 'Sample' else self.samcol] for row in self.data if \
				(not sample  or row[self.samcol] == sample) and \
				(not patient or 'Patient' in row and row.Patient == patient) and \
				(not group   or 'Group' in row and row.Group == group) and \
				(not batch   or 'Batch' in row and row.Batch == batch)  
		]
		
	def dataframe(self, data):
		from pandas import DataFrame
		rownames  = []
		dataframe = {}
		for row in data:
			for key, val in row.items():
				if key == self.samcol:
					rownames.append(val)
				elif key in dataframe:
					dataframe[key].append(val)
				else:
					dataframe[key] = [val]
		return DataFrame(data = dataframe, index = rownames)
		
