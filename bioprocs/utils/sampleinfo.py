from bioprocs.utils import alwaysList
from bioprocs.utils.tsvio import TsvReader, TsvRecord

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
		ret = []
		for row in self.data:
			if sample and isinstance(sample, list) and row[self.samcol] not in sample:
				continue
			if sample and not isinstance(sample, list) and row[self.samcol] != sample:
				continue
			if patient and isinstance(patient, list) and row.Patient not in patient:
				continue
			if patient and not isinstance(patient, list) and row.Patient != patient:
				continue
			if group and isinstance(group, list) and row.Group not in group:
				continue
			if group and not isinstance(group, list) and row.Group != group:
				continue
			if batch and isinstance(batch, list) and row.Batch not in batch:
				continue
			if batch and not isinstance(batch, list) and row.Batch != batch:
				continue
			if not get:
				ret.append(row)
			else:
				if not isinstance(get, list):
					get = [get]
				r = TsvRecord()
				for g in get:
					if g == 'Sample':
						g = self.samcol
					r[g] = row[g]

				ret.append(r)
		return ret

	def dataframe(self, data = None):
		data = data or self.data
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
