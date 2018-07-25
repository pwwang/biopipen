from bioprocs.utils import alwaysList
from bioprocs.utils.tsvio import TsvReader, TsvRecord

class SampleInfoException(Exception):
	pass

class SampleInfo (object):

	NORMAL_GROUP = ['NORMAL', 'HEALTH', 'BEFORE', 'BLOOD', 'PRE', 'CONTROL']
	TUMOR_GROUP  = ['TUMOR', 'DISEASE', 'AFTER', 'TISSUE', 'POST', 'TREATMENT']

	def __init__(self, sfile):
		reader    = TsvReader(sfile, ftype = 'head')
		self.samcol = reader.meta.keys()[0]
		if self.samcol == 'ROWNAMES':
			self.samcol = 'Sample'
			del reader.meta['ROWNAMES']
			reader.meta.prepend('Sample')

		self.data = reader.dump()
		self.nrow = len(self.data)
		self.ncol = len(reader.meta)
		self.colnames = reader.meta.keys()
		self.rownames = [row[self.samcol] for row in self.data]

		expectColnames = ['Sample', 'Patient', 'Group', 'Batch']
		if not set(expectColnames) & set(self.colnames):
			raise SampleInfoException('Unexpected column names: %s.' % str(self.colnames))

	def select(self, sample = None, patient = None, group = None, batch = None, get = None):
		"""
		Select records by different standards
		"""
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

	def toChannel(self, datadir, paired = False, raiseExc = True):
		from os import path
		if not paired:
			samples = self.select(get = 'Sample')
			return [path.join(datadir, r.Sample) for r in samples]
		else:
			data     = self.select(get = ['Sample', 'Patient'])
			patients = [r.Patient for r in data]
			patients = sorted(set(patients), key = patients.index)
			ret      = []
			for patient in patients:
				records = self.select(patient = patient, get = ['Sample', 'Group'])
				if len(records) != 2:
					if raiseExc:
						raise SampleInfoException('Expect 2 samples for paired channel for patient: %s' % patient)
					else:
						continue
				record1, record2 = records
				# make sure record1 is tumor and record2 is normal
				if  record1.Group.upper() in SampleInfo.TUMOR_GROUP + SampleInfo.NORMAL_GROUP and \
					record2.Group.upper() in SampleInfo.TUMOR_GROUP + SampleInfo.NORMAL_GROUP:
					if record1.Group.upper() in SampleInfo.NORMAL_GROUP:
						record1, record2 = record2, record1
				else:
					common = path.commonprefix([record1.Sample, record2.Sample])
					if common and record1.Sample < record2.Sample:
						record1, record2 = record2, record1
				ret.append((path.join(datadir, record1.Sample), path.join(datadir, record2.Sample)))
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


