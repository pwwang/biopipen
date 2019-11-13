from os import path
from pyppl.utils import alwaysList
from bioprocs.utils.tsvio2 import TsvReader, TsvRecord

class SampleInfoException(Exception):
	pass

class SampleInfo (object):

	NORMAL_GROUP = ['NORMAL', 'HEALTH', 'BEFORE', 'BLOOD', 'PRE', 'CONTROL']
	TUMOR_GROUP  = ['TUMOR', 'DISEASE', 'AFTER', 'TISSUE', 'POST', 'TREATMENT']

	def __init__(self, sfile):
		reader    = TsvReader(sfile)
		self.samcol = reader.cnames[0]
		if self.samcol == 'ROWNAMES':
			self.samcol = 'Sample'
			reader.cnames[0] = 'Sample'

		self.data = reader.dump()
		self.nrow = len(self.data)
		self.ncol = len(reader.cnames)
		self.colnames = reader.cnames
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


class SampleInfo2(object):

	def toChannel(self, datadir, paired = False, raiseExc = True):
		if not paired:
			samples = self.getSamples()
			return [path.join(datadir, sample) for sample in samples]

		paired_samples =  self.getPairedSamples()
		ret = []
		for sample1, sample2 in paired_samples:
			if paired is not True and paired.upper() not in sample1.upper():
				sample1, sample2 = sample2, sample1
			ret.append((path.join(datadir, sample1), path.join(datadir, sample2)))
		return ret

	def _read(self, sifile):
		standard_cnames = ["", "Sample", "Patient", "Group", "Batch"]
		reader          = TsvReader(sifile)

		self.cnames     = reader.cnames
		if not self.cnames:
			raise SampleInfoException('Headers for sample information file is required.')

		if any(cname not in standard_cnames for cname in  self.cnames):
			raise SampleInfoException('Headers should be a subset of {!r}'.format(', '.join(standard_cnames)))

		if "" in self.cnames:
			self.cnames[self.cnames.index("")] = "Sample"

		self.mat = reader.dump()

	def __init__(self, sifile, checkPaired = False):
		self.mat    = None
		self.cnames = []
		self._read(sifile)
		if checkPaired and "Patient" in self.cnames:
			for patient in self.allPatients():
				if len(self.getSamples(by = 'Patient', value = patient)) != 2:
					raise SampleInfoException('Expect paired comparisons, but Patient {!r} has # samples other than 2.'.format(patient))

	def allPatients(self):
		if 'Patient' not in self.cnames:
			return None
		patients = [r.Patient for r in self.mat]
		ret = []
		for patient in patients:
			if patient not in ret:
				ret.append(patient)
		return ret

	def allGroups(self):
		allgroups = [r.Group for r in self.mat]
		ret = []
		for group in allgroups:
			if group not in ret:
				ret.append(group)
		return ret

	def allSamples(self, unique = False, datadir = None):
		samples = self.getSamples()
		if unique:
			ret = []
			for sample in samples:
				if sample not in ret:
					ret.append(sample)
			return ret
		return [path.join(datadir, sample) for sample in samples] if datadir else samples

	def getSamples(self, by = None, value = None, returnAll = False):
		if by and by not in self.cnames:
			raise SampleInfoException('{!r} is not a valid column name.'.format(by))
		if not by:
			return [r if returnAll else r.Sample for r in self.mat]
		return [r if returnAll else r.Sample for r in self.mat if r[by] == value]

	def sampleInfo(self, sample, info = None):
		if not info:
			return [r for r in self.mat if r.Sample == sample]
		if info not in self.cnames:
			raise SampleInfoException('{!r} is not a valid column name.'.format(info))
		return [r[info] for r in self.mat if r.Sample == sample]

	def select(self, samples):
		self.mat = [r for r in self.mat if r.Sample in samples]

	def isPaired(self):
		allpatients = self.allPatients()
		if not allpatients:
			return False
		for patient in allpatients:
			if len(self.getSamples(by = 'Patient', value = patient)) != 2:
				return False
		return True

	def getPairedSample(self, sample):
		patient = self.sampleInfo(sample, info = 'Patient')[0]
		samples = self.getSamples(by = 'Patient', value = patient)
		return samples[1 - samples.index(sample)]

	def getPairedSamples(self, datadir = None):
		patients = self.allPatients()
		groups   = self.allGroups()

		ret = []
		for patient in patients:
			samples = self.getSamples(by = 'Patient', value = patient, returnAll = True)
			s1 = path.join(datadir, samples[0].Sample) if datadir else samples[0].Sample
			s2 = path.join(datadir, samples[1].Sample) if datadir else samples[1].Sample
			if samples[0].Group != groups[1]: # make sure it's (tumor, normal) pair
				ret.append((s1, s2))
			else:
				ret.append((s2, s1))
		return ret


