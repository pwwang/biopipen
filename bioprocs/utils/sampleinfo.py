"""Sample information handling
A typical sample information file is like
```
Sample  Patient Group   Batch
S1      P1      Tumor   Batch1
S2      P1      Normal  Batch1
...
```
Sample and Group columns are required.
"""

from os import path
from bioprocs.utils.tsvio2 import TsvReader, TsvRecord


class SampleInfoException(Exception):
    """When parsing or extracting fails"""


class SampleInfo:
    """SampleInfo 1.0"""
    # pylint:disable=invalid-name,missing-function-docstring
    NORMAL_GROUP = ['NORMAL', 'HEALTH', 'BEFORE', 'BLOOD', 'PRE', 'CONTROL']
    TUMOR_GROUP = ['TUMOR', 'DISEASE', 'AFTER', 'TISSUE', 'POST', 'TREATMENT']

    def __init__(self, sfile):
        reader = TsvReader(sfile)
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
            raise SampleInfoException(
                'Unexpected column names: %s.' % str(self.colnames))

    def select(self, # pylint: disable=too-many-arguments
               sample=None,
               patient=None,
               group=None,
               batch=None,
               get=None):
        """
        Select records by different standards
        """
        ret = []
        if not isinstance(sample, list):
            sample = [sample]
        if not isinstance(patient, list):
            patient = [patient]
        if not isinstance(group, list):
            group = [group]
        if not isinstance(batch, list):
            batch = [batch]
        for row in self.data:
            if (row[self.samcol] not in sample
                    or row.Patient not in patient
                    or row.Group not in group
                    or row.Batch not in batch):
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

    def toChannel(self, datadir, paired=False, raiseExc=True):
        if not paired:
            samples = self.select(get='Sample')
            return [path.join(datadir, r.Sample) for r in samples]

        data = self.select(get=['Sample', 'Patient'])
        patients = [r.Patient for r in data]
        patients = sorted(set(patients), key=patients.index)
        ret = []
        for patient in patients:
            records = self.select(patient=patient, get=['Sample', 'Group'])
            if len(records) != 2:
                if raiseExc:
                    raise SampleInfoException(
                        'Expect 2 samples for paired channel '
                        'for patient: %s' % patient
                    )
                continue
            record1, record2 = records[0], records[1]
            # make sure record1 is tumor and record2 is normal
            if ((record1.Group.upper() in SampleInfo.TUMOR_GROUP
                 or record1.Group.upper() in SampleInfo.NORMAL_GROUP)
                    and (record2.Group.upper() in SampleInfo.TUMOR_GROUP
                         or record2.Group.upper() in SampleInfo.NORMAL_GROUP)):
                if record1.Group.upper() in SampleInfo.NORMAL_GROUP:
                    record1, record2 = record2, record1
            else:
                common = path.commonprefix(
                    [record1.Sample, record2.Sample])
                if common and record1.Sample < record2.Sample:
                    record1, record2 = record2, record1
            ret.append((path.join(datadir, record1.Sample),
                        path.join(datadir, record2.Sample)))
        return ret

    def dataframe(self, data=None):
        data = data or self.data
        from pandas import DataFrame
        rownames = []
        dataframe = {}
        for row in data:
            for key, val in row.items():
                if key == self.samcol:
                    rownames.append(val)
                elif key in dataframe:
                    dataframe[key].append(val)
                else:
                    dataframe[key] = [val]
        return DataFrame(data=dataframe, index=rownames)


def _unique(array, keep_order=True):
    if not keep_order:
        return list(set(array))
    ret = []
    for element in array:
        if element not in ret:
            ret.append(element)
    return ret

class SampleInfo2:
    """SampleInfo 2.0"""
    def toChannel(self, datadir, paired=False): # pylint: disable=invalid-name
        """Convert to pyppl channel"""
        if not paired:
            samples = self.getSamples()
            return [path.join(datadir, sample) for sample in samples]

        paired_samples = self.getPairedSamples()
        ret = []
        for sample1, sample2 in paired_samples:
            if paired is not True and paired.upper() not in sample1.upper():
                sample1, sample2 = sample2, sample1
            ret.append((path.join(datadir, sample1),
                        path.join(datadir, sample2)))
        return ret

    to_channel = toChannel

    def _read(self, sifile):
        standard_cnames = ["", "Sample", "Patient", "Group", "Batch"]
        reader = TsvReader(sifile)

        self.cnames = reader.cnames
        if not self.cnames:
            raise SampleInfoException(
                'Headers for sample information file is required.')

        if any(cname not in standard_cnames for cname in self.cnames):
            raise SampleInfoException(
                'Headers should be a subset of '
                '{!r}'.format(', '.join(standard_cnames))
            )

        if "" in self.cnames:
            self.cnames[self.cnames.index("")] = "Sample"

        self.mat = reader.dump()

    def __init__(self, sifile, checkPaired=False):
        self.mat = None
        self.cnames = []
        self._read(sifile)
        if checkPaired and "Patient" in self.cnames:
            for patient in self.allPatients():
                if len(self.getSamples(by='Patient', value=patient)) != 2:
                    raise SampleInfoException(
                        'Expect paired comparisons, but Patient '
                        '{!r} has # samples other than 2.'.format(patient)
                    )

    def allPatients(self): # pylint: disable=invalid-name
        """Get all patients"""
        if 'Patient' not in self.cnames:
            return None
        patients = [r.Patient for r in self.mat]
        return _unique(patients)

    all_patients = allPatients

    def all_batches(self):
        """Get all batches"""
        if 'Batch' not in self.cnames:
            return None
        batches = [r.Batch for r in self.mat]
        return _unique(batches)


    def allGroups(self): # pylint: disable=invalid-name
        """Get all groups"""
        allgroups = [r.Group for r in self.mat]
        ret = []
        for group in allgroups:
            if group not in ret:
                ret.append(group)
        return ret

    all_groups = allGroups

    def allSamples(self, unique=False, datadir=None): # pylint: disable=invalid-name
        """Get all samples"""
        samples = self.getSamples()
        if unique:
            ret = []
            for sample in samples:
                if sample not in ret:
                    ret.append(sample
                               if datadir is None
                               else path.join(datadir, sample))
            return ret
        return [path.join(datadir, sample)
                for sample in samples] if datadir else samples

    all_samples = allSamples

    def getSamples(self,  # pylint: disable=invalid-name
                   by=None,
                   value=None,
                   return_all=False):
        """Get the samples by certain criteria"""
        if by and by not in self.cnames:
            raise SampleInfoException(
                '{!r} is not a valid column name.'.format(by))
        if not by:
            return [r if return_all else r.Sample for r in self.mat]
        return [r if return_all else r.Sample
                for r in self.mat if r[by] == value]

    get_samples = getSamples

    def sampleInfo(self, sample, info=None): # pylint: disable=invalid-name
        """Get the given information of the given sample"""
        if not info:
            return [r for r in self.mat if r.Sample == sample]
        if info not in self.cnames:
            raise SampleInfoException(
                '{!r} is not a valid column name.'.format(info))
        return [r[info] for r in self.mat if r.Sample == sample]

    sample_info = sampleInfo

    def select(self, samples):
        """Select the records with given samples"""
        self.mat = [r for r in self.mat if r.Sample in samples]

    def isPaired(self): # pylint: disable=invalid-name
        """Tell the sample information is paired or not"""
        allpatients = self.allPatients()
        if not allpatients:
            return False
        for patient in allpatients:
            if len(self.getSamples(by='Patient', value=patient)) != 2:
                return False
        return True

    is_paired = isPaired

    def getPairedSample(self, sample): # pylint: disable=invalid-name
        """Get the paired sample of the given sample"""
        patient = self.sampleInfo(sample, info='Patient')[0]
        samples = self.getSamples(by='Patient', value=patient)
        return samples[1 - samples.index(sample)]

    get_paired_sample = getPairedSample

    def getPairedSamples(self, datadir=None): # pylint: disable=invalid-name
        """Get sample pairs"""
        patients = self.allPatients()
        groups = self.allGroups()

        ret = []
        for patient in patients:
            samples = self.getSamples(
                by='Patient', value=patient, return_all=True)
            sam1 = path.join(
                datadir, samples[0].Sample) if datadir else samples[0].Sample
            sam2 = path.join(
                datadir, samples[1].Sample) if datadir else samples[1].Sample
            # make sure it's (tumor, normal) pair
            if samples[0].Group != groups[1]:
                ret.append((sam1, sam2))
            else:
                ret.append((sam2, sam1))
        return ret

    get_paired_samples = getPairedSamples
