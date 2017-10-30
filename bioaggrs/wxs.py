
from os import path
from pyppl import Aggr, Channel
from bioprocs.fastx import pFastqTrim, pFastq2Sam, pFastQC, pFastMC
from bioprocs.sambam import pBam2Fastq, pSam2Bam, pBam2Counts, pBamRecal, pBam2Gmut, pBamPair2Smut
from bioprocs.common import pFile2Proc
"""
@name:
	aEBam2Bam
@description:
	Convert external bam to pair-end fastq and then call-ready bam file.
@depends:
	```
	pBam2Fastq[*] -> pFastqTrim -> pFastq2Sam -> pSam2Bam -> pBamRecal[!]
	```
"""
aEBam2Bam = Aggr(
	pBam2Fastq,
	pFastqTrim,
	pFastq2Sam,
	pSam2Bam,
	pBamRecal
)
# delegates
aEBam2Bam.delegate('args.ref', 'pFastq2Sam, pBamRecal')
aEBam2Bam.delegate('args.gz',  'pBam2Fastq, pFastqTrim')
# args
aEBam2Bam.pSam2Bam.args.markdup = True

"""
@name:
	aEBam2BamQC
@description:
	Convert external bam to pair-end fastq and then call-ready bam file with QC for fastq files.
@depends:
	```
	pBam2Fastq[*] -> pFastqTrim -> pFastq2Sam -> pSam2Bam -> pBamRecal[!]
	                     \
	                        pFastQC -> pFastMC[!]
	```
"""
aEBam2BamQC = Aggr(
	pBam2Fastq,
	pFastqTrim,
	pFastQC,
	pFastMC,
	pFastq2Sam,
	pSam2Bam,
	pBamRecal
)
# depends
aEBam2BamQC.ends               = aEBam2BamQC.pFastMC, aEBam2BamQC.pSam2Bam
aEBam2BamQC.pFastq2Sam.depends = aEBam2BamQC.pFastqTrim
# delegates
aEBam2BamQC.delegate('args.ref', 'pFastq2Sam, pBamRecal')
# args
aEBam2BamQC.pSam2Bam.args.markdup = True

"""
@name:
	aFastq2Bam
@description:
	Convert fastq pair-end files to call-ready bam files.
@depends:
	```
	pFastqTrim[*] -> pFastq2Sam -> pSam2Bam -> pBamRecal[!]
	```
"""
aFastq2Bam = Aggr(
	pFastqTrim,
	pFastq2Sam,
	pSam2Bam,
	pBamRecal
)
# delegates
aFastq2Bam.delegate('args.ref', 'pFastq2Sam, pBamRecal')
# args
aFastq2Bam.pSam2Bam.args.markdup = True

"""
@name:
	aFastq2BamQC
@description:
	Convert fastq pair-end files to call-ready bam files and do QC for fastq files.
@depends:
	```
	pFastqTrim[*] -> pFastq2Sam -> pSam2Bam -> pBamRecal[!]
	    \
	       pFastQC -> pFastMC[!]
	```
"""
aFastq2Bam = Aggr(
	pFastqTrim,
	pFastQC,
	pFastMC,
	pFastq2Sam,
	pSam2Bam,
	pBamRecal
)
# delegates
aFastq2Bam.delegate('args.ref', 'pFastq2Sam, pBamRecal')
# depends
aFastq2Bam.pFastq2Sam.depends = aFastq2Bam.pFastqTrim
# args
aFastq2Bam.pSam2Bam.args.markdup = True

"""
@name:
	aBam2Mut
@description:
	Call germline and somatic mutations from call-read bams
@depends:
	pSampleInfo
"""
aBam2Mut = Aggr(
	pFile2Proc.copy(newid = 'pBamDir'),
	pFile2Proc.copy(newid = 'pSampleInfo'),
	pBam2Gmut,
	pBamPair2Smut,
	depends = False
)
# defaults
aBam2Mut.pBamDir.runner     = 'local'
aBam2Mut.pSampleInfo.runner = 'local'
# delegates
aBam2Mut.delegate('args.ref', 'pBam2Gmut, pBamPair2Smut')
# depends
aBam2Mut.starts                = aBam2Mut.pBamDir,   aBam2Mut.pSampleInfo
aBam2Mut.ends                  = aBam2Mut.pBam2Gmut, aBam2Mut.pBamPair2Smut
aBam2Mut.pBam2Gmut.depends     = aBam2Mut.pBamDir,   aBam2Mut.pSampleInfo
aBam2Mut.pBamPair2Smut.depends = aBam2Mut.pBamDir,   aBam2Mut.pSampleInfo
# input
aBam2Mut.pBam2Gmut.input     = lambda ch1, ch2: [ 
	path.join(ch1.get(), s + '.bam') for s in     \
	set(Channel.fromFile(ch2.get(), header=True).colAt(0).flatten())
]
aBam2Mut.pBamPair2Smut.input = lambda ch1, ch2: [ \
	(path.join(ch1.get(), c1.get() + '.bam'), path.join(ch1.get(), c2.get() + '.bam')) for ch in \
	[Channel.fromFile(ch2.get(), header=True)] for p in \
	set(ch.Patient.flatten()) for c1, c2 in \
	[ 
		(Channel.fromChannels(ch.colAt(0), ch.Patient, ch.Group).filter(lambda x: x[1] == p and x[2] == 'Tumor'), 
		 Channel.fromChannels(ch.colAt(0), ch.Patient, ch.Group).filter(lambda x: x[1] == p and x[2] == 'Normal'))
	]
]
