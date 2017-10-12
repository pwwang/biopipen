
from pyppl import Aggr
from bioprocs.fastx import pFastqTrim, pFastq2Sam, pFastQC, pFastMC
from bioprocs.sambam import pBam2Fastq, pSam2Bam, pBam2Counts
"""
@name:
	aEBam2Bam
@description:
	Convert external bam to pair-end fastq and then call-ready bam file.
"""
aEBam2Bam = Aggr(
	pBam2Fastq,
	pFastqTrim,
	pFastq2Sam,
	pSam2Bam
)
# delegates
aEBam2Bam.delegate('args.ref', 'pFastq2Sam')
# args
aEBam2Bam.pSam2Bam.args.markdup = True

"""
@name:
	aEBam2BamQC
@description:
	Convert external bam to pair-end fastq and then call-ready bam file with QC for fastq files.
"""
aEBam2BamQC = Aggr(
	pBam2Fastq,
	pFastqTrim,
	pFastQC,
	pFastMC,
	pFastq2Sam,
	pSam2Bam
)
# depends
aEBam2BamQC.pFastq2Sam.depends = aEBam2BamQC.pFastq2Sam.pFastqTrim
aEBam2BamQC.pFastq2Sam.ends    = aEBam2BamQC.pFastq2Sam.pFastMC, aEBam2BamQC.pFastq2Sam.pSam2Bam
# delegates
aEBam2BamQC.delegate('args.ref', 'pFastq2Sam')
# args
aEBam2BamQC.pSam2Bam.args.markdup = True

