
from pyppl import Aggr
from bioprocs.fastx import pFastqTrim, pFastq2Sam, pFastQC, pFastMC
from bioprocs.sambam import pBam2Fastq, pSam2Bam, pBam2Counts, pBamRecal
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
	pFastqTrim -> pFastq2Sam -> pSam2Bam -> pBamRecal[!]
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

