
from os import path
from pyppl import Aggr, Channel
from bioprocs.fastx import pFastqTrim, pFastq2Sam, pFastQC, pFastMC, pFastqSETrim, pFastqSE2Sam
from bioprocs.sambam import pBam2Fastq, pSam2Bam, pBam2Counts, pBamRecal, pBam2Gmut, pBamPair2Smut, pBam2Cnv
from bioprocs.common import pFile2Proc, pFiles2Dir
from bioprocs.vcf import pVcf2Maf
from bioprocs.vcfnext import pMafMerge
from bioprocs.utils.sampleinfo import SampleInfo

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
	aFastqSE2Bam
@description:
	Convert fastq single-end files to call-ready bam files.
@depends:
	```
	pFastqSETrim[*] -> pFastqSE2Sam -> pSam2Bam -> pBamRecal[!]
	```
"""
aFastqSE2Bam = Aggr(
	pFastqSETrim,
	pFastqSE2Sam,
	pSam2Bam,
	pBamRecal
)
# delegates
aFastqSE2Bam.delegate('args.ref', 'pFastqSE2Sam, pBamRecal')
# args
aFastqSE2Bam.pSam2Bam.args.markdup = True

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
aFastq2BamQC = Aggr(
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
	aFastqSE2BamQC
@description:
	Convert fastq single-end files to call-ready bam files and do QC for fastq files.
@depends:
	```
	pFastqSETrim[*] -> pFastqSE2Sam -> pSam2Bam -> pBamRecal[!]
	    \
	       pFastQC -> pFastMC[!]
	```
"""
aFastqSE2BamQC = Aggr(
	pFastqSETrim,
	pFastQC,
	pFastMC,
	pFastqSE2Sam,
	pSam2Bam,
	pBamRecal
)
# delegates
aFastqSE2BamQC.delegate('args.ref', 'pFastqSE2Sam, pBamRecal')
# depends
aFastqSE2BamQC.pFastqSE2Sam.depends = aFastqSE2Bam.pFastqSETrim
# args
aFastqSE2BamQC.pSam2Bam.args.markdup = True

"""
@name:
	aBam2MutCnv
@description:
	Call germline, somatic mutations and CNV from call-read bams
@depends:
	pSampleInfo
"""
aBam2MutCnv = Aggr(
	pFile2Proc.copy(id = 'pBamDir'),     # single job
	pFile2Proc.copy(id = 'pSampleInfo'), # single job
	pBam2Gmut,
	pBamPair2Smut,
	pBam2Cnv,
	depends = False
)
# defaults
aBam2MutCnv.pBamDir.runner     = 'local'
aBam2MutCnv.pSampleInfo.runner = 'local'
# delegates
aBam2MutCnv.delegate('args.ref', 'pBam2Gmut, pBamPair2Smut, pBam2Cnv')
# depends
aBam2MutCnv.starts                = aBam2MutCnv.pBamDir,   aBam2MutCnv.pSampleInfo
aBam2MutCnv.ends                  = aBam2MutCnv.pBam2Gmut, aBam2MutCnv.pBamPair2Smut, aBam2MutCnv.pBam2Cnv
aBam2MutCnv.pBam2Gmut.depends     = aBam2MutCnv.pBamDir,   aBam2MutCnv.pSampleInfo
aBam2MutCnv.pBamPair2Smut.depends = aBam2MutCnv.pBamDir,   aBam2MutCnv.pSampleInfo
aBam2MutCnv.pBam2Cnv.depends      = aBam2MutCnv.pBamDir,   aBam2MutCnv.pSampleInfo
# input
aBam2MutCnv.pBam2Gmut.input     = lambda ch1, ch2: SampleInfo(ch2.get()).toChannel(ch1.get())
aBam2MutCnv.pBam2Cnv.input      = lambda ch1, ch2: SampleInfo(ch2.get()).toChannel(ch1.get())
aBam2MutCnv.pBamPair2Smut.input = lambda ch1, ch2: SampleInfo(ch2.get()).toChannel(ch1.get(), paired = True)

"""
@name:
	aBam2Mut
@description:
	Call germline and somatic mutations from call-read bams
@depends:
	pSampleInfo
"""
aBam2Mut = Aggr(
	pFile2Proc.copy(id = 'pBamDir'),     # single job
	pFile2Proc.copy(id = 'pSampleInfo'), # single job
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
aBam2Mut.ends                  = aBam2Mut.pBam2Gmut, aBam2Mut.pBamPair2Smut, aBam2Mut.pBam2Cnv
aBam2Mut.pBam2Gmut.depends     = aBam2Mut.pBamDir,   aBam2Mut.pSampleInfo
aBam2Mut.pBamPair2Smut.depends = aBam2Mut.pBamDir,   aBam2Mut.pSampleInfo
# input
aBam2Mut.pBam2Gmut.input     = lambda ch1, ch2: SampleInfo(ch2.get()).toChannel(ch1.get())
aBam2Mut.pBamPair2Smut.input = lambda ch1, ch2: SampleInfo(ch2.get()).toChannel(ch1.get(), paired = True)


aVcfs2Maf = Aggr(
	pVcf2Maf,
	pFiles2Dir,
	pMafMerge
)
aVcfs2Maf.pFiles2Dir.input     = lambda ch: [ch.flatten()]
