
from os import path
from pyppl import Aggr, Channel
from bioprocs.fastx import pFastqTrim, pFastq2Sam, pFastQC, pFastMC, pFastqSETrim, pFastqSE2Sam
from bioprocs.sambam import pBam2Fastq, pSam2Bam, pBam2Counts, pBamRecal, pBam2Gmut, pBamPair2Smut
from bioprocs.common import pFile2Proc, pFiles2Dir
from bioprocs.vcf import pVcf2Maf
from bioprocs.vcfnext import pMafMerge
from bioprocs.cnvkit import pCNVkit2Vcf, pCNVkitCall, pCNVkitRef, pCNVkitFlatRef, pCNVkitFix, pCNVkitCov, pCNVkitHeatmap, pCNVkitDiagram, pCNVkitReport, pCNVkitScatter, pCNVkitSeg, pCNVkitPrepare
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

aBam2SCnv = Aggr(
	pFile2Proc.copy(id = 'pBamDir'),     # single job
	pFile2Proc.copy(id = 'pSampleInfo'), # single job
	#pCNVkitAccess,    # single job        , in: index/tag  , out: access file
	#pCNVkitTarget,    # single job        , in: access file, out: target file 
	pCNVkitPrepare,    # prepare files 
	pCNVkitCov,       # each sample       , in: bam file (target required), out: cnn file
	pFiles2Dir.copy(id = 'pCNNDir'),  # put all cnn files in a directory
	pCNVkitRef,       # all normal samples, in: all normal cnn files      , out: refcnn
	pCNVkitFix,       # each tumor sample , in: tumor cnn, refcnn         , out: cnr
	pCNVkitSeg,       # each tumor sample , in: cnr                       , out: cns
	pCNVkitCall,      # each tumor sample , in: cns                       , out: callcns
	pCNVkitScatter,   # each tumor sample , in: cnr, cns                  , out: plots
	pCNVkitDiagram,   # each tumor sample , in: cnr, cns                  , out: plots
	pCNVkitHeatmap,   # all tumor samples , in: cnfiles                   , out: plots
	pCNVkitReport,    # each tumor sample , in: cnr, cns                  , out: reports
	pCNVkit2Vcf,      # each tumor sample , in: callcns                   , out: vcf,
	depends = False
)
# defaults
aBam2SCnv.pBamDir.runner     = 'local'
aBam2SCnv.pSampleInfo.runner = 'local'
# delegates
aBam2SCnv.delegate('args.ref', 'pCNVkitPrepare, pCNVkitRef')
aBam2SCnv.delegate('args.cnvkit')
aBam2SCnv.delegate('args.nthread', 'pCNVkitCov, pCNVkitSeg')
# depends
aBam2SCnv.starts                 = aBam2SCnv.pBamDir,        aBam2SCnv.pSampleInfo
aBam2SCnv.ends                   = aBam2SCnv.pCNVkitScatter, aBam2SCnv.pCNVkitDiagram, aBam2SCnv.pCNVkitHeatmap, aBam2SCnv.pCNVkitReport, aBam2SCnv.pCNVkit2Vcf
#aBam2SCnv.pCNVkitAccess.depends  = aBam2SCnv.pSampleInfo
#aBam2SCnv.pCNVkitTarget.depends  = aBam2SCnv.pCNVkitAccess
aBam2SCnv.pCNVkitPrepare.depends = aBam2SCnv.pBamDir,        aBam2SCnv.pSampleInfo
aBam2SCnv.pCNVkitCov.depends     = aBam2SCnv.pBamDir,        aBam2SCnv.pSampleInfo,    aBam2SCnv.pCNVkitPrepare
aBam2SCnv.pCNNDir.depends        = aBam2SCnv.pCNVkitCov
aBam2SCnv.pCNVkitRef.depends     = aBam2SCnv.pCNNDir,        aBam2SCnv.pSampleInfo
aBam2SCnv.pCNVkitFix.depends     = aBam2SCnv.pCNNDir,        aBam2SCnv.pSampleInfo,    aBam2SCnv.pCNVkitRef
aBam2SCnv.pCNVkitSeg.depends     = aBam2SCnv.pCNVkitFix
aBam2SCnv.pCNVkitCall.depends    = aBam2SCnv.pCNVkitSeg
aBam2SCnv.pCNVkitScatter.depends = aBam2SCnv.pCNVkitFix, aBam2SCnv.pCNVkitSeg
aBam2SCnv.pCNVkitDiagram.depends = aBam2SCnv.pCNVkitFix, aBam2SCnv.pCNVkitSeg
aBam2SCnv.pCNVkitHeatmap.depends = aBam2SCnv.pCNVkitSeg
aBam2SCnv.pCNVkitReport.depends  = aBam2SCnv.pCNVkitFix, aBam2SCnv.pCNVkitSeg
aBam2SCnv.pCNVkit2Vcf.depends    = aBam2SCnv.pCNVkitCall
# input
#aBam2SCnv.pCNVkitAccess.input  = lambda ch: path.basename(ch.get()).split('.')[:1]
aBam2SCnv.pCNVkitPrepare.input = lambda ch_bamdir, ch_saminfo: [Channel.create(SampleInfo(ch_saminfo.get()).toChannel(ch_bamdir.get())).unique().flatten()]
aBam2SCnv.pCNVkitCov.input     = lambda ch_bamdir, ch_saminfo, ch_target: Channel.create(SampleInfo(ch_saminfo.get()).toChannel(ch_bamdir.get())).unique().cbind(ch_target)
aBam2SCnv.pCNNDir.input        = lambda ch: [ch.flatten()]
aBam2SCnv.pCNVkitRef.input     = lambda ch_covs, ch_saminfo: [Channel.create(SampleInfo(ch_saminfo.get()).toChannel(ch_covs.get(), paired = True, raiseExc = False)).colAt(1).unique().map(lambda x: (x[0].rpartition('.')[0] + '.target.cnn', x[0].rpartition('.')[0] + '.antitarget.cnn')).flatten()]
aBam2SCnv.pCNVkitFix.input     = lambda ch_covs, ch_saminfo, ch_ref: Channel.create(SampleInfo(ch_saminfo.get()).toChannel(ch_covs.get(), paired = True)).colAt(0).unique().map(lambda x: (x[0].rpartition('.')[0] + '.target.cnn', x[0].rpartition('.')[0] + '.antitarget.cnn')).cbind(ch_ref)
aBam2SCnv.pCNVkitHeatmap.input = lambda ch: [ch.flatten()]

aBam2GCnv = Aggr(
	pFile2Proc.copy(id = 'pBamDir'),     # single job
	pFile2Proc.copy(id = 'pSampleInfo'), # single job
	#pCNVkitAccess,    # single job        , in: index/tag  , out: access file
	#pCNVkitTarget,    # single job        , in: access file, out: target file 
	pCNVkitPrepare,    # prepare files 
	pCNVkitCov,       # each sample       , in: bam file (target required), out: cnn file
	pFiles2Dir.copy(id = 'pCNNDir'),  # put all cnn files in a directory
	pCNVkitFlatRef,   # flat reference    , in: target file               , out: refcnn
	pCNVkitFix,       # each tumor sample , in: tumor cnn, refcnn         , out: cnr
	pCNVkitSeg,       # each tumor sample , in: cnr                       , out: cns
	pCNVkitCall,      # each tumor sample , in: cns                       , out: callcns
	pCNVkitScatter,   # each tumor sample , in: cnr, cns                  , out: plots
	pCNVkitDiagram,   # each tumor sample , in: cnr, cns                  , out: plots
	pCNVkitHeatmap,   # all tumor samples , in: cnfiles                   , out: plots
	pCNVkitReport,    # each tumor sample , in: cnr, cns                  , out: reports
	pCNVkit2Vcf,      # each tumor sample , in: callcns                   , out: vcf,
	depends = False
)
# defaults
aBam2GCnv.pBamDir.runner     = 'local'
aBam2GCnv.pSampleInfo.runner = 'local'
# delegates
aBam2GCnv.delegate('args.ref', 'pCNVkitPrepare, pCNVkitFlatRef')
aBam2GCnv.delegate('args.cnvkit')
aBam2GCnv.delegate('args.nthread', 'pCNVkitCov, pCNVkitSeg')
# depends
aBam2GCnv.starts                 = aBam2GCnv.pBamDir,        aBam2GCnv.pSampleInfo
aBam2GCnv.ends                   = aBam2GCnv.pCNVkitScatter, aBam2GCnv.pCNVkitDiagram, aBam2GCnv.pCNVkitHeatmap, aBam2GCnv.pCNVkitReport, aBam2GCnv.pCNVkit2Vcf
#aBam2GCnv.pCNVkitAccess.depends  = aBam2GCnv.pSampleInfo
#aBam2GCnv.pCNVkitTarget.depends  = aBam2GCnv.pCNVkitAccess
aBam2GCnv.pCNVkitPrepare.depends = aBam2GCnv.pBamDir,        aBam2GCnv.pSampleInfo
aBam2GCnv.pCNVkitCov.depends     = aBam2GCnv.pBamDir,        aBam2GCnv.pSampleInfo,    aBam2GCnv.pCNVkitPrepare
aBam2GCnv.pCNNDir.depends        = aBam2GCnv.pCNVkitCov
aBam2GCnv.pCNVkitFlatRef.depends = aBam2GCnv.pCNVkitPrepare
aBam2GCnv.pCNVkitFix.depends     = aBam2GCnv.pCNNDir,        aBam2GCnv.pSampleInfo,    aBam2GCnv.pCNVkitFlatRef
aBam2GCnv.pCNVkitSeg.depends     = aBam2GCnv.pCNVkitFix
aBam2GCnv.pCNVkitCall.depends    = aBam2GCnv.pCNVkitSeg
aBam2GCnv.pCNVkitScatter.depends = aBam2GCnv.pCNVkitFix, aBam2GCnv.pCNVkitSeg
aBam2GCnv.pCNVkitDiagram.depends = aBam2GCnv.pCNVkitFix, aBam2GCnv.pCNVkitSeg
aBam2GCnv.pCNVkitHeatmap.depends = aBam2GCnv.pCNVkitSeg
aBam2GCnv.pCNVkitReport.depends  = aBam2GCnv.pCNVkitFix, aBam2GCnv.pCNVkitSeg
aBam2GCnv.pCNVkit2Vcf.depends    = aBam2GCnv.pCNVkitCall
# input
#aBam2GCnv.pCNVkitAccess.input  = lambda ch: path.basename(ch.get()).split('.')[:1]
aBam2GCnv.pCNVkitPrepare.input = lambda ch_bamdir, ch_saminfo: [Channel.create(SampleInfo(ch_saminfo.get()).toChannel(ch_bamdir.get())).unique().flatten()]
aBam2GCnv.pCNVkitCov.input     = lambda ch_bamdir, ch_saminfo, ch_target: Channel.create(SampleInfo(ch_saminfo.get()).toChannel(ch_bamdir.get())).unique().cbind(ch_target)
aBam2GCnv.pCNNDir.input        = lambda ch: [ch.flatten()]
aBam2GCnv.pCNVkitFix.input     = lambda ch_covs,   ch_saminfo, ch_ref: Channel.create(SampleInfo(ch_saminfo.get()).toChannel(ch_covs.get())).colAt(0).unique().map(lambda x: (x[0].rpartition('.')[0] + '.target.cnn', x[0].rpartition('.')[0] + '.antitarget.cnn')).cbind(ch_ref)
aBam2GCnv.pCNVkitHeatmap.input = lambda ch: [ch.flatten()]


aVcfs2Maf = Aggr(
	pVcf2Maf,
	pFiles2Dir,
	pMafMerge
)
aVcfs2Maf.pFiles2Dir.input     = lambda ch: [ch.flatten()]
