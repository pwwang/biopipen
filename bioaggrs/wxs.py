
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
	aPrepareBam
@description:
	Prepare call-ready bam files from external bams or fastq files
@procs:
	`pBam2Fastq`: Convert external bam files into fastq files
	`pFastqTrim`: Trim fastq reads.
	`pFastQC`   : Run fastqc
	`pFastMC`   : Run multiqc
	`pFastq2Sam`: Convert fastq files into sam/bam files.
	`pSam2Bam`  : Convert sam files into bam files and do some cleaning (sort, markdup)
	`pBamRecal` : Recalibrate bam file.
@modules
	`ebam` : The aggregation starts with external bam files (off by default)
		- `pBam2Fastq[*]` -> `pFastqTrim` -> ...
	`fastq`: The aggregation starts with external fastq files (on by default, mutually exclusive with `ebam`)
		- `pFastqTrim[*]` -> ...
	`qc`   : Do qc for the fastq files (off by default)
		- `pFastqTrim` -> `pFastQC, pFastMC`
"""
aPrepareBam = Aggr(
	pBam2Fastq,
	pFastqTrim,
	pFastQC,
	pFastMC,
	pFastq2Sam,
	pSam2Bam,
	pBamRecal,
	depends = False
)
# starts and ends
aPrepareBam.ends = 'pBamRecal'
# depends
aPrepareBam['pFastq2Sam',].depends = ['pFastqTrim']
aPrepareBam['pSam2Bam',].depends   = ['pFastq2Sam']
aPrepareBam['pBamRecal',].depends  = ['pSam2Bam']
# delegates
aPrepareBam.delegate('args.ref', 'pFastq2Sam, pBamRecal')
# args
aPrepareBam.pSam2Bam.args.markdup = True
# modules
aPrepareBam.module('ebam', starts = 'pBam2Fastq', depends = {'pFastqTrim': 'pBam2Fastq'})
aPrepareBam.module('fastq', starts = 'pFastqTrim')
aPrepareBam.module('qc', ends = 'pFastQC, pFastMC', depends = {'pFastQC, pFastMC': ['pFastqTrim'] * 2})
# inital modules
aPrepareBam.on('fastq')


"""
@name:
	aBam2SCNV
@description:
	Call CNV with paired control samples
@procs:
	`pBamDir`       : The directory with the bam files
	`pSampleInfo`   : The sample information
	`pCNVkitPrepare`: Prepare the region/coverage file for CNV calling
	`pCNVkitCov`    : Get the coverage file for each bam file.
	`pCNNDir`       : Put all coverage file into a directory.
	`pCNVkitRef`    : Aggregate the reference cnn file from controls.
	`pCNVkitFix`    : Fix the coverage for each sample using the reference cnn file.
	`pCNVkitSeg`    : Segment the regions
	`pCNVkitCall`   : Call CNVs.
	`pCNVkitScatter`: Draw scatter plot for CNVs.
	`pCNVkitDiagram`: Draw diagram plot for CNVs.
	`pCNVkitHeatmap`: Draw heatmap for CNVs.
	`pCNVkitReport` : Generate a report for the results.
	`pCNVkit2Vcf`   : Export the results to VCF file.
@modules
	`plots`: Draw plots
		- ... -> `pCNVkitFix, pCNVkitSeg` -> `pCNVkitScatter[!], pCNVkitDiagram[!], pCNVkitHeatmap[!], pCNVkitReport[!]`
"""
aBam2SCNV = Aggr(
	pFile2Proc.copy(id = 'pBamDir'),     # single job
	pFile2Proc.copy(id = 'pSampleInfo'), # single job
	pCNVkitPrepare,   # prepare files 
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
aBam2SCNV.pBamDir.runner     = 'local'
aBam2SCNV.pSampleInfo.runner = 'local'
# delegates
aBam2SCNV.delegate('args.ref', 'pCNVkitPrepare, pCNVkitRef')
aBam2SCNV.delegate('args.cnvkit', 'pCNVkitRef, pCNVkitFix, pCNVkitSeg, pCNVkitCall, pCNVkitScatter, pCNVkitDiagram, pCNVkitHeatmap, pCNVkitReport, pCNVkit2Vcf')
aBam2SCNV.delegate('args.nthread', 'pCNVkitCov, pCNVkitSeg')
# depends
aBam2SCNV.starts                 = aBam2SCNV.pBamDir,        aBam2SCNV.pSampleInfo
aBam2SCNV.ends                   = aBam2SCNV.pCNVkit2Vcf
aBam2SCNV.pCNVkitPrepare.depends = aBam2SCNV.pBamDir,        aBam2SCNV.pSampleInfo
aBam2SCNV.pCNVkitCov.depends     = aBam2SCNV.pBamDir,        aBam2SCNV.pSampleInfo,    aBam2SCNV.pCNVkitPrepare
aBam2SCNV.pCNNDir.depends        = aBam2SCNV.pCNVkitCov
aBam2SCNV.pCNVkitRef.depends     = aBam2SCNV.pCNNDir,        aBam2SCNV.pSampleInfo
aBam2SCNV.pCNVkitFix.depends     = aBam2SCNV.pCNNDir,        aBam2SCNV.pSampleInfo,    aBam2SCNV.pCNVkitRef
aBam2SCNV.pCNVkitSeg.depends     = aBam2SCNV.pCNVkitFix
aBam2SCNV.pCNVkitCall.depends    = aBam2SCNV.pCNVkitSeg
aBam2SCNV.pCNVkit2Vcf.depends    = aBam2SCNV.pCNVkitCall
# input
aBam2SCNV.pCNVkitPrepare.input = lambda ch_bamdir, ch_saminfo: [
	Channel.create(
		SampleInfo(ch_saminfo.get()).toChannel(ch_bamdir.get())
	).unique().flatten()]
aBam2SCNV.pCNVkitCov.input     = lambda ch_bamdir, ch_saminfo, ch_target: \
	Channel.create(
		SampleInfo(ch_saminfo.get()).toChannel(ch_bamdir.get())
	).unique().cbind(ch_target)
aBam2SCNV.pCNNDir.input        = lambda ch: [ch.flatten()]
aBam2SCNV.pCNVkitRef.input     = lambda ch_covs, ch_saminfo: [
	Channel.create(
		SampleInfo(ch_saminfo.get()).toChannel(ch_covs.get(), paired = True, raiseExc = False)
	).colAt(1).unique().map(
		lambda x: (x[0].rpartition('.')[0] + '.target.cnn', x[0].rpartition('.')[0] + '.antitarget.cnn')
	).flatten()]
aBam2SCNV.pCNVkitFix.input     = lambda ch_covs, ch_saminfo, ch_ref: \
	Channel.create(
		SampleInfo(ch_saminfo.get()).toChannel(ch_covs.get(), paired = True)
	).colAt(0).unique().map(
		lambda x: (x[0].rpartition('.')[0] + '.target.cnn', x[0].rpartition('.')[0] + '.antitarget.cnn')
	).cbind(ch_ref)
aBam2SCNV.pCNVkitHeatmap.input = lambda ch: [ch.flatten()]
# module
aBam2SCNV.module('plots', ends = 'pCNVkitScatter, pCNVkitDiagram, pCNVkitHeatmap, pCNVkitReport', depends = {
	'pCNVkitScatter': 'pCNVkitFix, pCNVkitSeg',
	'pCNVkitDiagram': 'pCNVkitFix, pCNVkitSeg',
	'pCNVkitHeatmap': 'pCNVkitSeg',
	'pCNVkitReport' : 'pCNVkitFix, pCNVkitSeg',
})

"""
@name:
	aBam2GCNV
@description:
	Call CNV with a flat reference coverage file.
@procs:
	`pBamDir`       : The directory with the bam files
	`pSampleInfo`   : The sample information
	`pCNVkitPrepare`: Prepare the region/coverage file for CNV calling
	`pCNVkitCov`    : Get the coverage file for each bam file.
	`pCNNDir`       : Put all coverage file into a directory.
	`pCNVkitFlatRef`: Generate a flat reference coverage file.
	`pCNVkitFix`    : Fix the coverage for each sample using the reference cnn file.
	`pCNVkitSeg`    : Segment the regions
	`pCNVkitCall`   : Call CNVs.
	`pCNVkitScatter`: Draw scatter plot for CNVs.
	`pCNVkitDiagram`: Draw diagram plot for CNVs.
	`pCNVkitHeatmap`: Draw heatmap for CNVs.
	`pCNVkitReport` : Generate a report for the results.
	`pCNVkit2Vcf`   : Export the results to VCF file.
@modules
	`plots`: Draw plots
		- ... -> `pCNVkitFix, pCNVkitSeg` -> `pCNVkitScatter[!], pCNVkitDiagram[!], pCNVkitHeatmap[!], pCNVkitReport[!]`
"""
aBam2GCNV = Aggr(
	pFile2Proc.copy(id = 'pBamDir'),     # single job
	pFile2Proc.copy(id = 'pSampleInfo'), # single job
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
aBam2GCNV.pBamDir.runner     = 'local'
aBam2GCNV.pSampleInfo.runner = 'local'
# delegates
aBam2GCNV.delegate('args.ref', 'pCNVkitPrepare, pCNVkitFlatRef')
aBam2GCNV.delegate('args.nthread', 'pCNVkitCov, pCNVkitSeg')
aBam2GCNV.delegate('args.cnvkit', 'pCNVkitFlatRef, pCNVkitFix, pCNVkitSeg, pCNVkitCall, pCNVkitScatter, pCNVkitDiagram, pCNVkitHeatmap, pCNVkitReport, pCNVkit2Vcf')
# depends
aBam2GCNV.starts                 = aBam2GCNV.pBamDir,        aBam2GCNV.pSampleInfo
aBam2GCNV.ends                   = aBam2GCNV.pCNVkit2Vcf
aBam2GCNV.pCNVkitPrepare.depends = aBam2GCNV.pBamDir,        aBam2GCNV.pSampleInfo
aBam2GCNV.pCNVkitCov.depends     = aBam2GCNV.pBamDir,        aBam2GCNV.pSampleInfo,    aBam2GCNV.pCNVkitPrepare
aBam2GCNV.pCNNDir.depends        = aBam2GCNV.pCNVkitCov
aBam2GCNV.pCNVkitFlatRef.depends = aBam2GCNV.pCNVkitPrepare
aBam2GCNV.pCNVkitFix.depends     = aBam2GCNV.pCNNDir,        aBam2GCNV.pSampleInfo,    aBam2GCNV.pCNVkitFlatRef
aBam2GCNV.pCNVkitSeg.depends     = aBam2GCNV.pCNVkitFix
aBam2GCNV.pCNVkitCall.depends    = aBam2GCNV.pCNVkitSeg
aBam2GCNV.pCNVkit2Vcf.depends    = aBam2GCNV.pCNVkitCall
# input
#aBam2GCNV.pCNVkitAccess.input  = lambda ch: path.basename(ch.get()).split('.')[:1]
aBam2GCNV.pCNVkitPrepare.input = lambda ch_bamdir, ch_saminfo: [
	Channel.create(
		SampleInfo(ch_saminfo.get()).toChannel(ch_bamdir.get())
	).unique().flatten()]
aBam2GCNV.pCNVkitCov.input     = lambda ch_bamdir, ch_saminfo, ch_target: \
	Channel.create(
		SampleInfo(ch_saminfo.get()).toChannel(ch_bamdir.get())
	).unique().cbind(ch_target)
aBam2GCNV.pCNNDir.input        = lambda ch: [ch.flatten()]
aBam2GCNV.pCNVkitFix.input     = lambda ch_covs, ch_saminfo, ch_ref: \
	Channel.create(
		SampleInfo(ch_saminfo.get()).toChannel(ch_covs.get())
	).colAt(0).unique().map(
		lambda x: (x[0].rpartition('.')[0] + '.target.cnn', x[0].rpartition('.')[0] + '.antitarget.cnn')
	).cbind(ch_ref)
aBam2GCNV.pCNVkitHeatmap.input = lambda ch: [ch.flatten()]
# module
aBam2GCNV.module('plots', ends = 'pCNVkitScatter, pCNVkitDiagram, pCNVkitHeatmap, pCNVkitReport', depends = {
	'pCNVkitScatter': 'pCNVkitFix, pCNVkitSeg',
	'pCNVkitDiagram': 'pCNVkitFix, pCNVkitSeg',
	'pCNVkitHeatmap': 'pCNVkitSeg',
	'pCNVkitReport' : 'pCNVkitFix, pCNVkitSeg',
})

"""
@name:
	aVcfs2Maf
@description:
	Convert vcf files to maf files and merge them.
@procs:
	`pVcf2Maf`  : Convert single vcf file to maf file
	`pFiles2Dir`: Put the maf files into a directory
	`pMafMerge` : Merge the maf files.
"""
aVcfs2Maf = Aggr(
	pVcf2Maf,
	pFiles2Dir,
	pMafMerge
)
aVcfs2Maf.pFiles2Dir.input = lambda ch: [ch.flatten()]
