import sys, os
sys.path.insert(0, "/home/m161047/tools/pyppl")     # remove this if you have pyppl installed
sys.path.insert(0, "/home/m161047/tools/bioprocs")  # remove this if you have bioprocs installed

from pyppl import pyppl, proc, channel
import unittest

from bioprocs.web import pDownloadGet
from bioprocs.wxsprep import pTrimmomaticPE, pTrimmomaticSE, pAlignPEByBWA, pAlignSEByBWA, pAlignPEByNGM, pAlignSEByNGM
from bioprocs.wxscall import pCNVnator
from bioprocs.wxsstat import pVcf2List, pCallRate
from bioprocs.cnvkit import pCNVkitAccess, pCNVkitTarget, pCNVkitCov, pCNVkitRef, pCNVkitFix, pCNVkitSeg, pCNVkitPlot, pCNVkitRpt, pCNVkit2Vcf
from bioprocs.picard import pSortSam, pMarkDuplicates, pIndexBam
from bioprocs.common import pFiles2Dir, pCbindList
from bioprocs.gatk import pRealignerTargetCreator, pIndelRealigner, pBaseRecalibrator, pHaplotypeCaller, pPrintReads, pSelectVariants, pVariantFiltration, pMuTect2, pMuTect2Interval


class TestWXSRelated (unittest.TestCase):

	def testGerm (self):
		pDownloadGet.input = {pDownloadGet.input: [
			# reference genome
			"https://raw.githubusercontent.com/BenLangmead/bowtie2/master/example/reference/lambda_virus.fa",
			# reads simulator
			"https://raw.githubusercontent.com/BenLangmead/bowtie2/master/example/reads/simulate.pl",
			# trimmometic adapter
			#"https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa"
		]}
		pDownloadGet.desc = 'Download the reference genome and the script to generate the reads'
		
		pChangeChrom = proc ()
		pChangeChrom.depends = pDownloadGet
		pChangeChrom.input = {"fafile:file": lambda ch: ch[0]}
		pChangeChrom.output = "outfile:file:{{ fafile | bn }}"
		pChangeChrom.script = """
		echo ">chrPhage" > "{{outfile}}"
		tail -n +2 "{{fafile}}" >> "{{outfile}}"
		"""
		
		pGeneratePEReads                  = proc ()
		pGeneratePEReads.depends          = [pDownloadGet, pChangeChrom]
		pGeneratePEReads.input            = {"fafile:file, script:file": lambda ch1, ch2: ch2.cbind(ch1[1])}
		pGeneratePEReads.output           = "read1:file:testPE_1.fq, read2:file:testPE_2.fq"
		pGeneratePEReads.script           = """
		> {{read1 | quote}}
		> {{read2 | quote}}
		perl {{script | quote}} --fasta {{fafile | quote}} --prefix {{read1 | prefix | [:-2]}}
		"""
		
		
		# single end
		pGenerateSEReads                  = proc ()
		pGenerateSEReads.depends          = [pDownloadGet, pChangeChrom]
		pGenerateSEReads.input            = {"fafile:file, script:file": lambda ch1, ch2: ch2.cbind(ch1[1])}
		pGenerateSEReads.output           = "read:file:testSE.fq"
		pGenerateSEReads.script           = """
		> {{read | quote}}
		perl {{script | quote}} --fasta {{fafile | quote}} --unpaired --prefix {{read | prefix}}
		"""

		pTrimmomaticPE.depends            = pGeneratePEReads
		pTrimmomaticPE.callfront          = lambda p: p.args.update({"params": "LEADING:3 TRAILING:3"})
		
		pTrimmomaticSE.depends            = pGenerateSEReads
		pTrimmomaticSE.callfront          = lambda p: p.args.update({"params": "LEADING:3 TRAILING:3"})
		
		pAlignPEByBWA.depends             = [pTrimmomaticPE, pChangeChrom]
		
		pAlignSEByBWA.depends             = [pTrimmomaticSE, pChangeChrom]
		
		pAlignPEByNGM.depends             = [pTrimmomaticPE, pChangeChrom]
		
		pAlignSEByNGM.depends             = [pTrimmomaticSE, pChangeChrom]

		pSortSam.depends                  = [pAlignPEByBWA, pAlignSEByBWA]
		pSortSam.input                    = {pSortSam.input: lambda ch1, ch2: ch1 + ch2}

		pMarkDuplicates.depends           = pSortSam
		
		#pIndexBam.depends                 = pMarkDuplicates

		#pRealignerTargetCreator.depends   = [pIndexBam, pDownloadGet]
		pRealignerTargetCreator.depends   = [pMarkDuplicates, pChangeChrom]
		pRealignerTargetCreator.args['params'] = ' -Xms2g -Xmx8g '
		
		pIndelRealigner.depends           = [pMarkDuplicates, pRealignerTargetCreator, pChangeChrom]
		
		pKnownSites                       = proc ()
		pKnownSites.output                = "outfile:file:knownSites.vcf"
		pKnownSites.lang                  = "python"
		pKnownSites.script                = """
		with open ({{outfile | quote}}, "w") as f:
			f.write('''##fileformat=VCFv4.2
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description="Consistency of the site with at most two segregating haplotypes">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">
##contig=<ID=gi|9626243|ref|NC_001416.1|,length=48502>
##reference=lambda_virus.fa
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
chrPhage	1	.	G	A	406.77	.	AC=2;AF=1.00;AN=2;DP=15;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=27.12;SOR=1.112	GT:AD:DP:GQ:PL	1/1:0,15:15:44:435,44,0
''')
		"""
		
		pBaseRecalibrator.depends         = [pIndelRealigner, pChangeChrom, pKnownSites]
		pBaseRecalibrator.input           = {pBaseRecalibrator.input: lambda ch1, ch2, ch3: ch1.cbind(ch2)}
		pBaseRecalibrator.callfront       = lambda p: p.args.update({'knownSites': pKnownSites.channel[0][0], 'params': '--maximum_cycle_value 10000'})
		
		pPrintReads.depends               = [pIndelRealigner, pBaseRecalibrator, pChangeChrom]
		
		pHaplotypeCaller.depends          = [pPrintReads, pChangeChrom]
		
		pSelectVariants.depends           = [pHaplotypeCaller, pChangeChrom]
		pSelectVariants.args['params']    = '-selectType SNP'
		
		pVariantFiltration.depends        = [pSelectVariants, pChangeChrom]
		
		pMuTect2.depends                  = [pPrintReads, pChangeChrom]
		pMuTect2.input                    = {pMuTect2.input: lambda ch1, ch2: ch1.unfold(2).cbind(ch2)}
		pMuTect2.args['params']           = " -Xms2g -Xmx8g "
		
		pMuTect2Interval.depends          = [pPrintReads, pChangeChrom]
		pMuTect2Interval.input            = {pMuTect2Interval.input: lambda ch1, ch2: ch1.unfold(2).cbind(ch2)}
		pMuTect2Interval.args['params']   = " -Xms2g -Xmx8g "
		
		pCNVnator.depends                 = [pPrintReads, pChangeChrom]
		pCNVnator.input                   = {pCNVnator.input: lambda ch1, ch2: ch1}
		pCNVnator.callfront               = lambda p: p.args.update({"binsize":10, "chrdir": os.path.dirname(pChangeChrom.channel[0][0])})
		
		pCNVkitAccess.depends             = pChangeChrom
		pCNVkitAccess.args['params']      = '-s 50'
		
		pGetCNVKitGenes                   = proc ()
		pGetCNVKitGenes.depends           = pCNVkitAccess
		pGetCNVKitGenes.input             = "infile:file"
		pGetCNVKitGenes.output            = "outfile:file:{{infile | fn}}.genes"
		pGetCNVKitGenes.script            = """
		bedtools random -g <(cut -f1,3 "{{infile}}") -n 100 > "{{outfile}}"
		"""
		
		pCNVkitTarget.depends             = [pCNVkitAccess, pGetCNVKitGenes]
		pCNVkitCov.depends                = [pSortSam, pCNVkitTarget]
		pCNVkitCov.input                  = {pCNVkitCov.input: lambda ch1, ch2: ch1}
		pCNVkitCov.callfront              = lambda p: p.args.update({"tgfile": pCNVkitTarget.channel[0][0]})
		
		pFiles2DirCNR                     = pFiles2Dir.copy('cnr')
		pFiles2DirCNS                     = pFiles2Dir.copy('cns')
		pFiles2DirVCF                     = pFiles2Dir.copy('vcf')
		pFiles2Dir.depends                = pCNVkitCov
		pFiles2Dir.input                  = {pFiles2Dir.input: lambda ch: [ch.toList()]}
		
		pCNVkitRef.depends                = pFiles2Dir
		
		pCNVkitFix.depends                = [pCNVkitCov, pCNVkitRef]
		pCNVkitFix.input                  = {pCNVkitFix.input: lambda ch1, ch2: ch1.cbind(ch2)}
		
		pCNVkitSeg.depends                = pCNVkitFix
		
		pFiles2DirCNR.depends             = pCNVkitFix
		pFiles2DirCNR.input               = {pFiles2DirCNR.input: lambda ch: [ch.toList()]}
		
		pFiles2DirCNS.depends             = pCNVkitSeg
		pFiles2DirCNS.input               = {pFiles2DirCNS.input: lambda ch: [ch.toList()]}
		
		pCNVkitPlot.depends               = [pFiles2DirCNR, pFiles2DirCNS]
		pCNVkitRpt.depends                = [pCNVkitFix, pCNVkitSeg]
		
		pCNVkit2Vcf.depends               = pCNVkitSeg
		
		pVcf2List.depends                 = pVariantFiltration
		pFiles2DirVCF.depends             = pVcf2List
		pFiles2DirVCF.input               = {pFiles2DirVCF.input: lambda ch: [ch.toList()]}
		pCbindList.depends                = pFiles2DirVCF
		pCallRate.depends                 = pCbindList
		
		pyppl({'proc':{'forks':4}}).starts(pDownloadGet, pKnownSites).flowchart().run()

		
	
	
if __name__ == '__main__':
	unittest.main()