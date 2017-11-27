import unittest
from os import path
from pyppl import PyPPL, Channel
from helpers import getfile, procOK, config, procOKIn

from bioprocs.fastx import pFastqSim, pFastq2Sam, pFastqSE2Sam
from bioprocs.sambam import pSam2Bam, pBamMarkdup, pBamRecal, pBamReadGroup, pBamReorder, pBamMerge, pBam2Gmut, pBamPair2Smut, pBam2Cnv, pBamStats, pBam2Fastq, pBam2FastqSE, pBam2Counts

class TestSambam (unittest.TestCase):
	
	@classmethod
	def setUpClass (self):
		
		pFastqSim.args.ref = getfile('ref.fa')
		pFastqSim.args.num = 100000
		pFastqSim.output   = 'fq1:file:read{{in.seed}}_{{args.tool}}_nogit_1.fastq, fq2:file:read{{in.seed}}_{{args.tool}}_nogit_2.fastq'
		pFastqSim.exdir    = getfile()

		pFastqSimWgsim  = pFastqSim.copy ()
		pFastqSimDwgsim = pFastqSim.copy ()
		pFastqSimWgsim.input       = [1]
		pFastqSimDwgsim.input      = [1]
		pFastqSimWgsim.args.tool   = 'wgsim'
		pFastqSimDwgsim.args.tool  = 'dwgsim'
		
		PyPPL(config).start(pFastqSimWgsim, pFastqSimDwgsim).run()

	
	def test2_pFastq2Sam (self):
		pFastq2Sam.args.nthread     = 2
		pFastq2Sam.forks            = 4
		pFastq2Sam.input            = Channel.fromPairs(getfile('*.fastq'))
		pFastq2Sam.args.ref         = getfile('ref.fa')
		pFastq2Sam.tag              = 'bwa'
		pFastq2SamNgm               = pFastq2Sam.copy('ngm')
		pFastq2SamBowtie2           = pFastq2Sam.copy('bowtie2')
		pFastq2Sam.args.tool        = 'bwa'
		pFastq2SamNgm.args.tool     = 'ngm'
		pFastq2SamBowtie2.args.tool = 'bowtie2'
		PyPPL().start(pFastq2Sam, pFastq2SamNgm, pFastq2SamBowtie2).run()
		procOKIn (pFastq2Sam, [
			'chr1_23166_22718_1_0_0_0_0:0:0_2:0:0_0	83	chr1	23166	60	100M	=	22718	-548	AGTGATATTTTTAGGCTTATCTACCAGTTTTAGACGCTCTTTAATATCTTCAGGAATTATTTTATTGTCATATTGTATCATGCTAAATGACAATTTGCTT	12632215/22221/32215233242222224243044222321222020352333214222271302321232.3222215215174423122525422	NM:i:0	MD:Z:100	AS:i:100	XS:i:0	RG:Z:read1_dwgsim_nogit.L0'
		], self)
		procOKIn (pFastq2Sam, [
			'chr2_2473_2109_1_0_0_0_4:0:0_1:0:0_a3ba	163	chr2	2109	60	100M	=	2473	464	GCGTACTCGGCATGGTTGGGAACTAAGAGTGAAGGGCACACAGCCCAGGTTTCAAGCGTTAAGCTCCATGCTTTCTAAGCGGAAGACGGTGATCGTGATC	224247222418616342241232340441220333131/422224245028222/00232226052036230333334440231247441042624322	NM:i:1	MD:Z:5G94	AS:i:95	XS:i:0	RG:Z:read1_dwgsim_nogit.L0'
		], self)
		procOKIn (pFastq2SamNgm, [
			'chr1_23166_22718_1_0_0_0_0:0:0_2:0:0_0	163	chr1	22718	60	100M	=	23166	548	TTGCTTTTAAGACTGAACGCATGAAATATGGTTTTTCGTCATGTTTTGAGTCTGCTGTTGATATTTCTAAAGTAGGTTTTTTTTCCTCGTTTTCTCTAAC	/222622224232432244202220112200332231232/222045/522222216122/440135212345302226354523614224214122340	RG:Z:read1_dwgsim_nogit.L0	AS:i:950	NM:i:2	NH:i:0	XI:f:0.98	X0:i:0	XE:i:21	XR:i:100	MD:Z:73C11T14'
		], self)
		procOKIn (pFastq2SamNgm, [
			'chr2_26823_27211_0_1_0_0_1:0:0_0:0:0_9fd2	99	chr2	26823	60	100M	=	27211	488	CAATCGCCGACAAGATGCGATTCGATATGCGCTAACGCCGGCAAGGCCACGTCTCCACAATGTAAGGTGGCCTTTTCAGATAGACCCCCCACACTGTTGC	130321132/435620113204223243252222.33034245435102112.1230228224304402134222222/231213/40132222-31232	RG:Z:read1_dwgsim_nogit.L0	AS:i:975	NM:i:1	NH:i:0	XI:f:0.99	X0:i:0	XE:i:25	XR:i:100	MD:Z:83G16'
		], self)
		procOKIn (pFastq2SamBowtie2, [
			'chr1_23166_22718_1_0_0_0_0:0:0_2:0:0_0	81	chr1	23166	42	100M	=	22718	-548	AGTGATATTTTTAGGCTTATCTACCAGTTTTAGACGCTCTTTAATATCTTCAGGAATTATTTTATTGTCATATTGTATCATGCTAAATGACAATTTGCTT	12632215/22221/32215233242222224243044222321222020352333214222271302321232.3222215215174423122525422	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:100	YS:i:-7	YT:Z:DP	RG:Z:read1_dwgsim_nogit.L0'
		], self)
		procOKIn (pFastq2SamBowtie2, [
			'chr2_2473_2109_1_0_0_0_4:0:0_1:0:0_a3ba	163	chr2	2109	42	100M	=	2473	464	GCGTACTCGGCATGGTTGGGAACTAAGAGTGAAGGGCACACAGCCCAGGTTTCAAGCGTTAAGCTCCATGCTTTCTAAGCGGAAGACGGTGATCGTGATC	224247222418616342241232340441220333131/422224245028222/00232226052036230333334440231247441042624322	AS:i:-4	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:5G94	YS:i:-12	YT:Z:CP	RG:Z:read1_dwgsim_nogit.L0'
		], self)
				
	def test2_FastqSE2Sam (self):
		pFastqSE2Sam.args.nthread = 2
		pFastqSE2Sam.forks        = 8
		pFastqSE2Sam.input        = Channel.fromPairs(getfile('*.fastq'))
		pFastqSE2Sam.args.ref     = getfile('ref.fa')
		pFastqSE2Sam1             = pFastqSE2Sam.copy('bwa')
		pFastqSE2Sam2             = pFastqSE2Sam.copy('ngm')
		pFastqSE2Sam3             = pFastqSE2Sam.copy('bowtie2')
		pFastqSE2Sam1.args.tool   = 'bwa'
		pFastqSE2Sam2.args.tool   = 'ngm'
		pFastqSE2Sam3.args.tool   = 'bowtie2'
		
		PyPPL().start(pFastqSE2Sam1, pFastqSE2Sam2, pFastqSE2Sam3).run()
		procOKIn (pFastqSE2Sam1, [
			'chr1_1933_2305_0_1_0_0_0:0:0_3:0:0_14	0	chr1	1933	60	100M	*	0	0	AAATGCGCGTATGGGGATGGGGGCCGGGTGAGGAAAGCTGGCTGATTGACCGGCAGATTATTATGGGCCGCCACGACGATGAACAGACGCTGCTGCGTGT	224102142332342226355262-231232231242440502423134242002/21332524234322222655223300414240143231243417	NM:i:0	MD:Z:100	AS:i:100	XS:i:0	RG:Z:read1_dwgsim_nogit_1.L0'
		], self)
		procOKIn (pFastqSE2Sam2, [
			'chr2_20086_19655_1_0_0_0_0:0:0_1:0:0_7ca4/1	16	chr2	20086	60	100M	*	0	0	AATCCCATAACGTAGGACAAGACATTTTCGCGCTGTTTGACGCCTTTATGGGAGGGGAAACGGGCGAACGAAGGGTGGATATGAGACTACTTCTCACAAG	124023213/2313543301232523223245222402524/4680216112122205423222142621035112236322202341202622224332	RG:Z:read1_dwgsim_nogit_1.L0	AS:i:1000	NM:i:0	NH:i:1	XI:f:1	X0:i:1	XE:i:30	XR:i:100	MD:Z:100'
		], self)
		procOKIn (pFastqSE2Sam3, [
			'chr2_30589_30191_1_0_0_0_2:0:0_3:0:0_a3a7/1	16	chr2	30589	42	100M	*	0	0	GTCGAGAAACGCATTCACGGCGGACTGGCGGACCAATGGGCTACACACGCAAGACAACGTCATCACATGTGCGCTTTTGTATATAATTCTTTTGCAAGTC	1020312232542210021123722235351231232313123345522325264325514326342234243531122/432433.21134-0/62332	AS:i:-6	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:40A26G32	YT:Z:UU	RG:Z:read1_dwgsim_nogit_1.L0'
		], self)
	
	def test3_1Sam2BamSambamba (self):
		pSam2Bam.input            = pFastq2Sam.channel
		pSam2Bam.forks            = 2
		pSam2BamSbb1              = pSam2Bam.copy('sambamba')
		pSam2BamSbb2              = pSam2Bam.copy('sambamba')
		pSam2BamSbb3              = pSam2Bam.copy('sambamba')
		pSam2BamSbb1.args.tool    = 'sambamba'
		pSam2BamSbb1.args.markdup = True
		pSam2BamSbb2.args.tool    = 'sambamba'
		pSam2BamSbb2.args.sort    = False
		pSam2BamSbb2.args.index   = False
		pSam2BamSbb2.args.markdup = False
		pSam2BamSbb2.args.rmdup   = False
		pSam2BamSbb3.args.tool    = 'sambamba'
		pSam2BamSbb3.args.sort    = True
		pSam2BamSbb3.args.index   = True
		pSam2BamSbb3.args.markdup = False
		pSam2BamSbb3.args.rmdup   = False
		
		PyPPL().start(pSam2BamSbb1, pSam2BamSbb2, pSam2BamSbb3).run()
			
	def test3_2Sam2BamSamtools (self):
		pSam2Bam.input            = pFastq2Sam.channel
		pSam2Bam.forks            = 2
		pSam2BamSmt1              = pSam2Bam.copy('samtools')
		pSam2BamSmt2              = pSam2Bam.copy('samtools')
		pSam2BamSmt3              = pSam2Bam.copy('samtools')
		pSam2BamSmt1.args.tool    = 'samtools'
		pSam2BamSmt1.args.markdup = True
		pSam2BamSmt2.args.tool    = 'samtools'
		pSam2BamSmt2.args.sort    = False
		pSam2BamSmt2.args.index   = False
		pSam2BamSmt2.args.markdup = False
		pSam2BamSmt2.args.rmdup   = False
		pSam2BamSmt3.args.tool    = 'samtools'
		pSam2BamSmt3.args.sort    = True
		pSam2BamSmt3.args.index   = True
		pSam2BamSmt3.args.markdup = False
		pSam2BamSmt3.args.rmdup   = False
		
		PyPPL().start(pSam2BamSmt1, pSam2BamSmt2, pSam2BamSmt3).run()
			
	def test3_3Sam2BamPicard (self):
		pSam2Bam.input            = pFastq2Sam.channel
		pSam2Bam.forks            = 2
		pSam2BamPcd1              = pSam2Bam.copy('picard')
		pSam2BamPcd2              = pSam2Bam.copy('picard')
		pSam2BamPcd3              = pSam2Bam.copy('picard')
		pSam2BamPcd1.args.tool    = 'picard'
		pSam2BamPcd1.args.markdup = True
		pSam2BamPcd2.args.tool    = 'picard'
		pSam2BamPcd2.args.sort    = False
		pSam2BamPcd2.args.index   = False
		pSam2BamPcd2.args.markdup = False
		pSam2BamPcd2.args.rmdup   = False
		pSam2BamPcd3.args.tool    = 'picard'
		pSam2BamPcd3.args.sort    = True
		pSam2BamPcd3.args.index   = True
		pSam2BamPcd3.args.markdup = False
		pSam2BamPcd3.args.rmdup   = False

		PyPPL().start(pSam2BamPcd1, pSam2BamPcd2, pSam2BamPcd3).run()

	
	def test3_4Sam2BamBiobambam (self):
		pSam2Bam.input                      = pFastq2Sam.channel
		pSam2Bam.forks                      = 2
		pSam2Bam.tag                        = 'biobambam'
		pSam2BamBbb2                        = pSam2Bam.copy('biobambam')
		pSam2BamBbb3                        = pSam2Bam.copy('biobambam')
		pSam2Bam.args.params.maxreadlen     = 5000
		pSam2BamBbb2.args.params.maxreadlen = 5000
		pSam2BamBbb3.args.params.maxreadlen = 5000
		pSam2Bam.args.tool                  = 'biobambam'
		pSam2Bam.args.sort                  = False
		pSam2Bam.args.index                 = False
		pSam2Bam.args.markdup               = False
		pSam2Bam.args.rmdup                 = False
		pSam2BamBbb2.args.tool              = 'biobambam'
		pSam2BamBbb2.args.markdup           = True
		pSam2BamBbb3.args.tool              = 'biobambam'
		pSam2BamBbb3.args.sort              = True
		pSam2BamBbb3.args.index             = True
		pSam2BamBbb3.args.markdup           = False
		pSam2BamBbb3.args.rmdup             = False
		
		PyPPL().start(pSam2Bam, pSam2BamBbb2, pSam2BamBbb3).run()

	def test4_1BamMarkdup (self):
		pBamMarkdup.args.nthread = 2
		pBamMarkdup.forks        = 2
		pBamMarkdup.input        = pSam2Bam.channel
		pBamMarkdup.tag          = 'sambamba'
		pBamMarkdup2             = pBamMarkdup.copy ('biobambam')
		pBamMarkdup3             = pBamMarkdup.copy ('picard')
		pBamMarkdup4             = pBamMarkdup.copy ('samtools')
		pBamMarkdup5             = pBamMarkdup.copy ('bamutil')
		pBamMarkdup.args.tool    = 'sambamba'
		pBamMarkdup2.args.tool   = 'biobambam'
		pBamMarkdup3.args.tool   = 'picard'
		pBamMarkdup4.args.tool   = 'samtools'
		pBamMarkdup5.args.tool   = 'bamutil'
		
		PyPPL().start(pBamMarkdup,pBamMarkdup2,pBamMarkdup4,pBamMarkdup5).run()
	
	def test4_2BamRecal (self):
		pBamRecalTest                                                 = pBamRecal.copy()
		pBamRecalTest.input                                           = pBamMarkdup.channel
		pBamRecalTest.args.ref                                        = getfile('ref.fa')
		#pBamRecalTest.args.knownSites                                 = params.dbsnp
		pBamRecalTest.args.paramsBaseRecalibrator.maximum_cycle_value = 5000
		pBamRecalTest.forks                                           = 2
		
		PyPPL().start(pBamRecalTest).run()
		
	def test4_BamReadGroup (self):
		pBamReadGroup1           = pBamReadGroup.copy()
		pBamReadGroup2           = pBamReadGroup.copy()
		pBamReadGroup1.input     = pBamMarkdup.channel
		pBamReadGroup2.input     = pBamMarkdup.channel
		pBamReadGroup2.args.tool = 'bamutil'
		pBamReadGroup1.forks     = 10
		pBamReadGroup2.forks     = 10
		
		PyPPL().start(pBamReadGroup1, pBamReadGroup2).run()
	
	def test4_BamReorder (self):
		
		pBamReorder.args.ref = getfile('ref.fa')
		pBamReorder.forks    = 10
		pBamReorder.input    = pBamMarkdup.channel
		
		PyPPL().start(pBamReorder).run()
	
	def test4_Mergebams (self):
		
		pBamMerge1 = pBamMerge.copy()
		pBamMerge2 = pBamMerge.copy()
		pBamMerge3 = pBamMerge.copy()
		pBamMerge4 = pBamMerge.copy()
		
		pBamMerge1.input = [pBamMarkdup.channel.flatten(0)]
		pBamMerge2.input = [pBamMarkdup.channel.flatten(0)]
		pBamMerge3.input = [pBamMarkdup.channel.flatten(0)]
		pBamMerge4.input = [pBamMarkdup.channel.flatten(0)]
		
		pBamMerge1.args.tool = 'bamutil'
		pBamMerge2.args.tool = 'samtools'
		pBamMerge3.args.tool = 'picard'
		pBamMerge4.args.tool = 'sambamba'
		
		PyPPL().start(pBamMerge1, pBamMerge2, pBamMerge3, pBamMerge4).run()
		
	def test4_Bam2Gmut (self):
		pBam2Gmut.args.ref     = getfile('ref.fa')
		pBam2Gmut.args.nthread = 16
		pBam2Gmut.forks        = 10
		pBam2Gmut.input        = pBamMarkdup.channel
		
		pBam2Gmut1         = pBam2Gmut.copy()
		pBam2Gmut2         = pBam2Gmut.copy()
		pBam2Gmut3         = pBam2Gmut.copy()
		pBam2Gmut4         = pBam2Gmut.copy()
		pBam2Gmut5         = pBam2Gmut.copy()

		pBam2Gmut1.args.tool = 'gatk'
		pBam2Gmut2.args.tool = 'snvsniffer'
		pBam2Gmut3.args.tool = 'platypus'
		pBam2Gmut4.args.tool = 'vardict'
		pBam2Gmut5.args.tool = 'strelka'
		
		PyPPL().start([
			pBam2Gmut1, 
			pBam2Gmut2, 
			pBam2Gmut3, 
			#pBam2Gmut4, # too slow for test
			pBam2Gmut5
		]).run()
		procOKIn(pBam2Gmut1, 'chr2	23306	.	G	C	1038.77	.	AC=2;AF=1.00;AN=2;BaseQRankSum=0.551;ClippingRankSum=0.000;DP=45;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQRankSum=0.000;QD=23.08;ReadPosRankSum=0.887;SOR=0.303	GT:AD:DP:GQ:PL	1/1:1,44:45:99:1067,138,0', self)
		procOKIn(pBam2Gmut2, 'chr2	27659	.	G	T	.	PASS	DP=199;VT=SNP;HC=4.350	GT	0|1', self)
		procOKIn(pBam2Gmut3, 'chr1	44447	.	G	GA	65	badReads	BRF=0.99;FR=0.9998;HP=1;HapScore=2;MGOF=41;MMLQ=15;MQ=60.0;NF=1;NR=1;PP=65;QD=47.0;SC=TGACGCTCTGGTGGTGCAATG;SbPval=0.74;Source=Platypus;TC=2;TCF=1;TCR=1;TR=2;WE=44455;WS=44437	GT:GL:GOF:GQ:NR:NV	1/1:-9.4,-0.3,0.0:41:5:2:2', self)
		procOKIn(pBam2Gmut5, 'chr2	34994	.	A	.	.	PASS	.	GT:GQX:DP:DPF	0/0:21:8:4', self)
		
	def test4_pBamPair2Smut (self):
		pBamPair2Smut.input        = pBamMarkdup.channel.colAt(0).unfold(2)
		pBamPair2Smut.args.ref     = getfile('ref.fa')
		pBamPair2Smut.args.nthread = 16
		pBamPair2Smut.forks        = 10
		
		pBamPair2Smut1 = pBamPair2Smut.copy()
		pBamPair2Smut2 = pBamPair2Smut.copy()
		pBamPair2Smut3 = pBamPair2Smut.copy()
		pBamPair2Smut4 = pBamPair2Smut.copy()
		pBamPair2Smut5 = pBamPair2Smut.copy()
		pBamPair2Smut6 = pBamPair2Smut.copy()
		
		pBamPair2Smut1.args.tool = 'gatk'
		pBamPair2Smut2.args.tool = 'somaticsniper'
		pBamPair2Smut3.args.tool = 'snvsniffer'
		pBamPair2Smut4.args.tool = 'strelka'
		pBamPair2Smut5.args.tool = 'virmid'
		pBamPair2Smut6.args.tool = 'vardict'
		pBamPair2Smut4.args.gz   = True
		
		PyPPL().start([
			pBamPair2Smut1, 
			pBamPair2Smut2, 
			pBamPair2Smut3, 
			pBamPair2Smut4, 
			#pBamPair2Smut5, # itself has error for test
			#pBamPair2Smut6, # too slow for test
		]).run()
		procOKIn(pBamPair2Smut1, '##INFO=<ID=TLOD,Number=1,Type=String,Description="Tumor LOD score">', self)
		procOKIn(pBamPair2Smut2, 'chr2	29452	.	G	A	.	.	.	GT:IGT:DP:DP4:BCOUNT:GQ:JGQ:VAQ:BQ:MQ:AMQ:SS:SSC	0/1:0/1:230:50,66,64,50:110,4,116,0:111:.:111:17,17:60:60,60:1:.	0/0:0/0:220:98,119,1,2:0,1,217,2:216:.:230:17:60:60:3:108', self)
		procOKIn(pBamPair2Smut3, 'chr2	34754	.	A	T	.	PASS	DP=234;VT=SNP;SS=SOMATIC;LC=0.000;PVAL=0.123333	GT:DP	0/0:121	0/1:113', self)
		
	def test5_pBam2Cnv (self):
		pBam2Cnv.input = pBamMarkdup.channel
		
		pBam2Cnv1           = pBam2Cnv.copy()
		pBam2Cnv1.forks     = 10
		pBam2Cnv1.args.tool = 'cnvkit'
		pBam2Cnv1.args.ref  = getfile('ref.fa')
		
		''' 
		#until cnvnator is fixed
		pBam2Cnv2           = pBam2Cnv.copy()
		pBam2Cnv2.forks     = 10
		pBam2Cnv2.args.tool = 'cnvnator'
		pBam2Cnv2.args.ref  = getfile('ref.fa')
		'''
		PyPPL().start([
			pBam2Cnv1,
			#pBam2Cnv2
		]).run()
		procOKIn(pBam2Cnv1, '##FORMAT=<ID=CNQ,Number=1,Type=Float,Description="Copy number genotype quality for imprecise events">', self)
		
	def test6_pBamStats (self):
		pBamStats.input = pBamMarkdup.channel
		pBamStats.forks = 10
		
		pBamStats1 = pBamStats.copy()
		pBamStats1.args.plot = True
		
		PyPPL().start(pBamStats1).run()
		procOK(pBamStats1, 'bamstats.txt', self)
		
	def test7_pBam2Fastq (self):
		
		pBam2FastqBmm       = pBam2Fastq.copy(newid = 'pBam2FastqTest', tag = 'biobambam')
		pBam2FastqBmm.input = pBamMarkdup.channel
		pBam2FastqBmm.forks = 10
		
		pBam2FastqBtl           = pBam2FastqBmm.copy(newid = 'pBam2FastqTest', tag = 'bedtools')
		pBam2FastqBtl.args.tool = 'bedtools'
		
		pBam2FastqStl           = pBam2FastqBmm.copy(newid = 'pBam2FastqTest', tag = 'samtools')
		pBam2FastqStl.args.tool = 'samtools'
		
		pBam2FastqPcd           = pBam2FastqBmm.copy(newid = 'pBam2FastqTest', tag = 'picard')
		pBam2FastqPcd.args.tool = 'picard'
		
		PyPPL().start('pBam2FastqTest').run()
		procOKIn(pBam2FastqBmm, '@chr2_34349_34706_0_1_0_0_2:0:0_1:0:0_8382/1', self)
		procOKIn(pBam2FastqBtl, 'GGGCCGCCCCAGCCATCCCTTCACAGCTTTTCTGCCATGCGAGCTTACTGCCTTGTCCTGGGCCCTCTTGAGCTGAACGTTCGTAAGTGTCACGGACAGC', self)
		procOKIn(pBam2FastqStl, 'ACCAGAATCGGAATCTTAAGGTCACCTTCATGCAATGCTCTGATATGGTCAAAGGCCTGGATAGTTTTGACATGCATGATGGAACGACCAGGCTCCGATT', self)
		procOKIn(pBam2FastqPcd, 'TATATCGCGGTCCCTGCGAGATGTATTACGTGATTCGGCCGTCCCGGATTGTTGGGCGCCGATCCAGAAGAGCTAAGAGCAGAGGAGTGTTTAGGCATCT', self)
		
	def test7_pBam2FastqSE (self):
		
		pBam2FastqSEBmm       = pBam2FastqSE.copy(newid = 'pBam2FastqSETest', tag = 'biobambam')
		pBam2FastqSEBmm.input = pBamMarkdup.channel
		pBam2FastqSEBmm.forks = 10
		
		pBam2FastqSEBtl           = pBam2FastqSEBmm.copy(newid = 'pBam2FastqSETest', tag = 'bedtools')
		pBam2FastqSEBtl.args.tool = 'bedtools'
		
		pBam2FastqSEStl           = pBam2FastqSEBmm.copy(newid = 'pBam2FastqSETest', tag = 'samtools')
		pBam2FastqSEStl.args.tool = 'samtools'
		
		#pBam2FastqSEPcd       = pBam2FastqSEBmm.copy(newid = 'pBam2FastqSETest', tag = 'picard')
		#pBam2FastqSEPcd.args.tool = 'picard'
		
		PyPPL().start('pBam2FastqSETest').run()
		procOKIn(pBam2FastqSEBmm, '@chr2_34788_34437_1_0_0_0_1:1:0_3:0:0_48f9/1', self)
		procOKIn(pBam2FastqSEBtl, 'TATATCGCGGTCCCTGCGAGATGTATTACGTGATTCGGCCGTCCCGGATTGTTGGGCGCCGATCCAGAAGAGCTAAGAGCAGAGGAGTGTTTAGGCATCT', self)
		procOKIn(pBam2FastqSEStl, 'TCCTGTCGTTTATTGGCCGCGATCCGAAGCTTGTATAGCAGATAACTCCCGCACCCGCGCGTGTACTTTAGGCGCTTGCTGCAACACCCGTATCATCCAT', self)


	def test8_pBam2Counts(self):
		pBam2Counts.input        = pBamMarkdup.channel
		pBam2Counts.args.refgene = getfile('gene.gtf')
		pBam2Counts.forks        = 10
		PyPPL().start(pBam2Counts).run()
		procOKIn(pBam2Counts, '__no_feature', self)

if __name__ == '__main__':
	unittest.main(failfast=True)