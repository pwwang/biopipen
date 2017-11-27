import unittest
from os import path
from pyppl import PyPPL
from helpers import getfile, procOK, config
from bioprocs.fastx import pFastqSim, pFastQC, pFastMC, pFastqTrim, pFastqSETrim, pFastqSE2Sam, pFastq2Sam, pFastq2Expr
from bioprocs.common import pFiles2Dir
from bioprocs import params

class TestFastx (unittest.TestCase):
	
	def test1_FastqPESim (self):
		
		pFastqSim.input      = [8525]
		pFastqSim.args.ref   = params.ref.value
		pFastqSim.args.num   = 1000
		pFastqSim.args.gz    = True
		pFastqSim.forks      = 2
		pFastqSim1           = pFastqSim.copy()
		pFastqSim2           = pFastqSim.copy()
		pFastqSim1.args.tool = 'wgsim'
		pFastqSim2.args.tool = 'dwgsim'
		
		PyPPL().start(pFastqSim1, pFastqSim2).run()
		procOK(pFastqSim1, 'fastqsim_wgsim_1.fastq.gz', self)
		procOK(pFastqSim2, 'fastqsim_dwgsim_1.fastq.gz', self)
		
	def test2_FastQC (self):
		pFastQC1 = pFastQC.copy()
		pFastQC2 = pFastQC.copy()
		pFastQC1.input = [getfile('fastqsim_dwgsim_1.fastq.gz', input = False)]
		pFastQC2.input = [getfile('fastqsim_wgsim_1.fastq.gz', input = False)]
		
		PyPPL().start(pFastQC1, pFastQC2).run()
		procOK(pFastQC1, 'fastqc-out1', self)
		procOK(pFastQC2, 'fastqc-out2', self)
		
	def test3_MultiQC (self):

		pFiles2Dir.input = [[getfile('fastqc-out1', input = False), getfile('fastqc-out2', input = False)]]
		pFastMC.depends = pFiles2Dir
		PyPPL().start(pFiles2Dir).run()
		procOK(pFastMC, 'multiqc-out', self)
	
	def test4_PETrim (self):
		pFastqTrim.input = (getfile('fastqsim_dwgsim_1.fastq.gz', input = False), getfile('fastqsim_dwgsim_2.fastq.gz', input = False))
		pFastqTrim.args.gz = True
		pFastqTrim.args.nthread = 2
		pFastqTrim1 = pFastqTrim.copy()
		pFastqTrim2 = pFastqTrim.copy()
		pFastqTrim3 = pFastqTrim.copy()
		pFastqTrim1.args.tool = 'trimmomatic'
		pFastqTrim2.args.tool = 'cutadapt'
		pFastqTrim3.args.tool = 'skewer'
		
		PyPPL().start(pFastqTrim1, pFastqTrim2, pFastqTrim3).run()
		procOK(pFastqTrim1, 'fastqtrim-pe-trimmometic_1.fastq.gz', self)
		procOK(pFastqTrim2, 'fastqtrim-pe-cutadapt_1.fastq.gz', self)
		procOK(pFastqTrim3, 'fastqtrim-pe-skewer_1.fastq.gz', self)
		
	def test4_SETrim (self):
		pFastqSETrim.input = [getfile('fastqsim_dwgsim_1.fastq.gz', input = False)]
		pFastqSETrim.args.gz = True
		pFastqSETrim.args.nthread = 2
		pFastqSETrim1 = pFastqSETrim.copy()
		pFastqSETrim2 = pFastqSETrim.copy()
		pFastqSETrim3 = pFastqSETrim.copy()
		pFastqSETrim1.args.tool = 'trimmomatic'
		pFastqSETrim2.args.tool = 'cutadapt'
		pFastqSETrim3.args.tool = 'skewer'
		
		PyPPL().start(pFastqSETrim1, pFastqSETrim2, pFastqSETrim3).run()
		procOK(pFastqSETrim1, 'fastqtrim-se-trimmometic_1.fastq.gz', self)
		procOK(pFastqSETrim2, 'fastqtrim-se-cutadapt_1.fastq.gz', self)
		procOK(pFastqSETrim3, 'fastqtrim-se-skewer_1.fastq.gz', self)
		
	def test2_pFastq2Sam (self):
		pFastq2Sam.args.nthread     = 2
		pFastq2Sam.input            = (getfile('fastqsim_dwgsim_1.fastq.gz', input = False), getfile('fastqsim_dwgsim_2.fastq.gz', input = False))
		pFastq2Sam.args.outfmt      = 'bam'
		pFastq2Sam.args.ref         = params.ref.value
		pFastq2SamBwa               = pFastq2Sam.copy('bwa')
		pFastq2SamNgm               = pFastq2Sam.copy('ngm')
		pFastq2SamBowtie2           = pFastq2Sam.copy('bowtie2')
		pFastq2SamBwa.args.tool     = 'bwa'
		pFastq2SamNgm.args.tool     = 'ngm'
		pFastq2SamBowtie2.args.tool = 'bowtie2'
		
		PyPPL().start(pFastq2SamBwa, pFastq2SamNgm, pFastq2SamBowtie2).run()
		procOK(pFastq2SamBwa, 'fastq2sam_bwa.bam', self)
		procOK(pFastq2SamNgm, 'fastq2sam_ngm.bam', self)
		procOK(pFastq2SamBowtie2, 'fastq2sam_bowtie2.bam', self)
		
	def test2_pFastqSE2Sam (self):
		pFastqSE2Sam.args.nthread     = 2
		pFastqSE2Sam.input            = [getfile('fastqsim_dwgsim_1.fastq.gz', input = False)]
		pFastqSE2Sam.args.outfmt      = 'bam'
		pFastqSE2Sam.args.ref         = params.ref.value
		pFastqSE2SamBwa               = pFastqSE2Sam.copy('bwa')
		pFastqSE2SamNgm               = pFastqSE2Sam.copy('ngm')
		pFastqSE2SamBowtie2           = pFastqSE2Sam.copy('bowtie2')
		pFastqSE2SamBwa.args.tool     = 'bwa'
		pFastqSE2SamNgm.args.tool     = 'ngm'
		pFastqSE2SamBowtie2.args.tool = 'bowtie2'
		
		PyPPL().start(pFastqSE2SamBwa, pFastqSE2SamNgm, pFastqSE2SamBowtie2).run()
		procOK(pFastqSE2SamBwa, 'fastqse2sam_bwa.bam', self)
		procOK(pFastqSE2SamNgm, 'fastqse2sam_ngm.bam', self)
		procOK(pFastqSE2SamBowtie2, 'fastqse2sam_bowtie2.bam', self)

	def test5_pFastq2Expr(self):
		pFastq2Expr.input = (getfile('fastqsim_dwgsim_1.fastq.gz', input = False), getfile('fastqsim_dwgsim_2.fastq.gz', input = False))
		pFastq2Expr.args.nthread = 10
		PyPPL().start(pFastq2Expr).run()
		procOK(pFastq2Expr, 'fastq2expr.txt', self)
	
if __name__ == '__main__':
	unittest.main(failfast=True)