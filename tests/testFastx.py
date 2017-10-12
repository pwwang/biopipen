import os, sys, unittest, addPath

from pyppl import PyPPL, Channel, Box
from bioprocs import params
from bioprocs.fastx import pFastqSim, pFastQC, pFastMC, pFastqTrim, pFastqSETrim, pFastqSE2Sam, pFastq2Sam, pFastq2Expr
from bioprocs.common import pFiles2Dir

params = params.toDict()

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

class TestFastx (unittest.TestCase):
	
	data = Box()
	def test1_FastqPESim (self):
		
		pFastqSim.input      = list(range(2))
		pFastqSim.args.ref   = params.ref
		pFastqSim.args.num   = 1000
		pFastqSim.forks      = 2
		pFastqSim1           = pFastqSim.copy()
		pFastqSim2           = pFastqSim.copy()
		pFastqSim1.args.tool = 'wgsim'
		pFastqSim2.args.tool = 'dwgsim'
		pFastqSim1.callback  = lambda p: self.data.update({'fqs': p.channel.fold(1).flatten()})
		
		PyPPL().start(pFastqSim1, pFastqSim2).run()
		
	def test2_FastQC (self):
		pFastQC.input = self.data.fqs
		pFastQC.forks = 4
		pFastQC.callback = lambda p: self.data.update({'qcdirs': p.channel.flatten()})
		PyPPL().start(pFastQC).run()
		
	def test3_MultiQC (self):
		pFiles2Dir.input = [self.data.qcdirs]
		pFastMC.depends = pFiles2Dir
		PyPPL().start(pFiles2Dir).run()
		
	def test4_PETrim (self):
		pFastqTrim.input = Channel.create(self.data.fqs).unfold(2)
		pFastqTrim.args.gz = False
		pFastqTrim.args.nthread = 2
		pFastqTrim.forks = 2
		pFastqTrim1 = pFastqTrim.copy()
		pFastqTrim2 = pFastqTrim.copy()
		pFastqTrim3 = pFastqTrim.copy()
		pFastqTrim1.args.tool = 'trimmomatic'
		pFastqTrim2.args.tool = 'cutadapt'
		pFastqTrim3.args.tool = 'skewer'
		
		PyPPL().start(pFastqTrim1, pFastqTrim2, pFastqTrim3).run()
		
	def test4_SETrim (self):
		pFastqSETrim.input = self.data.fqs
		pFastqSETrim.args.gz = False
		pFastqSETrim.args.nthread = 2
		pFastqSETrim.forks = 4
		pFastqSETrim1 = pFastqSETrim.copy()
		pFastqSETrim2 = pFastqSETrim.copy()
		pFastqSETrim3 = pFastqSETrim.copy()
		pFastqSETrim1.args.tool = 'trimmomatic'
		pFastqSETrim2.args.tool = 'cutadapt'
		pFastqSETrim3.args.tool = 'skewer'
		
		PyPPL().start(pFastqSETrim1, pFastqSETrim2, pFastqSETrim3).run()
		
	def test2_pFastq2Sam (self):
		pFastq2Sam.args.nthread     = 2
		pFastq2Sam.forks            = 4
		pFastq2Sam.input            = Channel.create(self.data.fqs).unfold(2)
		pFastq2Sam.args.ref         = params.ref
		pFastq2SamBwa               = pFastq2Sam.copy('bwa')
		pFastq2SamNgm               = pFastq2Sam.copy('ngm')
		pFastq2SamBowtie2           = pFastq2Sam.copy('bowtie2')
		pFastq2SamBwa.args.tool     = 'bwa'
		pFastq2SamNgm.args.tool     = 'ngm'
		pFastq2SamBowtie2.args.tool = 'bowtie2'
		
		PyPPL().start(pFastq2SamBwa, pFastq2SamNgm, pFastq2SamBowtie2).run()
		
	def test2_pFastqSE2Sam (self):
		pFastqSE2Sam.args.nthread     = 2
		pFastqSE2Sam.forks            = 8
		pFastqSE2Sam.input            = Channel.create(self.data.fqs)
		pFastqSE2Sam.args.ref         = params.ref
		pFastqSE2SamBwa               = pFastqSE2Sam.copy('bwa')
		pFastqSE2SamNgm               = pFastqSE2Sam.copy('ngm')
		pFastqSE2SamBowtie2           = pFastqSE2Sam.copy('bowtie2')
		pFastqSE2SamBwa.args.tool     = 'bwa'
		pFastqSE2SamNgm.args.tool     = 'ngm'
		pFastqSE2SamBowtie2.args.tool = 'bowtie2'
		
		PyPPL().start(pFastqSE2SamBwa, pFastqSE2SamNgm, pFastqSE2SamBowtie2).run()

	def test5_pFastq2Expr(self):
		pFastq2Expr.input = Channel.create(self.data.fqs).unfold(2)
		pFastq2Expr.forks = 10
		pFastq2Expr.args.nthread = 10
		PyPPL().start(pFastq2Expr).run()
		
	
if __name__ == '__main__':
	unittest.main(failfast=True)