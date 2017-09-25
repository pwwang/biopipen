import os, sys, unittest, addPath

from pyppl import PyPPL, Channel, Box
from bioprocs.fastx import pFastqPESim, pFastQC, pFastMC, pFastqPETrim, pFastqSETrim, pFastqSE2Sam, pFastqPE2Sam
from bioprocs.common import pFiles2Dir
from bioaggrs import params

params = params.toDict()

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

class TestFastx (unittest.TestCase):
	
	data = Box()
	def test1_FastqPESim (self):
		
		pFastqPESim.input      = list(range(2))
		pFastqPESim.args.ref   = params.ref
		pFastqPESim.args.num   = 1000
		pFastqPESim.forks      = 2
		pFastqPESim1           = pFastqPESim.copy()
		pFastqPESim2           = pFastqPESim.copy()
		pFastqPESim1.args.tool = 'wgsim'
		pFastqPESim2.args.tool = 'dwgsim'
		pFastqPESim1.callback  = lambda p: self.data.update({'fqs': p.channel.fold(1).flatten()})
		
		PyPPL().start(pFastqPESim1, pFastqPESim2).run()
		
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
		pFastqPETrim.input = Channel.create(self.data.fqs).unfold(2)
		pFastqPETrim.args.gz = False
		pFastqPETrim.args.nthread = 2
		pFastqPETrim.forks = 2
		pFastqPETrim1 = pFastqPETrim.copy()
		pFastqPETrim2 = pFastqPETrim.copy()
		pFastqPETrim3 = pFastqPETrim.copy()
		pFastqPETrim1.args.tool = 'trimmomatic'
		pFastqPETrim2.args.tool = 'cutadapt'
		pFastqPETrim3.args.tool = 'skewer'
		
		PyPPL().start(pFastqPETrim1, pFastqPETrim2, pFastqPETrim3).run()
		
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
		
	def test2_pFastqPE2Sam (self):
		pFastqPE2Sam.args.nthread     = 2
		pFastqPE2Sam.forks            = 4
		pFastqPE2Sam.input            = Channel.create(self.data.fqs).unfold(2)
		pFastqPE2Sam.args.ref         = params.ref
		pFastqPE2SamBwa               = pFastqPE2Sam.copy('bwa')
		pFastqPE2SamNgm               = pFastqPE2Sam.copy('ngm')
		pFastqPE2SamBowtie2           = pFastqPE2Sam.copy('bowtie2')
		pFastqPE2SamBwa.args.tool     = 'bwa'
		pFastqPE2SamNgm.args.tool     = 'ngm'
		pFastqPE2SamBowtie2.args.tool = 'bowtie2'
		
		PyPPL().start(pFastqPE2SamBwa, pFastqPE2SamNgm, pFastqPE2SamBowtie2).run()
		
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
		
	
if __name__ == '__main__':
	unittest.main(failfast=True)