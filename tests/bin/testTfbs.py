import unittest
from os import path
from bioprocs import params
from helpers import getfile, getbin, cmdOK, fileOK

pwmscan    = getbin('pwmscan')
tffile     = getfile('tfs.txt')
motiffile  = getfile('motifs.txt')
genefile   = getfile('genes.txt')
regionfile = getfile('regions.bed')
tfmotifs   = getfile('tf-motifs.txt')
class testTFBS (unittest.TestCase):
	def testPWMScanNoArgs(self):
		cmdOK([pwmscan], self, inerr = 'The target file. If target is region, using BED format.', testrc = False)
		cmdOK([pwmscan, '-tfile', tffile], self, inerr = 'ERROR: Option -infile is required.', testrc = False)

	def testPWMScanTfP(self):
		outdir = path.join(params.tmpdir.value, 'pwmscan-tfbstfp')
		cmd = [pwmscan, '-infile', tffile, '-tfile', genefile, '-outdir', outdir, '-consvp', '0', '-pval', '1e-3']
		cmdOK(cmd, self)
		fileOK(path.join(outdir, 'genes-promoters-human.etc.simread.bed'), getfile('aTfbsTfP.bed', False), self)
	
	def testPWMScanTfPC(self):
		outdir = path.join(params.tmpdir.value, 'pwmscan-tfbstfpc')
		cmd = [pwmscan, '-infile', tffile, '-tfile', genefile, '-outdir', outdir, '-consvp', '0.05', '-pval', '1e-3']
		cmdOK(cmd, self)
		fileOK(path.join(outdir, 'genes-promoters-human.etc.simread-consv.bed'), getfile('aTfbsTfPC.bed', False), self)

	def testPWMScanTfR(self):
		outdir = path.join(params.tmpdir.value, 'pwmscan-tfbstfr')
		cmd = [pwmscan, '-infile', tffile, '-target', 'region', '-tfile', regionfile, '-outdir', outdir, '-consvp', '0', '-pval', '1e-3']
		cmdOK(cmd, self)
		fileOK(path.join(outdir, 'regions-human.etc.simread.bed'), getfile('aTfbsTfR.bed', False), self)

	def testPWMScanTfRC(self):
		outdir = path.join(params.tmpdir.value, 'pwmscan-tfbstfrc')
		cmd = [pwmscan, '-infile', tffile, '-target', 'region', '-tfile', regionfile, '-outdir', outdir, '-consvp', '0.05', '-pval', '1e-3']
		cmdOK(cmd, self)
		fileOK(path.join(outdir, 'regions-human.etc.simread-consv.bed'), getfile('aTfbsTfRC.bed', False), self)

	def testPWMScanP(self):
		outdir = path.join(params.tmpdir.value, 'pwmscan-tfbsp')
		cmd = [pwmscan, '-infile', tfmotifs, '-input', 'motif', '-tfile', genefile, '-outdir', outdir, '-consvp', '0', '-pval', '1e-3']
		cmdOK(cmd, self)
		fileOK(path.join(outdir, 'genes-promoters-tf-motifs.bed'), getfile('aTfbsP.bed', False), self)
	
	def testPWMScanPC(self):
		outdir = path.join(params.tmpdir.value, 'pwmscan-tfbspc')
		cmd = [pwmscan, '-infile', tfmotifs, '-input', 'motif', '-tfile', genefile, '-outdir', outdir, '-consvp', '0.05', '-pval', '1e-3']
		cmdOK(cmd, self)
		fileOK(path.join(outdir, 'genes-promoters-tf-motifs-consv.bed'), getfile('aTfbsPC.bed', False), self)

	def testPWMScanR(self):
		outdir = path.join(params.tmpdir.value, 'pwmscan-tfbsr')
		cmd = [pwmscan, '-infile', tfmotifs, '-input', 'motif', '-target', 'region', '-tfile', regionfile, '-outdir', outdir, '-consvp', '0', '-pval', '1e-3']
		cmdOK(cmd, self)
		fileOK(path.join(outdir, 'regions-tf-motifs.bed'), getfile('aTfbsR.bed', False), self)

	def testPWMScanRC(self):
		outdir = path.join(params.tmpdir.value, 'pwmscan-tfbsrc')
		cmd = [pwmscan, '-infile', tfmotifs, '-input', 'motif', '-target', 'region', '-tfile', regionfile, '-outdir', outdir, '-consvp', '0.1', '-pval', '1e-3']
		cmdOK(cmd, self)
		fileOK(path.join(outdir, 'regions-tf-motifs-consv.bed'), getfile('aTfbsRC.bed', False), self)

if __name__ == '__main__':
	unittest.main(failfast = True, verbosity = 2)