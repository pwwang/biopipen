import unittest, addPath

import random
from os import path
from tempfile import gettempdir
from pyppl import PyPPL, Channel, Proc
from bioprocs.gsea import pTargetEnrichr, pEnrichr, pExpmat2Gct, pSSGSEA, pGSEA, pSampleinfo2Cls, pGMT2Mat

unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

tmpdir   = gettempdir()
genelist = ['Nsun3', 'Polrmt', 'Nlrx1', 'Sfxn5', 'Zc3h12c', 'Slc25a39', 'Arsg', 'Defb29', 'Ndufb6', 'Zfand1', 'Tmem77', '5730403B10Rik', 'RP23-195K8.6', 'Tlcd1', 'Psmc6', 'Slc30a6', 'LOC100047292', 'Lrrc40', 'Orc5l', 'Mpp7', 'Unc119b', 'Prkaca', 'Tcn2', 'Psmc3ip', 'Pcmtd2', 'Acaa1a', 'Lrrc1', '2810432D09Rik', 'Sephs2', 'Sac3d1', 'Tmlhe', 'LOC623451', 'Tsr2', 'Plekha7', 'Gys2', 'Arhgef12', 'Hibch', 'Lyrm2', 'Zbtb44', 'Entpd5', 'Rab11fip2', 'Lipt1', 'Intu', 'Anxa13', 'Klf12', 'Sat2', 'Gal3st2', 'Vamp8', 'Fkbpl', 'Aqp11', 'Trap1', 'Pmpcb', 'Tm7sf3', 'Rbm39', 'Bri3', 'Kdr', 'Zfp748', 'Nap1l1', 'Dhrs1', 'Lrrc56', 'Wdr20a', 'Stxbp2', 'Klf1', 'Ufc1', 'Ccdc16', '9230114K14Rik', 'Rwdd3', '2610528K11Rik', 'Aco1', 'Cables1', 'LOC100047214', 'Yars2', 'Lypla1', 'Kalrn', 'Gyk', 'Zfp787', 'Zfp655', 'Rabepk', 'Zfp650', '4732466D17Rik', 'Exosc4', 'Wdr42a', 'Gphn', '2610528J11Rik', '1110003E01Rik', 'Mdh1', '1200014M14Rik', 'AW209491', 'Mut', '1700123L14Rik', '2610036D13Rik', 'Cox15', 'Tmem30a', 'Nsmce4a', 'Tm2d2', 'Rhbdd3', 'Atxn2', 'Nfs1', '3110001I20Rik', 'BC038156', 'LOC100047782', '2410012H22Rik', 'Rilp', 'A230062G08Rik', 'Pttg1ip', 'Rab1', 'Afap1l1', 'Lyrm5', '2310026E23Rik', 'C330002I19Rik', 'Zfyve20', 'Poli', 'Tomm70a', 'Slc7a6os', 'Mat2b', '4932438A13Rik', 'Lrrc8a', 'Smo', 'Nupl2', 'Trpc2', 'Arsk', 'D630023B12Rik', 'Mtfr1', '5730414N17Rik', 'Scp2', 'Zrsr1', 'Nol7', 'C330018D20Rik', 'Ift122', 'LOC100046168', 'D730039F16Rik', 'Scyl1', '1700023B02Rik', '1700034H14Rik', 'Fbxo8', 'Paip1', 'Tmem186', 'Atpaf1', 'LOC100046254', 'LOC100047604', 'Coq10a', 'Fn3k', 'Sipa1l1', 'Slc25a16', 'Slc25a40', 'Rps6ka5', 'Trim37', 'Lrrc61', 'Abhd3', 'Gbe1', 'Parp16', 'Hsd3b2', 'Esm1', 'Dnajc18', 'Dolpp1', 'Lass2', 'Wdr34', 'Rfesd', 'Cacnb4', '2310042D19Rik', 'Srr', 'Bpnt1', '6530415H11Rik', 'Clcc1', 'Tfb1m', '4632404H12Rik', 'D4Bwg0951e', 'Med14', 'Adhfe1', 'Thtpa', 'Cat', 'Ell3', 'Akr7a5', 'Mtmr14', 'Timm44', 'Sf1', 'Ipp', 'Iah1', 'Trim23', 'Wdr89', 'Gstz1', 'Cradd', '2510006D16Rik', 'Fbxl6', 'LOC100044400', 'Zfp106', 'Cd55', '0610013E23Rik', 'Afmid', 'Tmem86a', 'Aldh6a1', 'Dalrd3', 'Smyd4', 'Nme7', 'Fars2', 'Tasp1', 'Cldn10', 'A930005H10Rik', 'Slc9a6', 'Adk', 'Rbks', '2210016F16Rik', 'Vwce', '4732435N03Rik', 'Zfp11', 'Vldlr', '9630013D21Rik', '4933407N01Rik', 'Fahd1', 'Mipol1', '1810019D21Rik', '1810049H13Rik', 'Tfam', 'Paics', '1110032A03Rik', 'LOC100044139', 'Dnajc19', 'BC016495', 'A930041I02Rik', 'Rqcd1', 'Usp34', 'Zcchc3', 'H2afj', 'Phf7', '4921508D12Rik', 'Kmo', 'Prpf18', 'Mcat', 'Txndc4', '4921530L18Rik', 'Vps13b', 'Scrn3', 'Tor1a', 'AI316807', 'Acbd4', 'Fah', 'Apool', 'Col4a4', 'Lrrc19', 'Gnmt', 'Nr3c1', 'Sip1', 'Ascc1', 'Fech', 'Abhd14a', 'Arhgap18', '2700046G09Rik', 'Yme1l1', 'Gk5', 'Glo1', 'Sbk1', 'Cisd1', '2210011C24Rik', 'Nxt2', 'Notum', 'Ankrd42', 'Ube2e1', 'Ndufv1', 'Slc33a1', 'Cep68', 'Rps6kb1', 'Hyi', 'Aldh1a3', 'Mynn', '3110048L19Rik', 'Rdh14', 'Proz', 'Gorasp1', 'LOC674449', 'Zfp775', '5430437P03Rik', 'Npy', 'Adh5', 'Sybl1', '4930432O21Rik', 'Nat9', 'LOC100048387', 'Mettl8', 'Eny2', '2410018G20Rik', 'Pgm2', 'Fgfr4', 'Mobkl2b', 'Atad3a', '4932432K03Rik', 'Dhtkd1', 'Ubox5', 'A530050D06Rik', 'Zdhhc5', 'Mgat1', 'Nudt6', 'Tpmt', 'Wbscr18', 'LOC100041586', 'Cdk5rap1', '4833426J09Rik', 'Myo6', 'Cpt1a', 'Gadd45gip1', 'Tmbim4', '2010309E21Rik', 'Asb9', '2610019F03Rik', '7530414M10Rik', 'Atp6v1b2', '2310068J16Rik', 'Ddt', 'Klhdc4', 'Hpn', 'Lifr', 'Ovol1', 'Nudt12', 'Cdan1', 'Fbxo9', 'Fbxl3', 'Hoxa7', 'Aldh8a1', '3110057O12Rik', 'Abhd11', 'Psmb1', 'ENSMUSG00000074286', 'Chpt1', 'Oxsm', '2310009A05Rik', '1700001L05Rik', 'Zfp148', '39509', 'Mrpl9', 'Tmem80', '9030420J04Rik', 'Naglu', 'Plscr2', 'Agbl3', 'Pex1', 'Cno', 'Neo1', 'Asf1a', 'Tnfsf5ip1', 'Pkig', 'AI931714', 'D130020L05Rik', 'Cntd1', 'Clec2h', 'Zkscan1', '1810044D09Rik', 'Mettl7a', 'Siae', 'Fbxo3', 'Fzd5', 'Tmem166', 'Tmed4', 'Gpr155', 'Rnf167', 'Sptlc1', 'Riok2', 'Tgds', 'Pms1', 'Pitpnc1', 'Pcsk7', '4933403G14Rik', 'Ei24', 'Crebl2', 'Tln1', 'Mrpl35', '2700038C09Rik', 'Ubie', 'Osgepl1', '2410166I05Rik', 'Wdr24', 'Ap4s1', 'Lrrc44', 'B3bp', 'Itfg1', 'Dmxl1', 'C1d']
reglist = ['Regulator1', 'Regulator2', 'Regulator3', 'Regulator4', 'Regulator5', 'Regulator6', 'Regulator7', 'Regulator8', 'Regulator9', 'Regulator10']
random.seed(8525)
relations = {}
for reg in reglist:
	relations[reg] = random.sample(genelist, 30)

class TestCommon (unittest.TestCase):

	def testpTargetEnrichrCol2 (self):
		infile = path.join(tmpdir, 'targetEnrichrCol2.txt')
		filestr = ''
		for reg, genes in relations.items():
			for g in genes:
				filestr += "%s\t%s\n" % (reg, g)
		ostr = ''
		if path.isfile(infile):
			with open(infile) as f:
				ostr = f.read()
		if not path.isfile(infile) or filestr != ostr:
			with open(infile, 'w') as f:
				f.write(filestr)

		pTargetEnrichrCol2        = pTargetEnrichr.copy()
		pTargetEnrichrCol2.input  = [infile]
		pTargetEnrichrCol2.errhow = 'terminate'
		PyPPL().start(pTargetEnrichrCol2).run()
	
	def testpTargetEnrichrCol3 (self):
		infile = path.join(tmpdir, 'targetEnrichrCol3.txt')
		filestr = ''
		for reg, genes in relations.items():
			for g in genes:
				filestr += "%s\t%s\t%s\n" % (reg, g, random.choice(['+', '-']))
		ostr = ''
		if path.isfile(infile):
			with open(infile) as f:
				ostr = f.read()
		if not path.isfile(infile) or filestr != ostr:
			with open(infile, 'w') as f:
				f.write(filestr)

		pTargetEnrichrCol3        = pTargetEnrichr.copy()
		pTargetEnrichrCol3.input  = [infile]
		pTargetEnrichrCol3.errhow = 'terminate'
		PyPPL().start(pTargetEnrichrCol3).run()

	def testpTargetEnrichrCol4 (self):
		infile = path.join(tmpdir, 'targetEnrichrCol4.txt')
		filestr = ''
		for reg, genes in relations.items():
			for g in genes:
				filestr += "%s\t%s\t%s\t%s\n" % (reg, g, random.choice(['+', '-']), random.choice(['+', '-']))
		ostr = ''
		if path.isfile(infile):
			with open(infile) as f:
				ostr = f.read()
		if not path.isfile(infile) or filestr != ostr:
			with open(infile, 'w') as f:
				f.write(filestr)

		pTargetEnrichrCol4        = pTargetEnrichr.copy()
		pTargetEnrichrCol4.input  = [infile]
		pTargetEnrichrCol4.errhow = 'terminate'
		PyPPL().start(pTargetEnrichrCol4).run()
	
	def testpTargetEnrichrCol5 (self):
		infile = path.join(tmpdir, 'targetEnrichrCol5.txt')		
		filestr = ''
		for reg, genes in relations.items():
			for g in genes:
				filestr += "%s\t%s\t%s\t%s\t%s\n" % (reg, g, random.choice(['+', '-']), random.choice(['+', '-']), random.choice(['+', '-']))
		ostr = ''
		if path.isfile(infile):
			with open(infile) as f:
				ostr = f.read()
		if not path.isfile(infile) or filestr != ostr:
			with open(infile, 'w') as f:
				f.write(filestr)

		pTargetEnrichrCol5        = pTargetEnrichr.copy()
		pTargetEnrichrCol5.input  = [infile]
		pTargetEnrichrCol5.errhow = 'terminate'
		PyPPL().start(pTargetEnrichrCol5).run()

	def testpEnrichr(self):
		infile = path.join(tmpdir, 'pEnrichr.txt')
		filestr = '\n'.join(genelist)
		ostr = ''
		if path.isfile(infile):
			with open(infile) as f:
				ostr = f.read()
		if not path.isfile(infile) or filestr != ostr:
			with open(infile, 'w') as f:
				f.write(filestr)

		pEnrichr.input  = [infile]
		pEnrichr.errhow = 'terminate'
		PyPPL().start(pEnrichr).run()

	def testpExpmat2Gct(self):
		infile  = path.join(tmpdir, 'pExpmat2Gct.txt')
		filestr = ''
		samples = ['Sample' + str(i+1) for i in list(range(10))]
		genes   = ['Gene' + str(i+1) for i in list(range(100))]
		filestr += "\t".join(samples) + '\n'
		for gene in genes:
			filestr += gene
			for sample in samples:
				filestr += '\t' + str(random.randint(1,100))
			filestr += '\n'

		ostr = ''
		if path.isfile(infile):
			with open(infile) as f:
				ostr = f.read()
		if not path.isfile(infile) or filestr != ostr:
			with open(infile, 'w') as f:
				f.write(filestr)
		pExpmat2GctCopy   = pExpmat2Gct.copy()
		pExpmat2GctCopy.input = [infile]
		PyPPL().start(pExpmat2GctCopy).run()	

	def testpSSGSEA(self):
		gctfile  = path.join(tmpdir, 'pSSGSEA.gct')
		gmtfile  = path.join(tmpdir, 'pSSGSEA.gmt')
		filestr  = ''
		samples  = ['Sample' + str(i+1) for i in list(range(1))]
		genes    = ['Gene' + str(i+1) for i in list(range(100))]
		pathws   = ['Pathway' + str(i+1) for i in list(range(20))]

		gctstr = "\t".join(samples) + '\n'
		for gene in genes:
			gctstr += gene
			for sample in samples:
				gctstr += '\t' + str(random.randint(1,100))
			gctstr += '\n'

		gmtstr = ''
		for pathw in pathws:
			gmtstr += pathw + '\t' + pathw + '\t' + '\t'.join(random.sample(genes, 20)) + '\n'

		gctostr = ''
		if path.isfile(gctfile):
			with open(gctfile) as f:
				gctostr = f.read()
		if not path.isfile(gctfile) or gctstr != gctostr:
			with open(gctfile, 'w') as f:
				f.write(gctstr)
		
		gmtostr = ''
		if path.isfile(gmtfile):
			with open(gmtfile) as f:
				gmtostr = f.read()
		if not path.isfile(gmtfile) or gmtstr != gmtostr:
			with open(gmtfile, 'w') as f:
				f.write(gmtstr)
		pExpmat2GctCopy2       = pExpmat2Gct.copy()
		pExpmat2GctCopy2.input = [gctfile]
		pSSGSEA.depends        = pExpmat2GctCopy2
		pSSGSEA.input          = lambda ch: ch.cbind(gmtfile)
		PyPPL().start(pExpmat2GctCopy2).run()


	def testpGSEA(self):
		gctfile  = path.join(tmpdir, 'pGSEA.gct')
		gmtfile  = path.join(tmpdir, 'pGSEA.gmt')
		sifile   = path.join(tmpdir, 'pGSEA.saminfo')
		filestr  = ''
		samples  = ['Sample' + str(i+1) for i in list(range(10))]
		genes    = ['Gene' + str(i+1) for i in list(range(100))]
		pathws   = ['Pathway' + str(i+1) for i in list(range(20))]
		sinfo    = {s:random.choice([0,1]) for s in samples}

		gctstr = "\t".join(samples) + '\n'
		for gene in genes:
			gctstr += gene
			for sample in samples:
				gctstr += '\t' + str(random.randint(1,100))
			gctstr += '\n'

		gmtstr = ''
		for pathw in pathws:
			gmtstr += pathw + '\t' + pathw + '\t' + '\t'.join(random.sample(genes, 20)) + '\n'

		sistr = 'Group\n'
		for sample, g in sinfo.items():
			sistr += '%s\tClass%s\n' % (sample, g)

		gctostr = ''
		if path.isfile(gctfile):
			with open(gctfile) as f:
				gctostr = f.read()
		if not path.isfile(gctfile) or gctstr != gctostr:
			with open(gctfile, 'w') as f:
				f.write(gctstr)
		
		gmtostr = ''
		if path.isfile(gmtfile):
			with open(gmtfile) as f:
				gmtostr = f.read()
		if not path.isfile(gmtfile) or gmtstr != gmtostr:
			with open(gmtfile, 'w') as f:
				f.write(gmtstr)
				
		siostr = ''
		if path.isfile(sifile):
			with open(sifile) as f:
				siostr = f.read()
		if not path.isfile(sifile) or sistr != siostr:
			with open(sifile, 'w') as f:
				f.write(sistr)

		pSampleinfo2Cls.input  = [sifile]
		pExpmat2GctCopy3       = pExpmat2Gct.copy()
		pExpmat2GctCopy3.input = [gctfile]
		pGSEA.depends          = pExpmat2GctCopy3, pSampleinfo2Cls
		pGSEA.args.nperm       = 1000
		pGSEA.args.nthread     = 4
		pGSEA.input            = lambda ch,   ch2: ch.cbind(ch2, gmtfile)
		PyPPL().start(pExpmat2GctCopy3, pSampleinfo2Cls).run()

	def testpGMT2Mat(self):
		gmtfile  = path.join(tmpdir, 'pGMT2Mat.gmt')
		genes    = ['Gene' + str(i+1) for i in list(range(100))]
		pathws   = ['Pathway' + str(i+1) for i in list(range(20))]

		gmtstr = ''
		for pathw in pathws:
			gmtstr += pathw + '\t' + pathw + '\t' + '\t'.join(random.sample(genes, 20)) + '\n'
		
		gmtostr = ''
		if path.isfile(gmtfile):
			with open(gmtfile) as f:
				gmtostr = f.read()
		if not path.isfile(gmtfile) or gmtstr != gmtostr:
			with open(gmtfile, 'w') as f:
				f.write(gmtstr)

		pGMT2Mat.input = [gmtfile]
		PyPPL().start(pGMT2Mat).run()
		matfile = pGMT2Mat.channel.outfile.get()
		with open(matfile) as f: lines = f.read().splitlines()
		self.assertEqual(len(lines), 101)
		self.assertEqual(lines[0], "\t" + "\t" .join(pathws))
		
if __name__ == '__main__':
	unittest.main(failfast=True)