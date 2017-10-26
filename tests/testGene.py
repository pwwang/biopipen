import unittest, addPath

from os import path
from pyppl import PyPPL, Channel, Box
from bioprocs import params
from bioprocs.gene import pGeneNameNorm, pGeneTss, pGenePromoters

params = params.toDict()
unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(x, y)

genes = [
	'Nsun3', 'Polrmt', 'Nlrx1', 'Sfxn5', 'Zc3h12c', 'Slc25a39', 'Arsg', 'Defb29', 'Ndufb6', 'Zfand1', 'Tmem77', '5730403B10Rik', 'RP23-195K8.6', 'Tlcd1', 'Psmc6', 'Slc30a6', 
	'LOC100047292', 'Lrrc40', 'Orc5l', 'Mpp7', 'Unc119b', 'Prkaca', 'Tcn2', 'Psmc3ip', 'Pcmtd2', 'Acaa1a', 'Lrrc1', '2810432D09Rik', 'Sephs2', 'Sac3d1', 'Tmlhe', 'LOC623451', 
	'Tsr2', 'Plekha7', 'Gys2', 'Arhgef12', 'Hibch', 'Lyrm2', 'Zbtb44', 'Entpd5', 'Rab11fip2', 'Lipt1', 'Intu', 'Anxa13', 'Klf12', 'Sat2', 'Gal3st2', 'Vamp8', 'Fkbpl', 'Aqp11', 
	'Trap1', 'Pmpcb', 'Tm7sf3', 'Rbm39', 'Bri3', 'Kdr', 'Zfp748', 'Nap1l1', 'Dhrs1', 'Lrrc56', 'Wdr20a', 'Stxbp2', 'Klf1', 'Ufc1', 'Ccdc16', '9230114K14Rik', 'Rwdd3', 
	'2610528K11Rik', 'Aco1', 'Cables1', 'LOC100047214', 'Yars2', 'Lypla1', 'Kalrn', 'Gyk', 'Zfp787', 'Zfp655', 'Rabepk', 'Zfp650', '4732466D17Rik', 'Exosc4', 'Wdr42a', 'Gphn', 
	'2610528J11Rik', '1110003E01Rik', 'Mdh1', '1200014M14Rik', 'AW209491', 'Mut', '1700123L14Rik', '2610036D13Rik', 'Cox15', 'Tmem30a', 'Nsmce4a', 'Tm2d2', 'Rhbdd3', 'Atxn2', 
	'Nfs1', '3110001I20Rik', 'BC038156', 'LOC100047782', '2410012H22Rik', 'Rilp', 'A230062G08Rik', 'Pttg1ip', 'Rab1', 'Afap1l1', 'Lyrm5', '2310026E23Rik', 'C330002I19Rik', 
	'Zfyve20', 'Poli', 'Tomm70a', 'Slc7a6os', 'Mat2b', '4932438A13Rik', 'Lrrc8a', 'Smo', 'Nupl2', 'Trpc2', 'Arsk', 'D630023B12Rik', 'Mtfr1', '5730414N17Rik', 'Scp2', 'Zrsr1', 
	'Nol7', 'C330018D20Rik', 'Ift122', 'LOC100046168', 'D730039F16Rik', 'Scyl1', '1700023B02Rik', '1700034H14Rik', 'Fbxo8', 'Paip1', 'Tmem186', 'Atpaf1', 'LOC100046254', 
	'LOC100047604', 'Coq10a', 'Fn3k', 'Sipa1l1', 'Slc25a16', 'Slc25a40', 'Rps6ka5', 'ABHD17AP4', 'ABHD17AP5'
]
genefile = path.join(params.tmpdir, 'genes.txt')
genestr  = ''
if path.isfile(genefile):
	with open(genefile) as f: genestr = f.read()

thestr  = ''
thestr += 'to be skipped\n'
thestr += '# skip again\n'
for i, g in enumerate(genes):
	thestr += '%s\t%s\n' % (i, g)

if thestr != genestr:
	with open(genefile, 'w') as f:
		f.write(thestr)

class TestGene (unittest.TestCase):
	def testGeneNameNorm(self):
		pGeneNameNorm.input     = [genefile]
		pGeneNameNorm.errhow    = 'terminate'
		pGeneNameNorm.args.skip = 1
		pGeneNameNorm.args.col  = 1
		PyPPL().start(pGeneNameNorm).run()

	def testGeneTss(self):
		pGeneTss.input = [genefile]
		pGeneTss.errhow    = 'terminate'
		pGeneTss.args.skip = 1
		pGeneTss.args.col  = 1
		PyPPL().start(pGeneTss).run()

	def testpGenePromoters(self):
		pGenePromoters.input = [genefile]
		pGenePromoters.errhow    = 'terminate'
		pGenePromoters.args.skip = 1
		pGenePromoters.args.col  = 1
		PyPPL().start(pGenePromoters).run()


if __name__ == '__main__':
	unittest.main(failfast=True)