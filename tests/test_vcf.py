import pytest
from pyppl import PyPPL
from bioprocs.vcf import pVcfFix
from . import assertInfile

def test_survival(rdata):
	pVcfFix1 = pVcfFix.copy()
	pVcfFix1.input = [rdata.get('vcf/badVcfToFix.vcf')]
	PyPPL().start(pVcfFix1).run()
	assertInfile(pVcfFix1.channel.outfile.get(), '##INFO=<ID=AF,Number=1,Type=String,Description="AF">')
	assertInfile(pVcfFix1.channel.outfile.get(), '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">')
	assertInfile(pVcfFix1.channel.outfile.get(), '##FILTER=<ID=LowEVS,Description="LOWEVS">')
	assertInfile(pVcfFix1.channel.outfile.get(), '##contig=<ID=chr99,length=999999999>')