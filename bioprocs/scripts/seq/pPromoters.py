import sys

from os import path
from pyppl import Box
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
from bioprocs.utils import log2pyppl

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
region   = {{args.region | repr}}
notfound = {{args.notfound | quote}}
inopts   = {{args.inopts | repr}}
outopts  = {{args.outopts | repr}}
refgene  = {{args.refgene | quote}}
genecol  = {{args.inopts.genecol | repr}}
ocnames  = {{args.outopts.cnames | repr}}
del inopts['genecol']
del outopts['cnames']
if region.down is None:
	region.down = region.up

# get all genes' TSS and strand
reader = TsvReader(refgene, cnames = False, delimit = '"')
genes  = {r[1]:r[0].split("\t")[:7] for r in reader}

reader = TsvReader(infile, **inopts)
writer = TsvWriter(outfile, **outopts)
if ocnames:
	writer.cnames = ['CHROM', 'START', 'END', 'NAME', 'SCORE', 'STRAND']
	writer.writeHead()
for r in reader:
	gene = r[genecol]
	if gene not in genes:
		msg = 'Gene does not exist: {}'.format(gene)
		if notfound == 'error':
			raise ValueError(msg)
		else:
			log2pyppl('Gene does not exist: {msg}', 'warning')
			continue
	chrom, _, _, start, end, _, strand = genes[gene]
	start, end = int(start), int(end)
	if strand == '-':
		record = [
			chrom, 
			min(start, end - region.down) if region.withbody else end - region.down, 
			end + region.up, 
			gene, 0, strand
		]
	else:
		record = [
			chrom, 
			start - region.up, 
			max(end, start + region.down) if region.withbody else start + region.down, 
			gene, 0, strand
		]
	writer.write(record)
writer.close()

