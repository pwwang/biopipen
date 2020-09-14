"""Script for seq.pPromoters"""

# pylint: disable=invalid-name,unused-import,not-a-mapping,undefined-variable

from diot import Diot
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
from bioprocs.utils import log2pyppl

infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
region = {{args.region | repr}}
notfound = {{args.notfound | quote}}
inopts = {{args.inopts | repr}}
refgene = {{args.refgene | quote}}
genecol = {{args.genecol | repr}}
base = {{args.base | repr}}
if region.down is None:
    region.down = region.up

# get all genes' TSS and strand
reader = TsvReader(refgene, cnames=False, delimit='"')
genes = {r[1]: r[0].split("\t")[:7] for r in reader}
log2pyppl(f"{len(genes)} loaded from database ...")

reader = TsvReader(infile, **inopts)
writer = TsvWriter(outfile)
for r in reader:
    gene = r[genecol]
    if gene not in genes:
        msg = 'Gene does not exist: {}'.format(gene)
        if notfound == 'error':
            raise ValueError(msg)

        log2pyppl('Gene does not exist: {msg}', 'warning')
        continue
    chrom, _, _, start, end, _, strand = genes[gene]
    start, end = int(start), int(end)
    if strand == '-':
        record = [
            chrom,
            max(base,
                min(start, end - region.down - 1 + base)
                if region.withbody
                else end - region.down - 1 + base),
            end + region.up,
            gene, 0, strand
        ]
    else:
        record = [
            chrom,
            max(base, start - region.up - 1 + base),
            (max(end, start + region.down)
             if region.withbody else start + region.down),
            gene, 0, strand
        ]
    writer.write(record)
writer.close()
