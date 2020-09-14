"""Script for gsea.pExprMat2GCT"""
# pylint: disable=undefined-variable,unused-import,invalid-name
from bioprocs.utils import logger
from bioprocs.utils.tsvio2 import TsvReader

expfile = {{i.expfile | quote}}
outfile = {{o.outfile | quote}}
reader = TsvReader(expfile)
ngenes = 0
nsamples = 0
rnames = {}
for row in reader:
    if nsamples == 0:
        nsamples = len(row) - 1
    ngenes += 1
    rnames[row[0]] = 1
reader.rewind()
samples = reader.cnames[-nsamples:]

with open(outfile, "w") as fout:
    fout.write("#1.2\n")
    fout.write("%s\t%s\n" % (ngenes, len(samples)))
    fout.write("NAME\tDescription\t%s\n" % ('\t'.join(samples)))
    for row in reader:
        rnames[row[0]] -= 1
        if rnames[row[0]] < 0:
            logger.warning('Duplicate gene found: %s, ignored!' % row[0])
            continue
        fout.write("\t".join([row[0]] + list(row)) + "\n")
