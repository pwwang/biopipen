"""Script for bed.pBed2Wiggle"""
# pylint: disable=invalid-name
from bioprocs.utils.tsvio2 import TsvReader

infile = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
span = {{args.span | int}}
step = {{args.step | int}}
inbase = {{args.inbase | int}}
addchr = {{args.addchr | bool}}
wigtobigwig = {{args.wigtobigwig | quote}}
bigwig = {{args.bigwig | bool}}
gsize = {{args.gsize | quote}}

reader = TsvReader(infile, cnames = False)

with open(outfile, 'w') as fout:
    for r in reader:
        fout.write(' '.join([
            "fixedStep",
            "chrom={}".format(r[0] if not addchr or "chr" in r[0]
                              else "chr" + r[0]),
            "start={}".format(r[1] if inbase == 1 else int(r[1]) + 1),
            "step={}".format(step),
            "span={}".format(span)
        ]) + "\n")

        signals = (sig.strip() for sig in r[-1].split(","))
        fout.write("\n".join(signals) + "\n")

reader.close()

if not bigwig:
    exit(0)

from bioprocs.utils import shell2 as shell
shell.mv(outfile, outfile + '.wig')
shell.load_config(wigtobigwig=wigtobigwig)
shell.wigtobigwig(outfile + '.wig', gsize, outfile, _sep="=", _prefix="-").fg
