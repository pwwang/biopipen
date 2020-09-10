from cyvcf2 import VCF, Writer
from bioprocs.utils import shell2 as shell, logger

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
gz      = {{args.gz | repr}}
if gz:
    outfile = outfile[:-3]
helper  = {{args.helper | repr}}
if not isinstance(helper, list):
    helper = [helper]
helper  = [line for line in helper if line]
exec('\n'.join(helper), globals())
vcfops    = {{args.vcf}}
recordops = {{args.record}}

vcf = VCF(infile)
if callable(vcfops):
    vcfops(vcf)

writer = Writer(outfile, vcf)
for var in vcf:
    recordops(var)
    try:
        writer.write_record(var)
    except Exception as ex:
        logger.error(str(ex) + '. Record: ' + str(var))
writer.close()
vcf.close()

if gz:
    shell.bgzip(outfile)
