"""Script for vcf.pVcfExplode"""
# pylint: disable=invalid-name,unused-import
from pathlib import Path
from diot import Diot
from cyvcf2 import VCF, Writer
from bioprocs.utils.parallel import Parallel

# pylint: disable=undefined-variable
infile = Path({{i.infile | quote}})
outdir = Path({{o.outdir | quote}})
nthread = {{args.nthread | repr}}
header_config = {{args.header | repr}}
# pylint: enable=undefined-variable

def one_contig(vcf, contig):
    """Do one contig
    Extract the records and save to new file"""

    outfile = outdir.parent.joinpath(infile.stem + '-' + contig + '.vcf')
    if header_config is False:
        # Needs test
        writer = Writer(outfile)
    elif header_config is True:
        writer = Writer(outfile, vcf)
    else:
        # Needs test
        headers = []
        for head in vcf.header_iter():
            if head['Type'] != 'contig' or head['ID'] == contig:
                headers.append(head)
        header_str = '\n'.join(str(head) for head in headers)
        writer = Writer.from_string(outfile, header_str)

    for rec in vcf(contig):
        writer.write_record(rec)

    writer.close()

vcf = VCF(infile)
para = Parallel(nthread, raiseExc=True)
para.run(one_contig, [(seq, ) for seq in vcf.seqnames]
