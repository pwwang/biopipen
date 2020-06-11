"""Script for bcftools.pMerge"""
# pylint: disable=not-a-mapping,invalid-name,expression-not-assigned
from pathlib import Path
from diot import Diot # pylint: disable=unused-import
from bioprocs.utils import shell2 as shell, gztype
from bioprocs.utils.reference import vcf_index

# pylint: disable=undefined-variable, unsubscriptable-object
infiles = {{i.infiles | repr}}
outfile = {{o.outfile | quote}}
nthread = {{args.nthread | repr}}
bcftools = {{args.bcftools | quote}}
gz = {{args.gz | bool}}
params = {{args.params | repr}}
# pylint: enable=undefined-variable

# check if infiles are bgipped and indexed
infile0 = infiles[0]
if gztype(infile0) != 'bgzip':
    raise ValueError(
        '`bcftools merge` requires vcf files to be bgzipped and indexed.'
    )
for infile in infiles:
    vcf_index(infile)

shell.load_config(bcftools=bcftools)

# bcftools merge [OPTIONS] A.vcf.gz B.vcf.gz […]

if len(infiles) < 20:
    params._ = infiles
else:
    # write files to a file and pass it by -l
    # to prevent the command from being too long
    vcflistfile = Path(outfile).parent.joinpath('vcflist.txt')
    vcflistfile.write_text('\n'.join(str(infile) for infile in infiles))
    params.l = vcflistfile
params.threads = nthread
params.O = 'z' if gz else 'v'
params.o = outfile

shell.bcftools.merge(**params).fg
