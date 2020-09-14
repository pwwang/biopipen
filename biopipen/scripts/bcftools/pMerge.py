"""Script for bcftools.pMerge"""
# pylint: disable=not-a-mapping,invalid-name,expression-not-assigned
from pathlib import Path
from diot import Diot # pylint: disable=unused-import
from bioprocs.utils import shell2 as shell, gztype, logger
from bioprocs.utils.reference import vcf_index
from bioprocs.utils.parallel import Parallel

# pylint: disable=undefined-variable, unsubscriptable-object
infiles = {{i.infiles | repr}}
outfile = {{o.outfile | quote}}
nthread = {{args.nthread | repr}}
bcftools = {{args.bcftools | quote}}
tabix = {{args.tabix | quote}}
gz = {{args.gz | bool}}
params = {{args.params | repr}}
bychr = {{args.bychr | bool}}
# pylint: enable=undefined-variable

# check if infiles are bgipped and indexed
infile0 = infiles[0]
if gztype(infile0) != 'bgzip':
    raise ValueError(
        '`bcftools merge` requires vcf files to be bgzipped and indexed.'
    )
for infile in infiles:
    vcf_index(infile, tabix)

shell.load_config(bcftools=bcftools, tabix=tabix)

# write files to a file and pass it by -l
# to prevent the command from being too long
vcflistfile = Path(outfile).parent.joinpath('vcflist.txt')
vcflistfile.write_text('\n'.join(str(infile) for infile in infiles))
if bychr:
    # get list of contigs
    # all infiles supposingly to have the same set of contigs
    contigs = shell.tabix(l=infile0).splitlines()

    merging_dir = Path(outfile).parent.joinpath('__merging__')
    merging_dir.mkdir()

    def merge_one_contig(contig):
        logger.info('Merge %s ...', contig)
        contig_vcf = merging_dir / f'{contig}.vcf.gz'
        params_contig = params.copy()
        params_contig.l = vcflistfile
        params_contig.o = contig_vcf
        params_contig.r = contig
        params_contig.O = 'z'
        shell.bcftools.merge(**params_contig).fg()
        vcf_index(contig_vcf)

    para = Parallel(nthread, 'process', True)
    para.run(merge_one_contig, [(contig,) for contig in contigs])

    if 'missing-to-ref' in params:
        del params['missing-to-ref']
    if '0' in params:
        del params['0']

    logger.info('Concat contigs ...')
    params.naive = True
    params.o = outfile
    params.O = 'z' if gz else 'v'
    params.threads = nthread
    params._ = list(merging_dir.glob('*.vcf.gz'))

    shell.bcftools.concat(**params).fg()

else:
    # bcftools merge [OPTIONS] A.vcf.gz B.vcf.gz […]

    params.l = vcflistfile
    params.threads = nthread
    params.O = 'z' if gz else 'v'
    params.o = outfile

    shell.bcftools.merge(**params).fg()
