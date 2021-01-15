"""Processes to process files with SNPs"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pRs2Bed = proc_factory(
    desc='Find coordinates for SNPs in BED format.',
    config=Diot(annotate="""
    @name:
        pRs2Bed
    @description:
        Find coordinates for SNPs in BED format.
    @input:
        `snpfile:file`: the snp file, each snp per line
    @output:
        `outfile:file`: the result file, could be a 3 or 6-col bed file.
    @args:
        `tool`    : The tool used to retrieve the coordinates. Default: `cruzdb`
            - Or retrieve from a dbsnp VCF file (dbsnp or local)
        `dbsnp`   : The dbsnp vcf file
        `notfound`: What to do if the snp is not found. Default: skip
        `inopts`  : The input options for input file. Default: `Diot(delimit = '\t', skip = 0, comment = '#')`
        `snpcol`  : The column where the snp is. Could be index (0-based) or the column name with `args.inopts.cnames = True`
        `ncol`    : How many columns to output. Possible values: 6, 8, 9. Default: `8`
            - `6`: The BED6 format
            - `8`: BED6 + `ref` + `alt`
            - `9`: BED6 + `ref` + `alt` + `allele frequences from 1000 genome project`
            - `args.ncol = 9` only avaliable when `args.tool = "dbsnp"`
        `vcftools`: The path to vcftools to extract snp from `args.dbsnp`
        `genome`  : The genome to locate the right database from `cruzdb`
        `dbsnpver`: The version of dbsnp from `cruzdb`
        `sortby`  : Sort the output file by coordinates (coord) or name. Default: `coord`
            - `False`: don't sort.
    """),
    input="snpfile:file",
    output="outfile:file:{{i.snpfile | fn}}.bed",
    lang=opts.python,
    args=Diot(
        tool='cruzdb',  # or local/dbsnp,
        notfound='skip',  # error,
        inopts=Diot(delimit='\t', skip=0, comment='#'),
        ncol=6,
        snpcol='',
        dbsnp=opts.dbsnp_all,
        vcftools=opts.vcftools,
        sortby='coord',  # name/False,
        genome=opts.genome,
        dbsnpver=opts.dbsnpver,
    )
)

pCoord2SnpBedx = proc_factory(
    desc='Find snps with coordinates.',
    config=Diot(annotate="""
    @name:
        pCoord2SnpBedx
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}-snps.bed',
    lang=opts.python,
    args=Diot(
        genome=opts.genome,
        dbsnpver=opts.dbsnpver,
        inmeta=['CHR', 'START', 'END'],
        inopts=Diot(delimit='\t', skip=0, comment='#'),
        xcols=[
            'refUCSC', 'alleles', 'alleleFreqs', 'alleleFreqCount'
        ]
    )
)
pCoord2SnpBedx.errhow = 'retry'

pSnp2Avinput = proc_factory(
    desc='Convert SNP list to avinput to ANNOVAR.',
    config=Diot(annotate="""
    @name:
        pSnp2Avinput
    @description:
        Convert SNP list to avinput to ANNOVAR.
    @input:
        `snpfile:file`: the snp file, each snp per line
    @output:
        `outfile:file`: the result avinput file
    @args:
        `genome`: default: hg19
        `snpver`: default: snp147
    @requires:
        [`python-cruzdb`](https://github.com/brentp/cruzdb)
    """),
    input="snpfile:file",
    output="outfile:file:{{i.snpfile | fn}}.avinput",
    lang=opts.python,
    args=Diot(
        genome=opts.genome,
        dbsnpver=opts.dbsnpver,
        notfound='skip',  # error,
        header=False,
        skip=0,
        comment='#',
        delimit='\t',
        col=0,
    )
)
pSnp2Avinput.errhow = 'retry'
