"""Gene related processes"""
from diot import Diot
from . import opts, proc_factory
from .seq import pPromoters

# pylint: disable=invalid-name

pGenePromoters = pPromoters.copy()

pGeneNameNorm = proc_factory(
    desc='Normalize gene names using MyGeneinfo.',
    config=Diot(annotate="""
    @name:
        pGeneNameNorm
    @description:
        Normalize gene names using MyGeneinfo.
    @input:
        `infile:file`: The input file
    @output:
        `outfile:file`: The output file
    @args:
        `inopts`: options for reading input file.
        `outopts`: options for writing output file.
            - `query` : Output the original query column? Default: `False`
            - `cnames`: Output headers? Default: `True`
        `notfound`: What if a symbol is not found. Default: ignore
            - skip  : skip the record(don't write it to output file)
            - ignore: use the original name;
            - error : report error
        `genecol` : the column index containing the gene names
        `frm`     : the original format. Default: 'symbol, alias'
        `to`      : the output gene name format. Default: 'symbol'
        `genome`  : the genome. Default: 'hg19'
        `cachedir`: The cache directory
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | bn}}',
    errhow='retry',
    lang=opts.python,
    args=Diot(
        notfound='ignore',
        inopts=Diot(skip=0, comment='#', delimit='\t'),
        outopts=Diot(delimit='\t', cnames=True, query=False),
        genecol='',
        frm='symbol, alias',
        to='symbol',
        genome=opts.genome,
        cachedir=opts.cachedir,
    )
)

pIPI = proc_factory(
    desc='Convert gene symbol to IPI protein accession and vice versa.',
    config=Diot(annotate="""
    @name:
        pIPI
    @description:
        Convert gene symbol to IPI protein accession and vice versa.
        One gene symbol could map to multiple IPIs, which will be separated by pipe (|)
    @input:
        `infile:file` : The input file
    @output:
        `outfile:file`: The output file
    @args:
        `notfound`: What if a record is not found: Default: `ignore`
            - `skip`  : skip the record(don't write it to output file)
            - `ignore`: use the original name;
            - `error` : report error
        `genecol`: The column index containing the gene/protein record
        `ipidb`: The IPI xref database (see http://ftp.ebi.ac.uk/pub/databases/IPI/last_release/current/).
        `fromipi`: Whether the input is IPI or genes
        `inopts`: The options for input file
        `outopts`: The options for output file
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | bn}}',
    errhow='retry',
    args=Diot(
        notfound='ignore',
        inopts=Diot(skip=0, comment='#', delimit='\t'),
        outopts=Diot(delimit='\t',
                     headDelimit='\t',
                     headPrefix='',
                     headTransform=None,
                     head=True,
                     query=False),
        genecol=None,
        fromipi=True,
        ipidb=opts.ipidb,
    ),
    lang=opts.python,
)

pGeneTss = proc_factory(
    desc='Get gene TSS in BED format',
    config=Diot(annotate="""
    @name:
        pGeneTss
    @description:
        Get gene TSS in BED format.
    @input:
        `infile:file`: The input file containing genes
    @output:
        `outfile:file`: The output BED file
    @args:
        `notfound`: What if the gene is not found. Default: skip.
            - error: report error
        `header`: Whether the input file contains header. Default: False
        `skip`: Skip N lines of input file. Default: 0
            - This has highest priority of header and comment
        `comment`: The comment line start sign. Default: #
        `delimit`: The delimit of input file if it has multiple column. Default: `\\t`
        `col`: The column index contains the genes. Default: 0
        `frm`: The format of the genes. Default: `symbol, alias`
        `tmpdir`: The tmpdir used to store mygene cache files.
        `genome`: In which genome to fetch the coordinates. Default: hg19
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}-tss.bedx',
    errhow='retry',
    args=Diot(
        notfound='skip',  # error,
        genecol='',
        inopts=Diot(skip=0, comment='#', delimit='\t'),
        outopts=Diot(delimit='\t',
                     headDelimit='\t',
                     headPrefix='',
                     headTransform=None,
                     head=False,
                     query=False,
                     ftype='bed'),
        frm='symbol, alias',
        cachedir=opts.cachedir,
        genome=opts.genome,
    ),
    #pGeneTss.envs.genenorm  = genenorm.py
    #pGeneTss.envs.writeBedx = write.bedx.py
    lang=opts.python,
)

pGeneBody = proc_factory(
    desc='Get gene body in BED format',
    config=Diot(annotate="""
    @name:
        pGeneBody
    @description:
        Get gene body region in BED format
    @input:
        `infile:file`: The input file containing genes
    @output:
        `outfile:file`: The gene body region
    @args:
        `notfound`: What if a gene is not found when transfer the gene names to gene symbols
            - error: report error
            - skip (default): skip it
        `inmeta`:   The metadata for input file, mainly to indicate where the GENE column is.
        `inopts`:   Input options for reading input file.
            - skip: number of lines to skip. Default: 0
            - comment: the starting string for comment lines. Default: #
            - delimit: The delimit for the input file. Default: '\\t'
        frm: The gene name format in the input file. Default: 'symbol, alias'
        tmpdir: The tmpdir to cache the gene name conversion.
        genome: The genome used to do the conversion.
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}-body.bed',
    args=Diot(
        inopts=Diot(cnames=False),
        notfound='skip',  # error,
        genecol='',
        refgene=opts.refgene,
    ),
    lang=opts.python,
)
