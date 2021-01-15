"""Commands of bcftools"""
from modkit import modkit
from diot import Diot
from .utils import fs2name
from . import opts, module_postinit

def p_query():
    """
    @input
        infile: The input VCF file. Currently only single VCF file is supported.
    @output
        outfile: The output file
    @args:
        bcftools (str): Path to bcftools.
        params (Diot): Other parameters for `bcftools query`.
            - See: https://samtools.github.io/bcftools/bcftools.html#query
    """
    return Diot(
        desc=("Extracts fields from VCF or BCF files and "
              "outputs them in user-defined format"),
        lang=opts.python,
        input='infile:file',
        output='outfile:file:{{i.infile | stem | stem}}.query.txt',
        args=Diot(bcftools=opts.bcftools,
                  params=Diot())
    )

def p_view():
    """
    @input
        infile: The input VCF file. Currently only single VCF file is supported.
        samfile: The sample name file to extract samples from infile
            - It also can be samples directly, separated by comma
            - If this is provided, `args.params.s` and `args.params.S`
              will be ignored
            - See https://samtools.github.io/bcftools/bcftools.html#view
    @output
        outfile: The output file
    @args:
        tabix    (str) : Path to tabix.
        bcftools (str) : Path to bcftools.
        gz       (bool): Whether output gzipped vcf file or not.
        nthread  (int) : Number of threads to use.
        params   (Diot) : Other parameters for `bcftools view`.
            - See: https://samtools.github.io/bcftools/bcftools.html#view
    """
    return Diot(
        desc=("View, subset and filter VCF or BCF files by position and "
              "filtering expression."),
        input='infile:file, samfile:var',
        output='''outfile:file:{{i.infile | stem2}}.vcf{{args.gz | ?
                                                                 | =:".gz"
                                                                 | !:""}}''',
        lang=opts.python,
        args=Diot(gz=False,
                  nthread=1,
                  tabix=opts.tabix,
                  bcftools=opts.bcftools,
                  params=Diot())
    )

def p_reheader():
    """
    @input:
        infile: The input VCF file
        hfile: The new header file
        samfile: The new sample name file
            - It also can be samples directly, separated by comma
            - If this is provided, `args.params.samples` will be ignored
            - See https://samtools.github.io/bcftools/bcftools.html#reheader
    @output:
        outfile: The output VCF file with new header.
    @args:
        bcftools (str): Path to bcftools.
        nthread  (int): Number of threads to use.
        params   (Diot): Other parameters for `bcftools view`.
            - See https://samtools.github.io/bcftools/bcftools.html#reheader
    """
    return Diot(
        desc="Modify header of VCF/BCF files, change sample names",
        input='infile:file, hfile:file, samfile:var',
        output='outfile:file:{{i.infile | bn}}',
        lang=opts.python,
        args=Diot(bcftools=opts.bcftools,
                  nthread=1,
                  params=Diot())
    )

def p_filter():
    """
    @input:
        infile: The input VCF file
    @output:
        outfile: The output filtered VCF file
        statfile: Statistics of variants for each FILTER.
    @args:
        bcftools (str) : Path to bcftools
        nthread  (int) : Number of threads to use.
        gz       (bool): Whether output gzipped vcf file or not.
            - Overwrite `args.params.O`
        keep     (bool): Whether should we keep the excluded variants or not.
            - If not, args about filter names will not work.
        include  (str|list|dict): include only sites for which EXPRESSION
                                  is true.
            - See: https://samtools.github.io/bcftools/bcftools.html#expressions
            - Use `args.params.include` or `args.params.exclude` only if you
              just have one filter.
            - If provided, `args.params.include/exclude` will be ignored.
            - If `str`/`list` used, The filter names will be `Filter%d`
            - A dict is used when keys are filter names and values are
              expressions
        exclude  (str|list|dict): exclude sites for which EXPRESSION is true.
            - See also `args.include`
        params   (Diot)          : Other parameters for `bcftools filter`
            - See: https://samtools.github.io/bcftools/bcftools.html#filter
    """
    return Diot(
        desc="Apply fixed-threshold filters to VCF files",
        lang=opts.python,
        input='infile:file',
        output=['outfile:file:{{i.infile | bn}}',
                'statfile:file:{{i.infile | stem | stem}}.filterstats.txt'],
        args=Diot(bcftools=opts.bcftools,
                  nthread=1,
                  gz=False,
                  stat=False,
                  keep=True,
                  include=None,
                  exclude=None,
                  # accumulating Filter names instead of repolacing
                  params=Diot(mode='+'))
    )

def p_annotate():
    """
    @input:
        infile: The input VCF file
    @output:
        outfile: The annotated output VCF file
    @args:
        tabix   (str)  : Path to tabix, used to index `annfile`
        bcftools (str) : Path to bcftools
        annfile  (file): The annotation file.
            - See: https://samtools.github.io/bcftools/bcftools.html#annotate
        nthread (int)     : Number of threads to use.
        cols    (str|list): Overwrite `-c/--columns`.
        header  (str|list): headers to be added.
        params  (Diot)     : Other parameters for `bcftools annotate`
    """
    return Diot(
        desc='Add or remove annotations from VCF files',
        lang=opts.python,
        input='infile:file',
        output='outfile:file:{{i.infile | bn}}',
        args=Diot(tabix=opts.tabix,
                  bcftools=opts.bcftools,
                  nthread=1,
                  annfile='',
                  cols=[],
                  header=[],
                  params=Diot())
    )

def p_concat():
    """
    @input:
        infiles: The input vcf files
    @output:
        outfile: The output merged vcf file
    @args:
        nthread  (int) : The number of threads to use
        bcftools (path): The path to bcftools
        tabix    (path): The path to tabix, used to index vcf files.
        params   (Diot) : Other parameters for `bcftools concat`
        gz       (bool): Whether output gzipped vcf or not.
    """
    return Diot(
        desc=('Concatenate or combine VCF/BCF files with same samples '
              'in the same order.'),
        input='infiles:files',
        output='''outfile:file:{{i.infiles | [0]
                                           | stem2
                                           | @append: "_etc.vcf"
                                }}{{ args.gz | ?
                                             | =:".gz"
                                             | !:""}}''',
        lang=opts.python,
        args=Diot(nthread=1,
                  bcftools=opts.bcftools,
                  tabix=opts.tabix,
                  params=Diot(),
                  gz=False)
    )

class PMerge:
    """
    @input:
        infiles: The input vcf files
    @output:
        outfile: The merged vcf file
    @args:
        nthread (int): The number of threads to use
        params (Diot) : Other parameters for `bcftools merge`
        bcftools (str): Path to bcftools
        gz (bool): Whether output gzipped vcf
        tabix (str): Path to tabix
        bychr (bool): Merge by chromosome and then concatenate
    """
    desc = ('Merge multiple VCF/BCF files from non-overlapping sample sets '
            'to create one multi-sample file using `bcftools merge`.')
    lang = opts.python
    input = 'infiles:files'
    output = ('outfile:file:{{i.infiles | fs2name}}.vcf'
              '{{args.gz | ?! :"" | ?= :".gz"}}')
    envs = Diot(fs2name=fs2name)
    args = Diot(nthread=1,
                bcftools=opts.bcftools,
                gz=False,
                tabix=opts.tabix,
                bychr=False,
                params=Diot())

modkit.alias(pBcftoolsMerge='pMerge',
             pBcftoolsAnnotate='pAnnotate',
             pBcftoolsConcat='pConcat',
             pBcftoolsFilter='pFilter',
             pBcftoolsQuery='pQuery',
             pBcftoolsReheader='pReheader',
             pBcftoolsView='pView')
modkit.postinit(module_postinit)
