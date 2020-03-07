"""Putative transcription factor binding sites analysis"""
from diot import Diot
from .utils import fs2name
from . import params, proc_factory

# pylint: disable=invalid-name

pMotifScan = proc_factory(
    desc='Scan motif along the given sequences.',
    config=Diot(annotate="""
    @input:
        `tffile:file`: The infile containing TF name and motif name.
            - If only one column is give, will be used as both TF and motif name
            - If there are 2+ columns, 1st column will be motif name, 2nd column will be TF name
        `sfile:file`: The sequence file
    @output:
        `outdir:file`: The output dir
    @args:
        `tool` : The tool used to scan the motif.
        `meme` : The path of MEME's fimo.
        fimo: Alias for `args.meme`
        `tfmotifs`: The motif database in MEME format.
        `cutoff` : The cutoff for the hits.
            - Could be a float that denote the pvalue threshold
            - Or a dict with keys `by` indicating pvalue (`p`) or qvalue (`q`) to filter, and `value` indicating the cutoff value.
        `cleanmname`: Whether to clean motif name.
        `ucsclink`: The ucsc link template.
        `nthread` : Number of threads used to scan, only available when you have multiple m-ids.
        `params` : Other parameters for `fimo`
    @requires:
        [`fimo` from MEME Suite](http://meme-suite.org/tools/fimo)
	"""),
    lang=params.python.value,
    input="tffile:file, sfile:file",
    output=[
        "outfile:file:{{i.sfile | fn}}-{{i.tffile | fn}}.fimo"
        "/{{i.sfile | fn}}-{{i.tffile | fn}}.bed",
        "outdir:dir:{{i.sfile | fn}}-{{i.tffile | fn}}.fimo"
    ],
    args=Diot(
        tool='meme',
        meme=params.fimo.value,
        fimo=params.fimo.value,
        params=Diot(),
        tfmotifs=params.tfmotifs.value,
        cutoff=1e-6,
        ucsclink='https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={}',
        nthread=1,
    )
)

pMotifSimilarity = proc_factory(
    desc='Compare the similarity between motifs.',
    config=Diot(annotate="""
    @name:
        pMotifSimilarity
    @description:
        Compare the similarity between motifs.
    @input:
        `mfile1:file`: The motif database in meme format
        `mfile2:file`: Another motif database in meme format.
            - If not provided, `mfile1` will be used.
    @output:
        `outfile:file`: The output file from `tomtom`
        `outdir:dir`  : The output directory from `tomtom`
    @args:
        `qval`   : The qvalue cutoff. Default: `0.5`
        `tomtom` : The path to `tomtom`
        `nthread`: Number of threads to use. Default: `1`
        `params` : parameters for `tomtom`:
            - `xalph`       : `True`
            - `no-ssc`      : `True`
            - `dist`        : `pearson`,
            - `min-overlap` : `1`,
            - `motif-pseudo`: `.1`,
            - `verbosity`   : `4`
	"""),
    lang=params.python.value,
    input='mfile1:file, mfile2:file',
    output=[
        'outfile:file:{{i.mfile1, i.mfile2, fn2 | *odfn \
                                                | path.join: "tomtom.tsv"}}',
        'outdir:dir:{{i.mfile1, i.mfile2, fn2 | *odfn}}'
    ],
    args=Diot(
        qval=0.5,
        tomtom=params.tomtom.value,
        params=Diot({
            'xalph': True,
            'no-ssc': True,
            'dist': 'pearson',
            'min-overlap': 1,
            'motif-pseudo': .1,
            'verbosity': 4
        }),
        nthread=1
    ),
    envs=Diot(
        path=__import__('os').path,
        odfn=lambda f1, f2, fn2: (
            fn2(f1),
            fn2(f1) + "-" + fn2(f2)
        )[int(bool(f2))] + ".tomtom"
    )
)

pMotifMerge = proc_factory(
    desc='Merge motif files in MEME format',
    config=Diot(annotate="""
    @name:
        pMotifMerge
	"""),
    input='infiles:files',
    output='outfile:file:{{i.infiles | fs2name}}.meme.txt',
    lang=params.python.value,
    envs=Diot(fs2name=fs2name)
)

pMotifFilter = proc_factory(
    desc='Filter motifs from a meme file.',
    config=Diot(annotate="""
    @name:
        pMotifFilter
    @description:
        Filter motifs from a meme file.
    @input:
        `infile:file`: The motif database in MEME format
    @output:
        `outfile:file`: the output file.
    @args:
        `filter`: A string of lambda function to filter the motifs.
            - Possible attributes for a motif are:
            - `E, URL, alength, altname, matrix, mtype, name, nsites, w`
            - Motifs will be filtered out when this function returns `True`.
	"""),
    input='infile:file',
    output='outfile:file:{{i.infile | bn}}',
    lang=params.python.value,
    args=Diot(filter=None)
)

pAtSnp = proc_factory(
    desc='Scan motifs on Snps to detect binding affinity changes.',
    config=Diot(annotate="""
    @name:
        pAtSnp
    @description:
        Scan motifs on Snps to detect binding affinity changes.
    @input:
        `tffile:file`:  The tf-motif file with 1st column the motif and 2nd the tf
        `snpfile:file`: The snp file.
            - Could be a bed file with first 6 columns the bed6 format and 7th the reference allele, and 8th the alternative alleles.
            - Alleles including reference allele should be seperated by `,`
            - It also could be a snp file required by `atSNP` package.
            - See: https://rdrr.io/github/kinsigne/atSNP_modified/man/LoadSNPData.html
    @output:
        `outfile:file`: The output file
        `outdir:dir`  : The output directory containing the output file and plots.
    @args:
        `tfmotifs`: The motif database. Defaut: `params.tfmotifs`
        `genome`  : The reference genome to get the sequences. Default: `params.genome`
        `fdr`     : Do fdr or not. Or could be the p-value adjustment method. Default: `True` (using `BH` method)
        `pval`    : The pvalue cutoff for output and plot.
        `plot`    : Do plot or not. Default: `True`
        `nthread` : How many threads to use. Default: `1`
        `depvars` : The device parameters for plotting. Default: `Diot(res = 300, width = 2000, height = 2000)`
    @requires:
        `r-atSNP`
	"""),
    lang=params.Rscript.value,
    input='tffile:file, snpfile:file',
    output=[
        'outfile:file:{{i.tffile | fn2}}-{{i.snpfile | fn2}}.atsnp'
        '/{{i.tffile | fn2}}-{{i.snpfile | fn2}}.atsnp.txt',
        'outdir:dir:{{i.tffile | fn2}}-{{i.snpfile | fn2}}.atsnp'
    ],
    args=Diot(
        tfmotifs=params.tfmotifs.value,
        genome=params.genome.value,
        fdr=True,
        pval=0.05,
        plot=True,
        nthread=1,
        devpars=Diot(res=300, width=2000, height=2000),
    )
)
