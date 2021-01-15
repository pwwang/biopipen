"""TCGA data analysis"""
from diot import Diot
from .utils import fs2name
from . import opts, proc_factory

# pylint: disable=invalid-name

pTCGADownload = proc_factory(
    desc='Download data with gdc-client and a menifest file.',
    config=Diot(echo_jobs=0, echo_types='stderr', annotate="""
    @name:
        pTCGADownload
    @description:
        Download TCGA use `gdc-client` and a manifest file
    @input:
        `manifile:file`: the manifest file
    @output:
        `outdir:file`:   the directory containing downloaded file
    @args:
        `params`    : other params for `gdc-client download`, default: `{'no-file-md5sum': True}`
        `gdc_client`: the executable file of `gdc-client`,    default: "gdc-client"
        `nthread`   : Number of threads to use. Default     : `1`
        `token`     : The token file if needed.
    """),
    input="manifile:file",
    output="outdir:dir:{{i.manifile | fn}}",
    lang=opts.python,
    args=Diot(
        params=Diot({'no-file-md5sum': True}),
        nthread=1,
        token=None,
        gdc_client=opts.gdc_client,
    )
)

pSample2SubmitterID = proc_factory(
    desc=('convert TCGA sample names with submitter id with '
          'metadata and sample containing folder'),
    config=Diot(annotate="""
    @name:
        pSample2SubmitterID
    @description:
        convert TCGA sample names with submitter id with metadata and sample containing folder
    @input:
        `indir:file`:    the directory containing the samples
        `mdfile:file`: the metadata file
    @output:
        `outdir:file`: the directory containing submitter-id named files
    @args:
        `method`: How the deal with the files.
            - We can also do `copy`
        `nthread`: Number threads to use.
        len: The length of the SubmitterID to keep.
            - The submitter id can be as long as "TCGA-78-7150-01A-21D-2035-01"
    """),
    input="indir:file, mdfile:file",
    output="outdir:dir:{{i.indir | fn}}",
    lang=opts.python,
    args=Diot(
        method='symlink',  # or copy,
        nthread=1,
        len=14 #TCGA-67-3771-10
    )
)

pGtFiles2Mat = proc_factory(
    desc='Convert genotype files to a matrix.',
    config=Diot(annotate="""
    @name:
        pGtFiles2Mat
    @description:
        Convert TCGA genotype files to a matrix.
    @input:
        `infiles:files`: The input genotypes files
    @output:
        `outfile:file`: The output matrix file
    @args:
        `rsmap`    : The rsid probe mapping file. If not provided, will use the probe id for matrix rownames.
        `fn2sample`: How to convert filename(without extension) to sample name.
        `confcut`  : The confidence cutoff. Genotype will be NA for lower confidence snps.
    """),
    input='infiles:files',
    output='outfile:file:{{i.infiles | fs2name}}.gt.txt',
    lang=opts.python,
    args=Diot(
        rsmap=opts.rsmap_gwsnp6,
        confcut=0.05,
        fn2sample="lambda fn: fn.split('.')[0]",
    ),
    envs=Diot(fs2name=fs2name)
)

pClinic2Survival = proc_factory(
    desc='Convert TCGA clinic data to survival data.',
    config=Diot(annotate="""
    @name:
        pClinic2Survival
    @description:
        Convert TCGA clinic data to survival data
        The clinic data should be downloaded as "BCR Biotab" format
    @input:
        `infile:file`: The clinic data file downloaded from TCGA
    @output:
        `outfile:file`: The output file
        `covfile:file`: The covariate file
    @args:
        `cols`: The column names:
            - `time_lastfollow`: The column names of last follow up. Default: `['days_to_last_followup']`
            - `time_death`: The column names of time to death. Default: `['days_to_death']`
            - `status`: The columns of vital status. Default: `['vital_status']`
            - `age`: The columns of days to birth. Default: `['days_to_birth']`
        `covs`: The covariates to output. Default:
            - `gender`, `race`, `ethnicity`, `age`
        `mat`: An expression or genotype matrix with samples as column names, used to get sample names for patient instead of short ones. Default: `None`
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | stem}}.survdata.txt, \
            covfile:file:{{i.infile | stem}}.survcov.txt',
    lang=opts.python,
    args=Diot(
        cols=Diot(time_lastfollow=['days_to_last_followup'],
                  time_death=['days_to_death'],
                  status=['vital_status'],
                  age=['days_to_birth'],
                  patient=['bcr_patient_barcode']),
        covs=['gender', 'race', 'ethnicity', 'age'],
        mat=None,
    )
)

pClinic2PlinkMeta = proc_factory(
    desc='Convert TCGA clinic data to plink meta file.',
    config=Diot(annotate="""
    @name:
        pClinic2PlinkMeta
    @description:
        Convert TCGA clinic data to plink meta file.
        Used with a `GTMat` file to convert to plink binary file (see vcfnext.pGTMat2Plink).
    @input:
        `infile:file`: The input TCGA patiant clinic file.
    @output:
        `outfile:file`: The meta file.
    @args:
        `suffix`: The suffix added to the barcode. Default: `''`
            - By default the barcode is like `TCGA-55-8087`, but in analysis, we want to be specific to samples
            - For example, for solid tumor, it'll be `TCGA-55-8087-01`, then suffix `-01` will be added
            - You can add multiple suffices by `-01, -10` or `['-01', '-10']`
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.meta.txt',
    lang=opts.python,
    args=Diot(suffix='')
)

pClinic2Cov = proc_factory(
    desc='Extract information from clinic data as covariance matrix',
    config=Diot(annotate="""
    @name:
        pClinic2Cov
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.cov.txt',
    lang=opts.Rscript,
    args=Diot(
        covs=[
            'gender', 'race', 'pathologic_stage',
            'age_at_initial_pathologic_diagnosis'
        ],
        asnum=True,
        suffix='',
    )
)
