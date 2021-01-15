"""A set of immunotherapy-related bioinformatics tools"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pTopiary = proc_factory(
    desc='Epitope prediction using Topiary',
    config=Diot(annotate="""
    @description:
        Predict mutation-derived cancer T-cell epitopes from somatic variants, tumor RNA expression data, and patient HLA type.
    @input:
        infile: The input VCF file.
            - Can be annotated with vcf_expression_annotator, as well as VEP
        afile: The HLA alleles in format of "HLA-A*02:01,HLA-B*07:02" or a file with one per line.
            - You can have star in them: "HLA-A*02:01,HLA-B*07:02"
    @output:
        outfile: The output file.
            - A html file will also be generated with suffix '.html'
        outdir: The output directory with output file and all intermediate files
    @args:
        topiary       (path): Path to topiary.
        netmhc        (path): Path to netmhc.
        netmhcpan     (path): Path to netmhcpan.
        netmhciipan   (path): Path to netmhciipan.
        netmhccons    (path): Path to netmhccons.
        smm           (path): Path to smm.
        smm_pmbec     (path): Path to smm-pmbec.
        mhc_predictor (path): The MHC binding predictor. Could be one of:
            - netmhc: Local NetMHC predictor (Topiary will attempt to automatically detect whether NetMHC 3.x or 4.0 is available)
            - netmhcpan: Local NetMHCpan predictor
            - netmhciipan: Local NetMHCIIpan predictor
            - netmhccons: Local NetMHCcons
            - random: Random IC50 values
            - smm: Local SMM predictor
            - smm-pmbec: Local SMM-PMBEC predictor
            - netmhcpan-iedb: Use NetMHCpan via the IEDB web API
            - netmhccons-iedb: Use NetMHCcons via the IEDB web API
            - smm-iedb: Use SMM via the IEDB web API
            - smm-pmbec-iedb: Use SMM-PMBEC via the IEDB web API
        params (Diot): Other parameters for topiary
        tmpdir (str): Temporary directory for running local MHC predictors
    """),
    lang=opts.python,
    input='infile:file, afile:var',
    output=[
        'outfile:file:{{i.infile | stem2}}.topiary/{{i.infile | stem2}}.txt',
        'outdir:dir:{{i.infile | stem2}}.topiary'
    ],
    args=Diot(topiary=opts.topiary,
              netmhc=opts.netmhc,
              netmhcpan=opts.netmhcpan,
              netmhciipan=opts.netmhciipan,
              netmhccons=opts.netmhccons,
              genome=opts.genome,
              smm=opts.smm,
              smm_pmbec=opts.smm_pmbec,
              tmpdir=opts.tmpdir,
              refall=opts.refall,
              mhc_predictor='netmhcpan',
              params=Diot()))

pVACseq = proc_factory(
    desc=('Identification of Variant Antigens using tumor '
          'mutation and expression data.'),
    config=Diot(annotate="""
    @description:
        pVACseq is a cancer immunotherapy pipeline for the identification of personalized Variant Antigens by Cancer Sequencing (pVACseq) that integrates tumor mutation and expression data (DNA- and RNA-Seq). It enables cancer immunotherapy research by using massively parallel sequence data to predicting tumor-specific mutant peptides (neoantigens) that can elicit anti-tumor T cell immunity. It is being used in studies of checkpoint therapy response and to identify targets for personalized cancer vaccines and adoptive T cell therapies. For more general information.
        See: http://www.genomemedicine.com/content/8/1/11
    """),
    lang=opts.python,
    input='infile:file, afile:var',
    output=[
        'outfile:file:{{i.infile | stem}}.pvacseq'
        '/{{i.infile | stem}}.epitopes.txt',
        'outdir:dir:{{i.infile | stem}}.pvacseq'
    ],
    args=Diot(
        # used to extract sample name
        bcftools=opts.bcftools,
        pvacseq=opts.pvacseq,
        allele='',
        bdtool='netmhc',
        netmhc=opts.netmhc,
        iedb_mhc_i=opts.iedb_mhc_i,
        nthread=1,
        params=Diot(epitope_length=9)
    )
)

pOptiType = proc_factory(
    desc=
    'Precision HLA typing from next-generation sequencing data using OptiType',
    config=Diot(annotate="""
    @description:
        Precision HLA typing from next-generation sequencing data using OptiType
        See: https://github.com/FRED-2/OptiType
    @input:
        fqfile1: The bam file or the fastq file of the 1st pair-end
        fqfile2: The fastq file of the 2nd pair-end
            - We will try to extract hla sequences from the bam/fastq file.
    @output:
        outfile: The output file with HLA types
        outdir: The output directory
    @args:
        optitype (str): Path to OptiTypePipeline.py
        picard   (str): Path to picard, used to convert bam files to fastq files.
            - Samtools not used, as it may mark paired reads with different names, which causes bwa crash.
        bwa      (str): Path to bwa
        hlaref   (str): Path to HLA reference file
        params   (Diot): Other parameters for optitype
        nthread  (int): Number of threads to use
    @requires:
        [OptiType](https://github.com/FRED-2/OptiType)
        bwa
        picard
    """),
    lang=opts.python,
    input='fqfile1:file, fqfile2:file',
    output=[
        'outfile:file:{% 	assign outdir = i.fqfile1, i.fqfile2 | \
                            ? :bool(_[1]) | \
                            = :__import__("os").path.commonprefix | \
                            ! :_[0] | stem2 | .rstrip: "._[]" | \
                            @append: ".optitype" %}{{outdir}}/{{outdir}}.txt',
        'outdir:dir:{{		i.fqfile1, i.fqfile2 | \
                            ? :bool(_[1]) | \
                            = :__import__("os").path.commonprefix | \
                            ! :_[0] | stem2 | .rstrip: "._[]"}}.optitype'
    ],
    args=Diot(optitype=opts.optitype,
              picard=opts.picard,
              bwa=opts.bwa,
              hlaref=opts.hlaref,
              nthread=opts.nthread,
              params=Diot(dna=True, bwa=Diot()))
)

pHLA_LA = proc_factory(
    desc='HLA typing using HLA-LA',
    config=Diot(annotate="""
    @input:
        infile: The bam file
    @output:
        outfile: The output file with HLA types
        outdir: The output directory
    @args:
        hla_la   (path): Path to HLA-LA.pl
        picard   (path): Path to picard
        bwa      (path): Path to bwa
        samtools (path): Path to samtools
        nthread  (int) : Number of threads to use
        params   (Diot) : Other parameters for `HLA-LA.pl`
    @requires:
        [HLA-LA](https://github.com/DiltheyLab/HLA-LA)
    """),
    lang=opts.python,
    input='infile:file',
    output=[
        'outfile:file:{{i.infile | stem2}}/{{i.infile | stem2}}.hlala.txt',
        'outdir:dir:{{i.infile | stem2}}'
    ],
    args=Diot(
        hla_la=opts.hla_la,
        picard=opts.picard,
        bwa=opts.bwa,
        samtools=opts.samtools,
        java=opts.java,
        params=Diot(),
        nthread=1,
    )
)

pVcf2ProteinSeq = proc_factory(
    desc='Generate protein sequences from VCF file using pVACseq',
    config=Diot(annotate="""
    @input:
        infile: The input VCF file
    @output:
        outdir: The output directory containing the protein sequences
    """),
    lang=opts.python,
    input='infile:file',
    output='outdir:dir:{{i.infile | stem2}}.protseqs',
    args=Diot(
        pvacseq=opts.pvacseq,
        lens=[15, 17, 19, 21],
        params=Diot(),
        nthread=1,
    )
)

pNetMHC = proc_factory(
    desc=('Run the netmhc to predict binding affinity '
          'between HLA-allele and peptides'),
    config=Diot(annotate="""
    @input:
        infile: Peptide sequences in fasta format or list of peptide sequences, optionally with 2nd columns of numeric values.
        afile: A list (separated by comma) or a file of alleles, one per line. For example:
            - HLA-A*03:79,HLA-A*04:03
    @output:
        outfile: The output file of binding affinity of all combinations (Peptide ~ HLA allele)
    @args:
        netmhc  (path)        : Path to netMHC
        isfa    (bool)        : Whether the input file is fasta sequence file (otherwise it is peptides)
        params  (Diot)         : Other parameters for netMHC
        lens    (int|str|list): Peptide length
        tmpdir  (path)        : Temporary directory for netMHC
        nthread (int)         : Number of threads to use.
            - netMHC itself does not multi-threading.
            - We will split the input file and put them into multiple threads
    """),
    input='infile:file, afile:var',
    output='outfile:file:{{i.infile | stem2}}.netmhc.xls',
    lang=opts.python,
    args=Diot(netmhc=opts.netmhc,
              isfa=True,
              params=Diot(v=True),
              lens=9,
              tmpdir=opts.tmpdir,
              nthread=1)
)
