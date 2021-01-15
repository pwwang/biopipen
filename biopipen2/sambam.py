"""
A set of processes to generate/process sam/bam files
"""
from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pSam2Bam = proc_factory(
    desc=('Deal with mapped sam/bam files, including '
          'sort, markdup, rmdup, and/or index.'),
    config=Diot(
        annotate="""
            @input:
                infile: The input file
            @output:
                outfile: The output bam file
            @args:
                tool             : The tool used to do the sort. Available: sambamba|picard|sambamba|biobambam|samtools
                sambamba         : The path of the sambamba.
                picard           : The path of the picard.
                biobambam_bamsort: The path of the biobambam's bamsort.
                samtools: The path of the samtools.
                sort    : Do sorting?
                    - If input is sam, tool is biobambam, this should be True
                index  : Do indexing?
                markdup: Do duplicates marking?
                    - `rmdup` for samtools will be called
                rmdup  : Do duplicates removing?
                tmpdir : The tmp dir used to store tmp files.
                sortby : Sort by coordinate or queryname.
                nthread: Number of threads to use.
                infmt  : The format of input file. Available: sam|bam
                params : Other parameters for `tool`.
                mem    : The max memory to use.
                    - Unit could be G/g/M/m
                    - Will be converted to -Xmx4G, and -Xms will be 1/8 of it
            @requires:
                [sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) if `args.tool` == samtools or reference used but not indexed.
                [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
                [biobambam](https://github.com/gt1/biobambam2)
                [samtools](https://github.com/samtools/samtools)
        """,
        runcmd_pre="""
            {{"reference.bash" | bashimport}}
            export elprep={{args.elprep | quote}}
            if [[ {{args.tool | quote}} == "elprep" ]]; then
                reference elprep {{args.ref | squote}}
            fi
        """
    ),
    lang=opts.python,
    input="infile:file",
    output=
    "outfile:file:{{i.infile | fn}}.bam, outidx:file:{{i.infile | fn}}.bam.bai",
    errhow='retry',
    args=Diot(tool="biobambam",
              sambamba=opts.sambamba,
              picard=opts.picard,
              biobambam=opts.biobambam_bamsort,
              samtools=opts.samtools,
              elprep=opts.elprep,
              steps=Diot(sort=True,
                         index=True,
                         markdup=True,
                         rmdup=True,
                         recal=True),
              tmpdir=opts.tmpdir,
              sortby="coordinate",
              nthread=1,
              params=Diot(),
              mem=opts.mem16G,
              ref=opts.ref,
              knownSites='')
)

pBamMarkdup = proc_factory(
    desc='Mark/remove duplicates for bam files.',
    config=Diot(annotate="""
    @name:
        pBamMarkdup
    @description:
        Mark/remove duplicates for bam files
    @input:
        `infile:file`: The input file
    @output:
        `outfile:file`: The output bam file
    @args:
        `tool`             : The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools|bamutil)
        `sambamba`         : The path of sambamba. Default: sambamba
        `picard`           : The path of picard. Default: picard
        `biobambam_bamsort`: The path of biobambam's bamsort. Default: bamsort
        `samtools`         : The path of samtools. Default: samtools
        `bamutil`          : The path of bamutil. Default: bam
        `rmdup`            : Do duplicates removing? Default: False
            - Samtools will anyway remove the duplicates
        `tmpdir`           : The tmp dir used to store tmp files. Default: <system default tmpdir>
        `nthread`          : Default: 1
            - Not available for samtools and picard
        `params`           : Other parameters for `tool`. Defaut: ""
        `mem`              : The max memory to use. Default: "16G"
            - Unit could be G/g/M/m
            - Will be converted to -Xmx4G, and -Xms will be 1/8 of it
    @requires:
        [sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
        [biobambam](https://github.com/gt1/biobambam2)
        [samtools](https://github.com/samtools/samtools)
        [bamutil](http://genome.sph.umich.edu/wiki/BamUtil#Programs)
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | fn}}.bam",
    lang=opts.python,
    args=Diot(
        tool="biobambam",
        sambamba=opts.sambamba,
        picard=opts.picard,
        biobambam_bamsort=opts.biobambam_bamsort,
        samtools=opts.samtools,
        bamutil=opts.bamutil,
        rmdup=False,
        tmpdir=opts.tmpdir,
        nthread=1,
        params=Diot(),
        mem=opts.mem16G,
    )
)

pBamRecal = proc_factory(
    desc='Recalibrate a bam file.',
    config=Diot(
        annotate="""
            @name:
                pBamRecal
            @description:
                Recalibrate a bam file
            @input:
                `infile:file`: The bam file
            @output:
                `outfile:file`: The output bam file
            @args:
                `tool`    : The tool used to recalibrate the bam file. Default: `gatk` (gatk|bamutil)
                `gatk`    : The path of gatk, including java path. Default: `gatk`
                `samtools`: The path of samtools. Default: `samtools`
                `bamutil` : The path of bamutil. Default: `bam`
                `picard`  : The path of picard. Default: `picard`
                `params`  : Other parameters for `bam recab`. Default         : ""
                    `RealignerTargetCreator` : Other parameters for `gatk RealignerTargetCreator`. Defaut: ""
                    `IndelRealigner`         : Other parameters for `gatk IndelRealigner`. Defaut: ""
                    `BaseRecalibrator`       : Other parameters for `gatk BaseRecalibrator`. Defaut: ""
                    `PrintReads`             : Other parameters for `gatk PrintReads`. Defaut: ""
                `mem`: The max memory to use. Default: "32G"
                `knownSites`: The known polymorphic sites to mask out. Default: "" (Required for GATK)
                `ref`: The reference file. Required.
                    - Will be converted to -Xmx4G, and -Xms will be 1/8 of it
            @requires:
                [gatk](https://software.broadinstitute.org/gatk)
                [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed, or bamutil is used for bam index file generation.
                [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
            """,
        runcmd_pre="""
            {{"reference.bash" | bashimport}}
            export samtools={{args.samtools | squote}}
            reference fasta {{args.ref | squote}}
        """
    ),
    input="infile:file",
    output="outfile:file:{{i.infile | fn}}.bam, \
            outidx:file:{{i.infile | fn}}.bam.bai",
    lang=opts.python,
    args=Diot(
        tool="bamutil",
        gatk=opts.gatk,
        samtools=opts.samtools,
        picard=opts.picard,
        bamutil=opts.bamutil,
        params=Diot(
            # For GATK
            # RealignerTargetCreator = Diot(),
            # IndelRealigner         = Diot(),
            # BaseRecalibrator       = Diot(),
            # PrintReads             = Diot()
        ),
        ref=opts.ref,
        tmpdir=opts.tmpdir,
        knownSites=opts.dbsnp,
        nthread=1,
        mem=opts.mem32G,
    )
)

pBamReadGroup = proc_factory(
    desc='Add or replace read groups of a bam file.',
    config=Diot(annotate="""
    @name:
        pBamReadGroup
    @description:
        Add or replace read groups of a bam file
    @input:
        `infile:file`: The bam file
    @output:
        `outfile:file`: The output bam file
    @args:
        `tool`                         : The tool used. Default: `picard` (picard|bamutil)
        `picard`                       : The path of picard. Default: `picard`
        `bamutil`                      : The path of bamutil. Default: `bam`
        `rg`                           : The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
            - `id` will be parsed from filename with "_LX_" in it if not given
            - `sm` will be parsed from filename
        `params`                       : Other parameters for `tool`. Defaut: ""
        `mem`                          : The max memory to use. Default: "4G"
            - Will be converted to -Xmx4G, and -Xms will be 1/8 of it
        `tmpdir`                       : The temporary directory. Default: <system tmpdir>
    @requires:
        [gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
        [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | bn}}",
    lang=opts.python,
    args=Diot(
        tool="bamutil",
        picard=opts.picard,
        bamutil=opts.bamutil,
        rg=Diot({
            'id': '',
            'pl': 'Illumina',
            'pu': 'unit1',
            'lb': 'lib1',
            'sm': ''
        }),
        params=Diot(),
        tmpdir=opts.tmpdir,
        mem=opts.mem4G,
    )
)

pBamReorder = proc_factory(
    desc='Reorder a sam/bam file by a given reference.',
    config=Diot(annotate="""
    @name:
        pBamReorder
    @description:
        Reorder a sam/bam file by a given reference file using `picard ReorderSam`
    @input:
        `infile:file`: The sam/bam file
    @output:
        `outfile:file`: The output bam file
    @args:
        `picard`                       : The path of picard. Default: `picard`
        `ref`                          : The reference file. Required
        `params`                       : Other parameters for `picard ReorderSam`. Defaut: ""
        `mem`                          : The max memory to use. Default: "4G"
            - Will be converted to -Xmx4G, and -Xms will be 1/8 of it
        `tmpdir`                       : The temporary directory. Default: <system tmpdir>
    @requires:
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
    """),
    input="infile:file",
    output="outfile:file:{{i.infile | bn}}",
    args=Diot(
        picard=opts.picard,
        params=Diot(),
        tmpdir=opts.tmpdir,
        mem=opts.mem4G,
        ref=opts.ref,
    )
)

pBamMerge = proc_factory(
    desc='Merges multiple SAM and/or BAM sorted files into a single file.',
    config=Diot(annotate="""
    @name:
        pBamMerge
    @description:
        Merges multiple SAM and/or BAM files (must be sorted by coordinate) into a single file.
    @input:
        `infiles:files`: Input sam/bam files to be merged
    @output:
        `outfile:file`: The merged bam file
    @args:
        `tool`     : The tool used to merge. Default: bamutil (picard|samtools|sambamba)
        `picard`   : The path of picard. Default: `picard`
        `bamutil`  : The path of bamutil. Default: `bam`
        `samtools` : The path of samtools. Default: `samtools`
        `sambamba` : The path of sambamba. Default: `sambamba`
        `params`   : Other parameters for `tool`. Defaut: ""
        `mem`      : The max memory to use. Default: "4G"
            - Will be converted to -Xmx4G, and -Xms will be 1/8 of it, just for picard
        `tmpdir`   : The temporary directory. Default: <system tmpdir>
        `nthread`  : # threads to use. Default: 1
            - For picard, if nthread>1, USE_THREADING=true, otherwise USE_THREADING=false
    @requires:
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
    """),
    input="infiles:files",
    output="outfile:file:{{i.infiles[0] | fn2}}.etc.bam",
    lang=opts.python,
    args=Diot(
        tool="picard",
        picard=opts.picard,
        bamutil=opts.bamutil,
        samtools=opts.samtools,
        sambamba=opts.sambamba,
        params=Diot(),
        tmpdir=opts.tmpdir,
        nthread=1,
        mem=opts.mem4G,
    )
)

pBam2Gmut = proc_factory(
    desc='Call germline (snps and indels) from a call-ready bam file.',
    lang=opts.python,
    config=Diot(
        annotate="""
            @name:
                pBam2Gmut
            @description:
                Call germline (snps and indels) from a call-ready bam file.
            @input:
                `infile:file`: The input bam file
            @output:
                `outfile:file`: The vcf file containing the mutations
            @args:
                `tool`      : The tool used to call mutations. Default: gatk (vardict, snvsniffer, platypus, strelka)
                `gatk`      : The path of gatk. Default: gatk
                `vardict`   : The path of vardict. Default: vardict
                `snvsniffer`: The path of snvsniffer. Default: SNVSniffer
                `samtools`  : The path of samtools. Default: samtools (used to generate reference index)
                `platypus`  : The path of platypus. Default: platypus
                `strelka`   : The path of strelka. Default: configureStrelkaGermlineWorkflow.py
                `cfgParams` : The params for `strelka` configuration. Default: ""
                `picard`    : The path of picard. Default: picard
                `mem`       : The memory to be used. Default: 32G
                    - will be converted to -Xms4G -Xmx32G for java programs
                `ref`:          The reference file. Required.
                `gz`:           Gzip output file? Default: False
                `tmpdir`:       The temporary directory. Default: <system tmpdir>
                `params`:       Other params for `tool`. Default: ""
            @requires:
                [gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
                [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
                [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
                [vardict](https://github.com/AstraZeneca-NGS/VarDict)
                [snvsniffer](http://snvsniffer.sourceforge.net/homepage.htm#latest)
                [platypus](http://www.well.ox.ac.uk/platypus)
                [strelka@2.7.1+](https://github.com/Illumina/strelka)
            """,
        runcmd_pre="""
            {{"reference.bash" | bashimport}}
            export samtools={{args.samtools | squote}}
            export picard={{args.picard | squote}}
            reference fasta {{args.ref | squote}}
            reference picard {{args.ref | squote}}
        """
    ),
    input="infile:file",
    output="outfile:file:{{i.infile | fn}}.vcf{{args.gz|?|=:'.gz'|!:''}}",
    args=Diot(
        tool="strelka",
        gatk=opts.gatk,
        vardict=opts.vardict,
        snvsniffer=opts.snvsniffer,
        # required by SNVSniffer to generate a bam header file,
        samtools=opts.samtools,
        platypus=opts.platypus,
        picard=opts.picard,
        strelka=opts.strelka_germ,
        mem=opts.mem24G,
        ref=opts.ref,
        tmpdir=opts.tmpdir,
        cfgParams=Diot(),  # only for strelka,
        params=Diot(),
        gz=False,
        nthread=1,  # for gatk and platypus,
    )
)


pBamPair2Smut = proc_factory(
    desc='Call somatic mutations from tumor-normal bam pair.',
    config=Diot(
        annotate="""
            @name:
                pBamPair2Smut
            @description:
                Call somatic mutations from tumor-normal bam pair.
            @input:
                `tumor:file`: The tumor bam file
                `normal:file`: The normal bam file
            @output:
                `outfile:file`: The vcf file
            @args:
                `tool`: The tool used to call mutations. Default: gatk (somaticsniper, strelka, snvsniffer, virmid, varidct)
                `gatk`: The path to gatk. Default: gatk
                `somaticsniper`: The path to gatk. Default: bam-somaticsniper
                `strelka`: The path to gatk. Default: configureStrelkaSomaticWorkflow.py
                `snvsniffer`: The path to gatk. Default: SNVSniffer
                `virmid`: The path to gatk. Default: virmid
                `vardict`: The path to gatk. Default: vardict
                `samtools`: The path to gatk. Default: samtools
                `picard`: The path to gatk. Default: picard
                `configParams`: The configuration parameters for `configureStrelkaSomaticWorkflow.py`. Default: `{}`
                `params`: The parameters for main programs. Default: `{}`
                `meme`: The memory. Default: 24G
                `ref`: The reference genom. Default: `opts.ref`
                `gz`: Whether gzip the output vcf file. Default: False
                `nthread`: The number of threads to use. Default: 1
                `tmpdir`: The temporary directory. Default: `opts.tmpdir`
            @requires:
                [gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
                [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
                [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
                [vardict](https://github.com/AstraZeneca-NGS/VarDict)
                [snvsniffer](http://snvsniffer.sourceforge.net/homepage.htm#latest)
                [platypus](http://www.well.ox.ac.uk/platypus)
                [strelka@2.7.1+](https://github.com/Illumina/strelka)
            """,
        runcmd_pre="""
            {{"reference.bash" | bashimport}}
            export samtools={{args.samtools | squote}}
            export picard={{args.picard | squote}}
            reference fasta {{args.ref | squote}}
            reference picard {{args.ref | squote}}
        """
    ),
    lang=opts.python,
    input="tumor:file, normal:file",
    output="outfile:file:{{i.tumor | fn}}-{{ \
        i.normal | fn}}.vcf{{args.gz|?|=:'.gz'|!:''}}",
    args=Diot(
        tool='strelka',
        gatk=opts.gatk,  # required for strelka,
        somaticsniper=opts.somaticsniper,
        strelka=opts.strelka_soma,  # @2.7.1,
        snvsniffer=opts.snvsniffer,
        virmid=opts.virmid,
        vardict=opts.vardict,
        samtools=opts.samtools,
        picard=opts.picard,
        configParams=Diot(),  # only for strelka,
        params=Diot(),
        mem=opts.mem24G,
        ref=opts.ref,
        gz=False,
        nthread=1,
        tmpdir=opts.tmpdir,
    )
)

pBamStats = proc_factory(
    desc='Get read depth from bam files.',
    config=Diot(
        annotate="""
            @name:
                pBamStats
            @description:
                Get read depth from bam files.
            @input:
                `infile:file`: The input bam file
            @output:
                `outfile:file`: The output statistic file
                `outdir:dir`:   The directory containing result files and figures.
            @args:
                `tool`: The tool used to do the job. Default: bamstats
                `bamstats`: The path to bamstats. Default: bamstats
                `params`: Other params to main program. Default: `{}`
                `mem`: The memory to be used. Default: 16G
                `plot`: Whether plot the result. Default: True
            """,
        runcmd_pre="""
            if [[ "{{args.plot | R}}" == "TRUE" && {{proc.forks}} -lt 2 ]]; then
                echo "Plots can only be done with proc.forks >= 2." 1>&2
                exit 1
            fi
        """
    ),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}/{{i.infile | fn}}.stat.txt, \
            outdir:dir:{{i.infile | fn}}',
    lang=opts.Rscript,
    args=Diot(
        tool='bamstats',
        bamstats=opts.bamstats,
        params=Diot(),
        mem=opts.mem16G,
        plot=True,
        histplotggs=Diot(xlab={0: 'Read depth'},
                         ylab={0: '# samples'}),
        boxplotggs=Diot(xlab={0: 'Counts'}),
        devpars=Diot(res=300, width=2000, height=2000),
        cap=500,
        cutoff=0,
        nfeats=40,
        feature='wgs',
    )
)

pBam2Fastq = proc_factory(
    desc='Convert sam/bam files to pair-end fastq files.',
    config=Diot(annotate="""
    @input:
        infile: The sam/bam file.
            - Sam files only available for biobambam, picard
    @output:
        fqfile1: The 1st match of paired reads
        fqfile2: The 2nd match of paired reads
    @args:
        `tool`     : The tool to use. Default: biobambam (bedtools, samtools, picard)
        `biobambam`: The path of bamtofastq of biobambam. Default: bamtofastq
        `bedtools` : The path of bedtools. Default: bedtools
        `samtools` : The path of samtools. Default: samtools
        `picard`   : The path of picard. Default: picard
        `mem`      : The memory to be used by picard. Default: 8G
        `gz`       : Whether gzip the output files. Default: True
        `params`:  : Other params for `tool`. Default: ''
        `tmpdir`   : The tmpdir. Default: `__import__('tempfile').gettempdir()`
    @requires:
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
        [biobambam](https://github.com/gt1/biobambam2)
        [samtools](https://github.com/samtools/samtools)
        [bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
    """),
    input="infile:file",
    output=[
        "fqfile1:file:{{ i.infile | fn }}_1.fastq{{args.gz|?|=:'.gz'|!:''}}",
        "fqfile2:file:{{ i.infile | fn }}_2.fastq{{args.gz|?|=:'.gz'|!:''}}"
    ],
    lang=opts.python,
    args=Diot(
        tool='biobambam',
        biobambam=opts.biobambam_bamtofastq,
        bedtools=opts.bedtools,
        samtools=opts.samtools,
        picard=opts.picard,
        mem=opts.mem8G,  # only for picard
        gz=False,
        params=Diot(),
        tmpdir=opts.tmpdir)
)

pBam2FastqSE = proc_factory(
    desc='Convert bam files to single-end fastq files.',
    config=Diot(annotate="""
    @name:
        pBam2FastqSE
    @description:
        Convert sam/bam files to single-end fastq files.
    @input:
        `infile:file`: The sam/bam file.
            - Sam files only available for biobambam, picard
    @output:
        `fqfile:file`: The fastq file
    @args:
        `tool`     : The tool to use. Default: biobambam (bedtools, samtools, picard)
        `biobambam`: The path of bamtofastq of biobambam. Default: bamtofastq
        `bedtools` : The path of bedtools. Default: bedtools
        `samtools` : The path of samtools. Default: samtools
        `picard`   : The path of picard. Default: picard
        `mem`      : The memory to be used by picard. Default: 8G
        `gz`       : Whether gzip the output files. Default: True
        `params`:  : Other params for `tool`. Default: ''
        `tmpdir`   : The tmpdir. Default: `__import__('tempfile').gettempdir()`
    @requires:
        [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
        [biobambam](https://github.com/gt1/biobambam2)
        [samtools](https://github.com/samtools/samtools)
        [bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
    """),
    input="infile:file",
    output="fqfile:file:{{i.infile | fn }}.fq{% if args.gz %}.gz{% endif %}",
    lang=opts.python,
    args=Diot(
        tool='biobambam',
        biobambam=opts.biobambam_bamtofastq,
        bedtools=opts.bedtools,
        samtools=opts.samtools,
        picard=opts.picard,
        mem=opts.mem8G,  # only for picard,
        gz=False,
        params=Diot(),
        tmpdir=opts.tmpdir,
    )
)

pBam2Counts = proc_factory(
    desc='Extract read counts from RNA-seq bam files.',
    config=Diot(annotate="""
    @name:
        pBam2Counts
    @description:
        Extract read counts from RNA-seq bam files.
    @input:
        `infile:file`: The input bam files
    @outfile:
        `outfile:file`: The count file
    @args:
        `tool`: The tool used to extract counts. Default: ht-seq
        `htseq`: The path of htseq-count.
        `params`: Other params for main program.
        `refgene`: The reference gene in GTF format.
    @requires:
        [`htseq`](https://htseq.readthedocs.io/)
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | fn}}.counts',
    lang=opts.python,
    args=Diot(
        tool='htseq',
        htseq=opts.htseq_count,
        params=Diot(),
        refgene=opts.refexon,
    )
)

pBamIndex = proc_factory(
    desc='Index bam files',
    config=Diot(annotate="""
    @name:
        pBamIndex
    @description:
        Index bam files.
    @input:
        `infile:file`: The input bam file
    @output:
        `outfile:file`: The symbolic link to the input file
        `outidx:file` : The index file
    @args:
        `samtools`: Path to samtools. Default: `params.samtools`
        `params`  : Other parameters for samtools. Default: `Diot(b = True)`
        `nthread` : # threads to use. Default: `1`
    """),
    input='infile:file',
    output='outfile:file:{{i.infile | bn}}, outidx:file:{{i.infile | bn}}.bai',
    lang=opts.python,
    args=Diot(
        tool='samtools',
        samtools=opts.samtools,
        params=Diot(b=True),
        nthread=1,
    )
)
