"""ChIP-seq data analysis"""

from diot import Diot
from . import opts, proc_factory

# pylint: disable=invalid-name

pPeakToRegPotential = proc_factory(
    desc='Convert peaks to regulatory potential score for each gene.',
    config=Diot(annotate="""
    @name:
        pPeakToRegPotential
    @description:
        Convert peaks to regulatory potential score for each gene
        The formula is:
        ```
                        -(0.5 + 4*di/d0)
        PC = sum (pi * e                  )
        ```
        Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4489297/
    @input:
        `peakfile:file`: The BED/peak file for peaks
        `genefile:file`: The BED file for gene coordinates
    @output:
        `outfile:file`: The regulatory potential file for each gene
    @args:
        `signal`: `pi` in the formula. Boolean value, whether use the peak intensity signale or not, default: `True`,
        `genefmt`: The format for `genefile`, default: `ucsc+gz`. It could be:
            - ucsc or ucsc+gz: typically, you can download from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
            - bed or bed+gz: [format](https://genome.ucsc.edu/FAQ/FAQformat#format1), 4th column required as gene identity.
        `peakfmt`: The format for `peakfile`, default: `peak`. It could be:
            - peak or peak+gz: (either [narrowPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) or [broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13), the 7th column will be used as intensity
            - bed or bed+gz: [format](https://genome.ucsc.edu/FAQ/FAQformat#format1), 5th column will be used as intensity.
        `window`: `2 * d0` in the formula. The window where the peaks fall in will be consided, default: `100000`.
            ```
                    |--------- window ----------|
                    |---- d0 -----|
                    |--- 50K --- TSS --- 50K ---|
                        ^ (peak center)
                        |-- di --|
            ```
    """),
    input="peakfile:file, genefile:file",
    output="outfile:file:{{peakfile | fn}}.rp.txt",
    lang=opts.python,
    args=Diot(
        signal=True,
        genefmt='ucsc+gz',
        peakfmt='peak',
        window=100000,
    )
)
