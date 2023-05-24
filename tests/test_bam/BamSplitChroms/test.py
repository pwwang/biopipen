from pathlib import Path
from biopipen.ns.web import Download
from biopipen.ns.bam import BamMerge, BamSplitChroms
from biopipen.core.testing import get_pipeline


# TOOL = "samtools"
TOOL = "sambamba"
BAM_URL = (
    "https://github.com/VCCRI/SVPV/raw/master/example/NA12877_S1.partial.bam"
)


BamSplitChroms.requires = Download
BamSplitChroms.input_data = lambda ch: [ch.iloc[0, 0]] * 2
BamSplitChroms.envs.update({"ncores": 2, "tool": TOOL})


def bammerge_input_data(ch):
    # ch is a list of directories with bam files split by chromosomes
    # we try to merge them by chromosomes
    # get chromosomes
    outdirs = [Path(p) for p in ch.iloc[:, 0]]
    bamfiles = outdirs[0].glob("*.bam")
    chroms = [bamfile.stem for bamfile in bamfiles]
    out = []
    for chrom in chroms:
        if chrom not in ["chr1", "chr2"]:
            # only merge chr1 and chr2
            continue
        out.append(
            [str(outdir.joinpath(f"{chrom}.bam")) for outdir in outdirs]
        )
    return out


BamMerge.requires = BamSplitChroms
BamMerge.input_data = bammerge_input_data
BamMerge.envs.update({"ncores": 2, "tool": TOOL})


def pipeline():
    return (
        get_pipeline(__file__, plugins=["no:report"])
        .set_start(Download)
        .set_data([BAM_URL])
    )


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
