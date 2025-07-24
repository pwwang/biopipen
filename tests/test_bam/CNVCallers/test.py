from biopipen.ns.web import Download
from biopipen.core.proc import Proc
from biopipen.ns.bam import (
    CNVpytor as CNVpytor_,
    ControlFREEC as ControlFREEC_,
    CNAClinic as CNAClinic_,
)
from biopipen.core.testing import get_pipeline


BAM_URL = [
    "https://raw.githubusercontent.com/etal/cnvkit/refs/heads/master/test/formats/na12878-chrM-Y-trunc.bam",  # noqa: E501
    # "https://raw.githubusercontent.com/VCCRI/SVPV/refs/heads/master/example/NA12878_S1.partial.bam",  # noqa: E501
    # "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam",  # noqa
]


class CNVpytor(CNVpytor_):
    requires = Download
    envs = {
        "binsizes": [1000],
        "filters": {
            # # https://github.com/abyzovlab/CNVpytor/blob/master/
            # # GettingStarted.md#predicting-cnv-regions
            # "CNVsize": [1, "inf"],  # CNV size
            # "eval1": [0, 0.9],  # filter on e-val1
            # "eval2": [],  # filter on e-val2
            # "eval3": [],  # filter on e-val3
            # "eval4": [],  # filter on e-val4
            # "q0": [-1, 0.5],  # filter on Q0
            # "pN": [0, 0.5],  # filter on pN
            # "dG": [1, "inf"],  # filter on dG
        },
    }


class ControlFREEC(ControlFREEC_):
    requires = Download


class CNAClinicMeta(Proc):
    """Metadata for CNAClinic."""
    requires = Download
    input = "bamfile:file"
    output = "outfile:file:cnaclinic_meta.txt"
    script = """
        echo "Bam" > {{out.outfile}}
        echo "{{in.bamfile}}" >> {{out.outfile}}
    """


class CNAClinic(CNAClinic_):
    requires = CNAClinicMeta
    envs = {"binsize": 10000}


def pipeline():
    return (
        get_pipeline(__file__)
        # get_pipeline(__file__, enable_report=True)
        .set_start(Download)
        .set_data(BAM_URL)
    )


def testing(pipen):
    # assert pipen._succeeded
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
