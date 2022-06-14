from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.tcr import TCRClustering, TCRClusteringStats
from biopipen.core.testing import get_pipeline


class PrepareImmdata(Proc):
    """Prepare the data

    Requires:
        - name: immunarch
          check: |
            {{proc.lang}} <(echo "library(immunarch)")
    """

    input = "name"
    input_data = ["immdata"]
    output = "outfile:file:{{in.name}}.RDS"
    lang = config.lang.rscript
    script = """
        library(immunarch)
        data(immdata)
        saveRDS(immdata, {{out.outfile | quote}})
    """


class TCRClusteringGIANA(TCRClustering):
    requires = PrepareImmdata
    envs = {"tool": "GIANA"}


class TCRClusteringClusTCR(TCRClustering):
    requires = PrepareImmdata
    envs = {"tool": "ClusTCR"}


class TCRClusteringStatsGIANA(TCRClusteringStats):
    requires = TCRClusteringGIANA


class TCRClusteringStatsClusTCR(TCRClusteringStats):
    requires = TCRClusteringClusTCR


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareImmdata)


def testing(pipen):
    for proc in pipen.procs[-2:]:
        outfile = (
            proc.workdir.joinpath(
                "0",
                "output",
                "immdata.tcrclusters_stats",
                "SharedClusters",
                "shared_clusters.txt",
            )
        )
        assert outfile.is_file()


if __name__ == "__main__":
    pipen = pipeline()
    pipen.run()
    testing(pipen)
