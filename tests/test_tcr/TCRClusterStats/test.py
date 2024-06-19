from pipen import Proc
from biopipen.core.config import config
from biopipen.ns.tcr import TCRClustering, TCRClusterStats
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
        set.seed(8525)
        for (name in names(immdata$data)) {
            immdata$data[[name]]$Barcode = paste0(
                "Cell-", 1:nrow(immdata$data[[name]])
            )
        }
        saveRDS(immdata, {{out.outfile | quote}})
    """


class TCRClusteringGIANA(TCRClustering):
    requires = PrepareImmdata
    envs = {"tool": "GIANA"}


class TCRClusteringClusTCR(TCRClustering):
    requires = PrepareImmdata
    envs = {"tool": "ClusTCR"}


class TCRClusterStatsGIANA(TCRClusterStats):
    requires = TCRClusteringGIANA
    order = 98


class TCRClusterStatsClusTCR(TCRClusterStats):
    requires = TCRClusteringClusTCR
    order = 99
    envs = {
        "shared_clusters": {
            "cluster_rows": False,
            "heatmap_meta": "Sex",
            "sample_order": ["MS3", "MS5", "MS1", "MS2"],
        }
    }


def pipeline():
    return get_pipeline(__file__).set_starts(PrepareImmdata)


def testing(pipen):
    # assert pipen._succeeded
    for proc in pipen.procs[-2:]:
        outfile = (
            proc.workdir.joinpath(
                "0",
                "output",
                "immdata.tcrclusters_stats",
                "SharedClusters",
                "DEFAULT",
                "shared_clusters.txt",
            )
        )
        assert outfile.is_file(), outfile


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
