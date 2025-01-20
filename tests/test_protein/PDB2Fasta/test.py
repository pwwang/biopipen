from biopipen.ns.web import Download
from biopipen.ns.protein import MMCIF2PDB as MMCIF2PDB_, PDB2Fasta as PDB2Fasta_
from biopipen.core.testing import get_pipeline


class DownloadCIF(Download):
    input_data = ["https://files.rcsb.org/download/1A2M.cif"]


class DownloadPDB(Download):
    input_data = ["https://files.rcsb.org/download/3NIR.pdb"]


class MMCIF2PDB(MMCIF2PDB_):
    requires = DownloadCIF
    input_data = lambda ch: [ch.iloc[0, 0]]


class PDB2Fasta1A2M(PDB2Fasta_):
    requires = MMCIF2PDB
    envs = {"chains": "B"}


class PDB2Fasta3NIR(PDB2Fasta_):
    requires = DownloadPDB


def pipeline():
    return get_pipeline(__file__).set_starts(DownloadCIF, DownloadPDB)


def testing(pipen):
    ...


if __name__ == "__main__":
    pipen = pipeline()
    assert pipen.run()
    testing(pipen)
